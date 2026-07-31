import logging
import os

from rdkit.Chem.rdchem import MolSanitizeException

from ..generation import _telemetry
from . import clash, clean_geometry, embed, om
from .utils.compute_chg_and_bo_pulp import clear_pulp_cache

logger = logging.getLogger(__name__)

# Structural, non-transient embed failures: the assembled complex is chemically
# invalid or its geometry cannot be built, so *every* attempt fails identically
# regardless of seed or scale. Distinct from a transient per-attempt failure
# (a single bad embed), which is retried. Confirmed dominant case on real
# bucket-B molecules: a dative donor drawn covalently over-valences its atom
# (e.g. a tertiary amine N -> metal makes N 4-valent), raising RDKit's
# AtomValenceException (a MolSanitizeException). TypeError/IndexError cover the
# originally-hypothesised under-coordination path (a None binding_site or an
# empty direction_vector in embed's cmap construction).
_STRUCTURAL_EMBED_ERRORS = (MolSanitizeException, TypeError, IndexError)


class StructuralAssemblyError(ValueError):
    """Every embed attempt failed with the same structural (non-transient) error.

    Raised when the attempt pool exhausted and no conformer was produced because
    the complex could not be assembled at all -- most often an over-valent dative
    donor (an ``AtomValenceException`` from RDKit sanitization), or the
    originally-hypothesised under-coordination (a ``TypeError``/``IndexError``
    while filling the coordination geometry). These repeat for every seed and
    scale, so retrying cannot help; surfacing them by *reason* (with the
    underlying cause chained) lets the harness report *why* instead of a generic
    ``no_conformers``.

    Subclasses ``ValueError`` so existing callers that catch ``ValueError`` keep
    working (mirrors ``UncoordinatedFragmentError`` in ``metallogen_adapter``).
    """


class BudgetExhaustedError(ValueError):
    """The generation budget ran out before any conformer was produced.

    Raised only when ``OIN_ENFORCE_BUDGET`` is on. It exists so a downstream
    consumer -- the sweep harness, the next release's regression triage -- can tell
    **"we stopped because we ran out of time"** from **"this molecule cannot be
    assembled"**. Without the distinction v0.4.10 cannot separate its own regressions
    from v0.4.9's intended behaviour: both would arrive as the same generic
    ``MetalloGen failed to generate any conformers``, which
    ``tools/classify_failures.py`` buckets as ``no_conformers``.

    The message carries the budget **and** the wall-clock actually spent, because a
    bound is only interpretable next to the number it was given: ``ULODUU`` assembles
    under a 60 s cap and not under 30 s, which is exactly why the blanket boron
    fast-fail was refuted.

    Raised for an **empty pool only.** A budget that expires with usable conformers in
    hand returns them -- a bound should stop work, not discard an answer it already
    has. ``StructuralAssemblyError`` takes priority over it: a uniformly-structural
    failure would have happened at any budget, so it is the better diagnosis.

    Subclasses ``ValueError``, like its two siblings.
    """


def calculate_heavy_atom_rmsd(mol1, mol2):
    """Calculate the heavy atom rmsd."""
    import numpy as np
    from scipy.spatial.transform import Rotation

    c1 = np.array([a.get_coordinate() for a in mol1.atom_list if a.get_atomic_number() > 1])
    c2 = np.array([a.get_coordinate() for a in mol2.atom_list if a.get_atomic_number() > 1])
    if len(c1) == 0 or len(c1) != len(c2):
        return float("inf")
    c1 -= c1.mean(axis=0)
    c2 -= c2.mean(axis=0)
    try:
        rot, rmsd = Rotation.align_vectors(c1, c2)
        return rmsd
    except Exception:
        return float("inf")


def _dihedral_deg(p0, p1, p2, p3):
    """Signed dihedral (degrees) of the p0-p1-p2-p3 sequence."""
    import math

    import numpy as np

    b1 = np.asarray(p2, float) - np.asarray(p1, float)
    n = np.linalg.norm(b1)
    if n == 0:
        return 0.0
    b1 /= n
    b0 = np.asarray(p0, float) - np.asarray(p1, float)
    b2 = np.asarray(p3, float) - np.asarray(p2, float)
    v = b0 - np.dot(b0, b1) * b1
    w = b2 - np.dot(b2, b1) * b1
    return math.degrees(math.atan2(np.dot(np.cross(b1, v), w), np.dot(v, w)))


def _complex_stereo_targets(metal_complex):
    """Map each ligand's carried C=C stereo to complex-global atom indices.

    The generator's embed uses a random seed, so a pendant alkene's dihedral is
    not guaranteed to land on the requested side every time (newer RDKit made
    distance geometry honor the constraint less reliably). Returning the target
    tuples lets ``generate_3d_structures`` reject conformers that embedded the
    wrong E/Z. Each tuple is ``(i, j, stereo, ref_a, ref_b)`` in the same index
    space as ``metal_complex.get_position()``.
    """
    targets = []
    try:
        idx_map = metal_complex.get_atom_indices_for_each_ligand()
        ligands = metal_complex.ligands
    except Exception:
        return targets
    for li, ligand in enumerate(ligands):
        if li >= len(idx_map):
            continue
        amap = idx_map[li]
        stereo_bonds = getattr(getattr(ligand, "molecule", None), "stereo_bonds", []) or []
        for si, sj, stereo, sra, srb in stereo_bonds:
            if max(si, sj, sra, srb) >= len(amap):
                continue
            targets.append((amap[si], amap[sj], stereo, amap[sra], amap[srb]))
    return targets


def _complex_chiral_targets(metal_complex):
    """Map each ligand's carried sp3 chirality to complex-global atom indices.

    ``enforceChirality`` makes the embed honor a chiral volume constraint, but a
    distorted or fallback embed can still land on the wrong enantiomer. These
    targets let ``generate_3d_structures`` reject such conformers. Each tuple is
    ``(center, (n0, n1, n2, n3), ChiralType)`` in ``get_position()`` index space.
    """
    targets = []
    try:
        idx_map = metal_complex.get_atom_indices_for_each_ligand()
        ligands = metal_complex.ligands
    except Exception:
        return targets
    for li, ligand in enumerate(ligands):
        if li >= len(idx_map):
            continue
        amap = idx_map[li]
        centers = getattr(getattr(ligand, "molecule", None), "chiral_centers", []) or []
        for center, nbrs, tag in centers:
            if center >= len(amap) or max(nbrs) >= len(amap):
                continue
            targets.append((amap[center], tuple(amap[k] for k in nbrs), tag))
    return targets


def _chiral_targets_satisfied(positions, targets):
    """True if every carried sp3 centre embedded with the requested handedness.

    For neighbours in the order the tag was recorded against, the signed volume
    ``(n1-n0) . ((n2-n0) x (n3-n0))`` is positive for CHI_TETRAHEDRAL_CW and
    negative for CHI_TETRAHEDRAL_CCW (measured against RDKit; pinned by
    ``test_generator_atom_chirality.test_signed_volume_sign_convention``). Like
    the E/Z check this is self-consistent with the stored neighbour order, so it
    needs no CIP perception. A near-planar centre (|volume| tiny) is treated as
    unresolved and rejected rather than guessed.
    """
    import numpy as np
    from rdkit import Chem

    for center, nbrs, tag in targets:
        try:
            p0, p1, p2, p3 = (np.asarray(positions[i], dtype=float) for i in nbrs)
        except Exception:
            continue
        volume = float(np.dot(p1 - p0, np.cross(p2 - p0, p3 - p0)))
        if abs(volume) < 0.1:  # degenerate/planar embed -- not a resolved centre
            return False
        wants_positive = tag == Chem.ChiralType.CHI_TETRAHEDRAL_CW
        if wants_positive != (volume > 0):
            return False
    return True


def _stereo_targets_satisfied(positions, targets):
    """True if every carried C=C stereo is cleanly reproduced by ``positions``.

    CIS/STEREOZ means the two reference atoms sit on the same side of the double
    bond (dihedral ~0); TRANS/STEREOE means opposite sides (dihedral ~180). The
    check is self-consistent with the stored reference atoms, so it needs no CIP
    perception.

    A genuine alkene is planar, so a correct conformer's dihedral is close to 0
    or 180. We accept only clearly-resolved geometry (|dihedral| <= 60 for cis,
    >= 120 for trans) and reject the ambiguous middle band: a distorted embed
    that lands near 90 reads as cis to a naive sign test but is perceived as the
    opposite isomer by downstream CIP (from the H positions), which is exactly
    how a wrong E/Z used to leak through.
    """
    from rdkit import Chem

    for i, j, stereo, ref_a, ref_b in targets:
        try:
            ang = abs(_dihedral_deg(positions[ref_a], positions[i], positions[j], positions[ref_b]))
        except Exception:
            continue
        wants_same_side = stereo in (Chem.BondStereo.STEREOCIS, Chem.BondStereo.STEREOZ)
        if wants_same_side:
            if ang > 30.0:  # not a clean, near-planar cis
                return False
        else:
            if ang < 150.0:  # not a clean, near-planar trans
                return False
    return True


def _greedy_enabled(ff_params):
    """SL3 greedy-placement gate: ``ff_params["greedy"]`` or ``OIN_GREEDY_PLACEMENT``.

    Default OFF -> the embed pool is byte-identical to pristine. Mirrors the
    ``_oin_direct_enabled`` (SL2) / ``clash.VDW_ACCEPTANCE_ENABLED`` env+flag pattern.
    """
    return bool(ff_params and ff_params.get("greedy")) or (
        os.environ.get("OIN_GREEDY_PLACEMENT", "0") != "0"
    )


#: ``accept_fn`` verdict meaning "acceptable as a FALLBACK, but keep looking for something
#: better". Introduced by v0.4.15 Lane 2 (``OIN_ACCEPT_STRING_EXACT``), whose predicate is
#: strictly narrower than the round-trip key: a conformer can carry the requested key and still
#: be the wrong enantiomer, because ``compare._parse_vertex_colors`` folds reflection. Rejecting
#: such a conformer outright would be a real regression -- an accepted conformer is returned as
#: the SOLE pool member, so rejecting it lets the energy-sorted pool hand back a DIFFERENT mol.
#: Returning this instead keeps the incumbent available, which is what makes a stricter
#: acceptance predicate cost latency rather than accuracy.
#:
#: A sentinel rather than ``True``/``False`` because both of those already mean something here,
#: and because a bool-returning ``accept_fn`` must keep working untouched.
ACCEPT_INCUMBENT = object()


def generate_3d_structures(
    m_smiles,
    num_conformers=1,
    optimizer=None,
    ff_params=None,
    uff_pool_size=50,
    rmsd_threshold=0.5,
    energy_threshold=2.0,
    timeout=None,
    seed=42,
    embed_time_budget=None,
    metal_complex=None,
    accept_fn=None,
    enforce_budget=None,
    incumbent_bound=None,
):
    """Generate 3D structures from an m-SMILES string.

    Returns a list of successful geometries as molecule objects, sorted by energy if an optimizer is
    used.

    accept_fn: optional predicate ``accept_fn(molecule) -> bool`` (SL1 acceptance-gate). When
    provided, the attempt loop STOPS building the pool the moment a freshly-accepted conformer
    satisfies it, and that single conformer is returned directly (skipping the pool
    sort/dedup/clash-rerank and the optimizer pass -- re-optimization could perturb it off the
    matched key). The default ``None`` leaves both attempt loops byte-identical to pristine.

    ``accept_fn`` may also return :data:`ACCEPT_INCUMBENT` -- "good enough to fall back on, not
    good enough to stop for". The pool then keeps filling, and if it exhausts without a full
    acceptance the FIRST such conformer is returned instead of the energy-sorted pool. That makes
    a stricter acceptance predicate non-regressive by construction: the incumbent is exactly the
    conformer a looser predicate would have stopped on and returned. A plain ``bool`` return is
    unaffected, so every pre-existing caller keeps its behaviour.

    ff_params: optional dict of TMCOptimizer convergence knobs
    (ff_max_iters, ff_force_tol, ff_energy_tol, d_converge, num_relaxation).

    seed: base distance-geometry seed. The embed was previously seeded from an
    unseeded ``random.randint``, so repeated runs of the same m-SMILES returned
    different structures -- and every sp3 stereocentre landed on a random
    enantiomer. Defaulting to 42 matches the rest of the project.

    embed_time_budget: optional wall-clock cap (seconds) on the attempt loop. The
    FF-only path had no bound at all -- ``timeout`` was consumed only by the ASE
    optimizer -- so a molecule whose embed never satisfies the validity checks ran
    the full ``max_attempts`` budget (ZIHGEE_comp_0: 250 attempts, ~1696 s) before
    returning nothing. The deadline is checked between attempts (never mid-embed;
    the in-flight attempt always finishes), so a molecule that DOES embed is
    unaffected -- it fills the pool and breaks first. On exhaustion the pool built
    so far is returned; an empty pool becomes the same ``[]`` as before, only fast.
    None (the default) preserves the prior unbounded behavior for direct callers.

    enforce_budget: turn ``embed_time_budget`` from a between-attempts hint into a
    BOUND (``OIN_ENFORCE_BUDGET``; ``None`` reads the lever, default OFF). Two
    changes, both no-ops when off:

      * the deadline is threaded into ``embed.get_embedding``, which checks it inside
        its two nested Python loops over ``AllChem.EmbedMolecule``. **That is where
        the time actually is** -- profiled on ``FOSNEI_comp_0`` (the corpus's 759.9 s
        worst case) at a 300 s budget, ``get_embedding`` is 61.5 s of SELF time across
        10 calls, while the CBC solve is 1.74 s (2.1%) and the ``accept_fn`` re-encode
        0.63 s (0.8%). The charter predicted the latter two; they are noise. Because
        the sink is a Python loop rather than one long native call, a deadline check
        is real enforcement and ``eps`` is one ``EmbedMolecule``, not one attempt --
        no ``fork``/``RLIMIT_CPU`` machinery is required.
      * an empty pool at the deadline raises ``BudgetExhaustedError`` instead of
        returning ``[]``.

    incumbent_bound: cap how far the pool keeps filling AFTER the first
    :data:`ACCEPT_INCUMBENT` verdict, counted in ``accept_fn`` evaluations
    (``OIN_STRING_EXACT_BOUND``; ``None`` = unbounded, the v0.4.15 behaviour). This is
    the runtime half of the ``OIN_ACCEPT_STRING_EXACT`` trade: the lever's whole cost is
    that it declines to stop the pool, and on the frozen v0.4.15 arm the 317 molecules
    that never gain consume 93.9% of its bill.

    **Bounding is answer-neutral except where it is the point.** ``incumbent_hit`` is the
    FIRST such conformer and is returned as the sole pool member regardless of how much
    longer the pool fills, so truncating changes the answer only for a molecule whose
    string-exact hit lies beyond the bound. ``incumbent_bound=0`` therefore reproduces
    the lever-OFF answer byte-for-byte -- which is the wiring gate, not a degenerate case.

    ⚠ Answer-neutral **only once an incumbent exists**. A molecule that never records one
    falls through to the energy-sorted pool below, where a shorter pool genuinely can
    select differently. The bound is armed by ``incumbent_hit``, so it cannot fire before
    then; a caller deriving a recovered-vs-bound curve must still exclude the molecules
    that never record an incumbent, and say how many it excluded.
    """
    import time

    # The PuLP/CBC bond-order solve is re-run on the identical topology across
    # every conformer attempt; a topology-keyed memo collapses those redundant
    # CBC subprocesses. Scope it to this one generation so a long in-process
    # sweep never accumulates stale topologies (the win is entirely within a
    # single molecule's attempt loop).
    clear_pulp_cache()
    # Per-generation memo for get_alternative_molecule, owned here (one fresh dict
    # per generation) and threaded into the embed calls; direct embed callers pass
    # None and recompute. Same rationale as the PuLP memo.
    alt_cache = {}

    # ``metal_complex`` may be supplied pre-built (SL2 oin-direct-winding path,
    # which constructs the complex directly from ParsedOIN with winding attached);
    # otherwise parse it from the m-SMILES string as before.
    if metal_complex is None:
        try:
            metal_complex = om.get_om_from_modified_smiles(m_smiles)
        except Exception as e:
            logger.debug(f"Failed to parse m-SMILES: {e}")
            return []

    clean_ff_params = {
        k: v
        for k, v in (ff_params or {}).items()
        if k
        not in (
            "max_attempts",
            "embed_num_threads",
            "optimize_num_workers",
            "use_kabsch",
            "kabsch_only",
            "embed_no_progress_attempts",
            "oin_direct",
            "greedy",
        )
    }
    cleaner = clean_geometry.TMCOptimizer(**clean_ff_params)
    # option 3 is the A4 rigid-placement (kabsch) embed. It is OPT-IN so the default
    # pool stays byte-identical to pristine: unset -> [0,1,2] exactly as before;
    # use_kabsch -> add 3 to the pool; kabsch_only -> isolate 3 for a clean A/B.
    options = [0, 1, 2]  # Added one more option to increase pool variety
    if ff_params:
        if ff_params.get("kabsch_only"):
            options = [3]
        elif ff_params.get("use_kabsch"):
            options = [0, 1, 2, 3]
    # SL3 greedy placement (opt-in via OIN_GREEDY_PLACEMENT / ff_params["greedy"],
    # default OFF -> byte-identical). Greedy is a *variant of the kabsch embed*
    # (option 3), so enabling it ensures option 3 is in the pool; by default it
    # enters COMPETITIVELY alongside the DG embeds ([0,1,2,3]) and the existing
    # clash-ranked pool selection keeps it only when it wins. kabsch_only/use_kabsch
    # still shape the pool as before -- greedy just changes how option 3 places.
    greedy_enabled = _greedy_enabled(ff_params)
    if greedy_enabled and 3 not in options:
        options = options + [3]
    scales = [0.8, 0.9, 1.0, 1.1, 1.2]

    # Target number of initial structures to generate
    target_pool = uff_pool_size
    successful_mols = []
    # Carried C=C (E/Z) and sp3 chirality targets, and a pool of otherwise-valid
    # conformers that embedded the wrong stereochemistry (kept only as a
    # last-resort fallback so a stubborn embed never hard-fails generation).
    stereo_targets = _complex_stereo_targets(metal_complex)
    chiral_targets = _complex_chiral_targets(metal_complex)
    stereo_rejects = []

    import itertools

    combinations = list(itertools.product(scales, options))

    default_max = max(target_pool * 5, 250)
    max_attempts = ff_params.get("max_attempts", default_max) if ff_params else default_max

    # Opt-in C++-parallel embed. num_threads == 1 (the default) keeps the serial
    # attempt loop below, which is byte-identical to pristine. num_threads != 1
    # (e.g. 0 = all cores) switches to the batched EmbedMultipleConfs path: RDKit
    # embeds a whole wave of conformers per (scale, option) in parallel C++,
    # releasing the GIL. That path is faster but NOT byte-identical (it samples
    # conformers differently), so it is gated behind this explicit opt-in and
    # validated by the accuracy gate rather than by byte-identity.
    num_threads = int(ff_params.get("embed_num_threads", 1)) if ff_params else 1

    # Bound the work the stereo filters can add, keyed on REJECTIONS rather than on
    # attempts. A rejected conformer does not fill the pool, so without a cap a
    # molecule whose embed rarely satisfies its constraints never reaches
    # target_pool and runs the FULL attempt budget: AFECIZ (a chelate imine whose
    # C=N the embed seldom lands planar, because useBasicKnowledge=False gives
    # distance geometry no double-bond planarity term) went from 565s unconstrained
    # to >27 min, all of it spent rejecting.
    #
    # Counting rejections, not attempts, leaves healthy molecules untouched -- they
    # reject nothing, so the cap never fires -- and does not penalise a molecule
    # that burns attempts for unrelated reasons (a valence exception, a haptic scale
    # that will not embed). Once the cap is hit we keep whatever satisfied the
    # constraints; if nothing did, the stereo_rejects fallback below returns the same
    # conformer the unfiltered code would have returned, so this is never worse.
    reject_budget = max(target_pool * 2, 25)

    # Wall-clock deadline for the attempt loop (see docstring). Checked between
    # attempts so a molecule that embeds cleanly is never interrupted; it exists
    # only to keep a molecule whose embed never validates from running the full
    # max_attempts budget (ZIHGEE ~1696 s) before giving up.
    deadline = (time.monotonic() + embed_time_budget) if embed_time_budget else None

    # OIN_ENFORCE_BUDGET. Read through lever_enabled(), never os.environ.get -- "0" is a
    # non-empty string and would ENABLE the lever; that trap cost 23 test failures across
    # two promotions and test_levers::TestNoTestUnsetsAPromotedLever lints for it. An
    # explicit argument wins over the environment so the unit tests do not have to mutate
    # process state.
    if enforce_budget is None:
        from ..oin.levers import lever_enabled

        enforce_budget = lever_enabled("OIN_ENFORCE_BUDGET")
    # The deadline only becomes a bound if there IS one; with no budget there is nothing
    # to enforce and this stays byte-identical to pristine.
    embed_deadline = deadline if (enforce_budget and deadline is not None) else None

    # Opt-in no-acceptance-progress cutoff (v0.4.4 SL4), gated OFF by default so the
    # pool-fill path stays byte-identical. When set, the loop gives up after this many
    # consecutive attempts that produced an embed but grew the pool by nothing -- the
    # OSIHUU pattern where the vdW acceptance gate rejects every candidate, so a
    # molecule already at its embed budget would otherwise burn the full max_attempts /
    # wall-clock deadline before returning empty. Keyed on acceptance progress (not
    # wall-clock), so it is deterministic and never cuts off a still-progressing embed.
    # Gated below on had_nonstructural_embed so a uniformly-structural pool still
    # surfaces via StructuralAssemblyError rather than this generic cutoff.
    no_progress_limit = ff_params.get("embed_no_progress_attempts") if ff_params else None
    if no_progress_limit is None and os.environ.get("OIN_EMBED_NO_PROGRESS"):
        no_progress_limit = int(os.environ["OIN_EMBED_NO_PROGRESS"])
    no_progress = 0

    def _try_accept(positions, scale):
        """Clean one embedded conformer, stereo-filter it, and file the result.

        Files into ``successful_mols`` or ``stereo_rejects``. Shared verbatim by the
        serial and batched fill loops so both apply identical acceptance. Returns the
        accepted ``Molecule`` (the one appended to ``successful_mols``) or ``None`` -- the
        SL1 early-exit hook keys off this to test the fresh conformer against ``accept_fn``.
        """
        if positions is None:
            return None
        tmp_complex = metal_complex.copy()
        tmp_complex.set_position(positions)

        # cleaner.clean_geometry will print logs, could be silenced later
        success = cleaner.clean_geometry(tmp_complex, scale)

        if success:
            position = tmp_complex.get_position()
            wrong_ez = stereo_targets and not _stereo_targets_satisfied(position, stereo_targets)
            wrong_chirality = chiral_targets and not _chiral_targets_satisfied(
                position, chiral_targets
            )
            if wrong_ez or wrong_chirality:
                # Right topology, wrong stereochemistry -- keep as fallback only.
                stereo_rejects.append(tmp_complex.get_molecule())
                return None
            mol = tmp_complex.get_molecule()
            successful_mols.append(mol)
            return mol
        return None

    # Surface *why* an exhausted pool failed (see StructuralAssemblyError) only
    # when the complex is *uniformly* unassemblable -- i.e. EVERY attempt raised a
    # structural error. `had_nonstructural_embed` records that at least one attempt
    # got past sanitization (an embed returned, or failed only transiently); when
    # it is True the emptiness is not purely structural (some option assembled but
    # its embed did not validate), so we degrade to a generic no_conformers.
    last_structural_error = None
    had_nonstructural_embed = False

    # SL1 accept-first early-exit. ``early_hit`` holds the first conformer that satisfies
    # ``accept_fn`` (re-encodes to the requested key); once set, both fill loops stop and it is
    # returned directly. ``accept_fn is None`` (default) => this never fires => byte-identical.
    early_hit = None

    # The first conformer ``accept_fn`` called ACCEPT_INCUMBENT: acceptable to fall back on, not
    # to stop for. Returned only if the pool exhausts with no ``early_hit``. Stays None unless an
    # accept_fn actually returns the sentinel, so this is byte-identical for every other caller.
    incumbent_hit = None

    # OIN_STRING_EXACT_BOUND bookkeeping. ``_since_incumbent`` counts accept_fn evaluations made
    # AFTER the incumbent was recorded; ``bounded_stop`` says the bound ended the search.
    #
    # 🔴 A SEPARATE FLAG, NOT ``early_hit``. Stopping by assigning ``early_hit`` would look like it
    # works -- both loops break on it -- and would return the conformer through the `return
    # [early_hit]` path at the top, which means "a conformer INDEPENDENTLY re-encoded to the
    # requested key and we stopped for it". That is precisely the claim the bound must NOT make:
    # the conformer it stops on is the incumbent, i.e. the one whose string did NOT match. The
    # verdict would be indistinguishable in the return value and wrong in the telemetry, which is
    # this project's recurring failure mode rather than a hypothetical.
    _since_incumbent = 0
    bounded_stop = False

    def _file_and_maybe_stop(positions, scale):
        """File a conformer, then test it against ``accept_fn``; True => stop building the pool."""
        nonlocal early_hit, incumbent_hit, _since_incumbent, bounded_stop
        accepted = _try_accept(positions, scale)
        if accept_fn is not None and accepted is not None and early_hit is None:
            # OIN_ENFORCE_BUDGET: do not START an accept_fn call we cannot afford.
            #
            # ⚠ THE SECOND MEASUREMENT CORRECTED THE FIRST. On FOSNEI_comp_0 this predicate
            # is 0.63 s of 82.4 s (0.8%) and looked like noise. On VAFMIA_comp_0 with a 30 s
            # budget it is 48.5 s of 78.6 s -- the single largest consumer -- because
            # build_contract_mol -> _template_sp3_label lands in
            # chirality._reparse_cip_label_once at ~2.4 s a call. The cost is BIMODAL, not
            # small: which function dominates is a property of the molecule, not of the
            # pipeline. A bound threaded only into whichever function profiled expensive
            # last is not a bound, which is exactly what the first version of this lever was.
            if embed_deadline is not None and time.monotonic() > embed_deadline:
                logger.debug("budget reached; skipping the accept_fn re-encode")
                return False
            # Whether an incumbent existed BEFORE this evaluation. The distinction is the whole
            # off-by-one: the evaluation that RECORDS the incumbent is not one of the extra
            # evaluations the bound is meant to cap, so it must not advance the counter.
            _had_incumbent = incumbent_hit is not None
            try:
                verdict = accept_fn(accepted)
                if verdict is ACCEPT_INCUMBENT:
                    # Keep filling, but remember the FIRST one -- it is the conformer a looser
                    # predicate would have stopped on, so returning it later reproduces the
                    # pre-lever answer exactly.
                    if incumbent_hit is None:
                        incumbent_hit = accepted
                        _telemetry.record("pool.accept_incumbent_recorded")
                elif verdict:
                    early_hit = accepted
                    # THE ORDINAL IS WHY ONE RUN SUFFICES. `min_bound` is the smallest
                    # incumbent_bound at which this molecule is still recovered, recorded so the
                    # whole recovered-vs-bound curve is arithmetic over a single unbounded run
                    # rather than one full run per candidate bound.
                    #
                    # It is recorded as the ANSWER (the minimum bound) rather than as the raw
                    # counter precisely so no consumer has to re-derive the +1 -- the reader of a
                    # frozen JSON six months from now cannot see this loop.
                    _telemetry.record(
                        "pool.string_exact_hit",
                        min_bound=int(_since_incumbent + 1) if _had_incumbent else 0,
                        had_incumbent=bool(_had_incumbent),
                    )
                    return True
                # Count only evaluations made once an incumbent ALREADY existed -- before that
                # there is nothing to fall back on and the bound must not be able to fire.
                if incumbent_hit is not None:
                    if _had_incumbent:
                        _since_incumbent += 1
                    if incumbent_bound is not None and _since_incumbent >= incumbent_bound:
                        # bound=0 lands here on the recording evaluation itself (0 >= 0), which
                        # is what makes it byte-identical to the lever-OFF arm.
                        bounded_stop = True
                        _telemetry.record("pool.incumbent_bound_stop", bound=int(incumbent_bound))
                        return True
            except Exception:
                logger.debug("accept_fn raised on a conformer; ignoring", exc_info=True)
        return False

    # Attempts actually spent, for the eta runtime question. A PLAIN COUNTER, not a
    # degradation site: the `adapter.early_exit_*` counters cannot answer it because
    # `_select_by_geometry(..., early_exit=False)` is the default and that block never runs, so a
    # telemetry sweep came back with every site at zero
    # (docs/agentic-notes/v0.4.6/V046_HFAITHFUL_FINDINGS.md).
    #
    # What this discriminates: flat attempts across eta and size-matched non-eta molecules means the
    # eta 3-6x penalty is COST-PER-ATTEMPT (a profiling target); systematically higher attempts for
    # eta means it is ACCEPTANCE-LIMITED (the embed rarely produces the requested ring face).
    # Recorded via _telemetry, so it is a no-op unless OIN_TELEMETRY=1 and a collecting() context is
    # active -- byte-identical otherwise, and it consumes no randomness.
    _attempts_spent = 0

    def _record_attempts(n_accepted):
        """Observation only: a no-op unless OIN_TELEMETRY=1 with a collecting() context active.

        Called at EVERY return of this function, not just the last one. The first version recorded
        only at the final return and never fired: the early-exit path returns at
        `return [early_hit]` well before it -- exactly the molecules whose count matters most.
        """
        _telemetry.record(
            "pool.attempts_spent",
            attempts=int(_attempts_spent),
            accepted=int(n_accepted),
            target_pool=int(target_pool),
        )

    if num_threads == 1:
        # Serial attempt loop -- byte-identical to pristine.
        for i in range(max_attempts):
            _attempts_spent = i
            if len(successful_mols) >= target_pool:
                break
            if deadline is not None and time.monotonic() > deadline:
                logger.debug(
                    f"Embed wall-clock budget ({embed_time_budget}s) reached after "
                    f"{i} attempt(s) with {len(successful_mols)} conformer(s) "
                    f"({len(stereo_rejects)} stereo-reject fallback(s)); stopping rather "
                    f"than exhausting {max_attempts} attempts."
                )
                break
            if len(stereo_rejects) >= reject_budget:
                logger.debug(
                    f"Stereo reject budget ({reject_budget}) reached with "
                    f"{len(successful_mols)} conformer(s) satisfying the requested stereo; "
                    f"stopping rather than exhausting {max_attempts} attempts."
                )
                break
            if no_progress_limit and no_progress >= no_progress_limit and had_nonstructural_embed:
                logger.debug(
                    f"No acceptance progress in {no_progress} consecutive attempt(s) "
                    f"with {len(successful_mols)} conformer(s); the acceptance gate is "
                    f"rejecting every embed -- stopping rather than exhausting "
                    f"{max_attempts} attempts."
                )
                break
            scale, option = combinations[i % len(combinations)]

            # A single scale/option combo can raise inside the embed (e.g. an
            # RDKit valence exception on a dative donor); skip it rather than
            # letting one bad combo abort the whole pool.
            try:
                # Distinct-but-reproducible seed per attempt: the pool keeps its
                # variety, but the same m-SMILES always yields the same conformers.
                positions = embed.get_embedding(
                    metal_complex,
                    scale,
                    option,
                    align=True,
                    seed=seed + i * 1009,
                    alt_cache=alt_cache,
                    greedy=greedy_enabled,
                    deadline=embed_deadline,
                )
                had_nonstructural_embed = True  # got past sanitization (positions or None)
            except _STRUCTURAL_EMBED_ERRORS as e:
                # Structural / non-transient (an over-valent dative donor or an
                # unfilled geometry -- see _STRUCTURAL_EMBED_ERRORS). Repeats for
                # every seed/scale, so retrying cannot help -- remember it, keep
                # trying the remaining (scale, option) combos in case another
                # assembles, and let a *uniformly* structural pool surface WHY
                # (StructuralAssemblyError below).
                logger.debug(f"Structural embed failure (scale={scale}, option={option}): {e}")
                last_structural_error = e
                positions = None
            except Exception as e:
                logger.debug(f"Embedding failed (scale={scale}, option={option}): {e}")
                _telemetry.record("pool.blanket_exception", exc=type(e).__name__)
                had_nonstructural_embed = True
                positions = None
            pool_before = len(successful_mols)
            if _file_and_maybe_stop(positions, scale):
                break
            no_progress = 0 if len(successful_mols) > pool_before else no_progress + 1
    else:
        # Batched, C++-parallel embed path (opt-in via embed_num_threads != 1).
        # Each (scale, option) combo embeds a wave of conformers in one
        # EmbedMultipleConfs call, parallelized across num_threads cores. Batch
        # size is chosen so a few feasible combos fill the pool (which also skips
        # the infeasible high-scale combos most of the time). Haptic complexes fall
        # back to the serial embed per combo (their cmap needs the scales_for_haptic
        # sweep).
        batch_k = max(4, -(-target_pool // 3))  # ceil(target_pool / 3), >= 4
        for i in range(max_attempts):
            if len(successful_mols) >= target_pool:
                break
            if deadline is not None and time.monotonic() > deadline:
                logger.debug(
                    f"Embed wall-clock budget ({embed_time_budget}s) reached after "
                    f"{i} batch(es) with {len(successful_mols)} conformer(s) "
                    f"({len(stereo_rejects)} stereo-reject fallback(s)); stopping rather "
                    f"than exhausting {max_attempts} batches."
                )
                break
            if len(stereo_rejects) >= reject_budget:
                logger.debug(
                    f"Stereo reject budget ({reject_budget}) reached with "
                    f"{len(successful_mols)} conformer(s) satisfying the requested stereo; "
                    f"stopping rather than exhausting {max_attempts} batches."
                )
                break
            if no_progress_limit and no_progress >= no_progress_limit and had_nonstructural_embed:
                logger.debug(
                    f"No acceptance progress in {no_progress} consecutive batch(es) "
                    f"with {len(successful_mols)} conformer(s); the acceptance gate is "
                    f"rejecting every embed -- stopping rather than exhausting "
                    f"{max_attempts} batches."
                )
                break
            scale, option = combinations[i % len(combinations)]
            pool_before = len(successful_mols)
            try:
                batch = embed.get_embeddings_batch(
                    metal_complex,
                    scale,
                    option,
                    num_confs=batch_k,
                    num_threads=num_threads,
                    align=True,
                    seed=seed + i * 1009,
                    alt_cache=alt_cache,
                )
                had_nonstructural_embed = True  # got past sanitization (positions or haptic None)
            except _STRUCTURAL_EMBED_ERRORS as e:
                # Structural (see the serial path). get_embeddings_batch returns
                # None -- it does not raise -- to signal a non-batchable haptic
                # complex, so an exception here is a genuine failure, not the
                # haptic fallback.
                logger.debug(f"Structural embed failure (scale={scale}, option={option}): {e}")
                last_structural_error = e
                batch = None
            except Exception as e:
                logger.debug(f"Embedding failed (scale={scale}, option={option}): {e}")
                had_nonstructural_embed = True
                batch = None
            if batch is None:
                # Haptic complex -- not batchable; fall back to the serial embed.
                try:
                    positions = embed.get_embedding(
                        metal_complex,
                        scale,
                        option,
                        align=True,
                        seed=seed + i * 1009,
                        alt_cache=alt_cache,
                        greedy=greedy_enabled,
                        deadline=embed_deadline,
                    )
                    had_nonstructural_embed = True  # got past sanitization
                except _STRUCTURAL_EMBED_ERRORS as e:
                    # Structural (see the serial path) -- record and continue.
                    logger.debug(f"Structural embed failure (scale={scale}, option={option}): {e}")
                    last_structural_error = e
                    positions = None
                except Exception as e:
                    logger.debug(f"Embedding failed (scale={scale}, option={option}): {e}")
                    had_nonstructural_embed = True
                    positions = None
                if _file_and_maybe_stop(positions, scale):
                    break
                no_progress = 0 if len(successful_mols) > pool_before else no_progress + 1
                continue
            for positions in batch:
                if len(successful_mols) >= target_pool:
                    break
                if _file_and_maybe_stop(positions, scale):
                    break
            no_progress = 0 if len(successful_mols) > pool_before else no_progress + 1
            # ⚠ `bounded_stop` as well as `early_hit`. The inner `for positions in batch` loop
            # above breaks only ITSELF, so without this the outer loop would embed a fresh batch
            # and the bound would cap nothing on the batched path -- a lever that fires, logs, and
            # saves no time, which is exactly the shape of v0.4.12's ETA_EARLY_EXIT.
            if early_hit is not None or bounded_stop:
                break

    if early_hit is not None:
        # SL1 accept-first early-exit fired: a conformer independently re-encoded to the
        # requested key. Return it directly, bypassing the pool sort/dedup/clash-rerank AND the
        # optimizer pass below -- re-optimization could perturb it off the matched key, and the
        # accept stamp was taken on this exact geometry. Unreachable when accept_fn is None.
        _record_attempts(1)
        return [early_hit]

    if incumbent_hit is not None:
        # The pool exhausted without a full acceptance, but a conformer met the OLD, looser
        # predicate. Return it on exactly the same terms `early_hit` would have been returned on
        # -- which is the point: this is byte-identical to what the pre-lever run produced, so a
        # stricter acceptance predicate can only cost latency here, never accuracy.
        logger.debug("accept_fn found no full match; returning the recorded incumbent conformer")
        _telemetry.record("pool.accept_incumbent_returned")
        _record_attempts(1)
        return [incumbent_hit]

    if not successful_mols and stereo_rejects:
        # No embed reproduced the requested stereochemistry within the attempt
        # budget; return the best available rather than nothing (non-regressive
        # vs. the prior, unfiltered behavior).
        logger.debug("WARNING: no conformer reproduced the requested stereo; using best available.")
        _telemetry.record("pool.stereo_fallback_wrong_isomer", n_rejects=len(stereo_rejects))
        successful_mols = stereo_rejects

    if not successful_mols:
        # OIN_ENFORCE_BUDGET: the budget ran out with NOTHING to show. Say so, typed,
        # rather than degrading to the generic `[]` that becomes "no_conformers" three
        # frames up. StructuralAssemblyError is checked FIRST because a uniformly
        # structural failure would have happened at any budget, so it is the better
        # diagnosis; "out of time" is what is left when nothing more specific applies.
        # Both sit after the stereo_rejects fallback above, so a molecule that produced
        # only wrong-stereo conformers still returns them.
        if last_structural_error is not None and not had_nonstructural_embed:
            # EVERY attempt failed structurally (no option ever got past
            # sanitization) and nothing filled the pool: the complex is uniformly
            # unassemblable, so surface WHY instead of a generic no_conformers.
            # Chain the underlying error for the full traceback.
            raise StructuralAssemblyError(
                "Could not assemble a valid 3D complex: every embed attempt failed "
                f"with the same structural error ({type(last_structural_error).__name__}: "
                f"{last_structural_error})"
            ) from last_structural_error
        if enforce_budget and deadline is not None and time.monotonic() > deadline:
            spent = embed_time_budget + (time.monotonic() - deadline)
            raise BudgetExhaustedError(
                f"Generation budget exhausted: {embed_time_budget}s requested, {spent:.1f}s "
                f"spent, {_attempts_spent} attempt(s), pool empty "
                f"({len(stereo_rejects)} stereo-reject(s)). This is a BOUND, not an "
                "assembly failure -- the molecule may well assemble under a larger budget "
                "(ULODUU does at 60s and does not at 30s)."
            )
        return []

    # Sort by UFF energy if available (handle None values safely)
    successful_mols.sort(
        key=lambda m: (
            getattr(m, "energy", None) if getattr(m, "energy", None) is not None else float("inf")
        )
    )

    # Deduplicate
    dedup_mols = []
    for mol in successful_mols:
        is_unique = True
        for acc_mol in dedup_mols:
            rmsd = calculate_heavy_atom_rmsd(mol, acc_mol)
            e1 = getattr(mol, "energy", None)
            e1 = e1 if e1 is not None else 0.0
            e2 = getattr(acc_mol, "energy", None)
            e2 = e2 if e2 is not None else 0.0
            if rmsd < rmsd_threshold and abs(e1 - e2) <= energy_threshold:
                is_unique = False
                break
        if is_unique:
            dedup_mols.append(mol)

    successful_mols = dedup_mols

    # Re-rank the deduplicated pool by whole-complex vdW clash -- gated OFF by default
    # (clash.VDW_ACCEPTANCE_ENABLED). Stable, so UFF energy still orders conformers within
    # equal clash. Keeps the least-clashing conformers inside the ``num_conformers`` slice
    # that feeds the optimizer/selection. Off by default because on the current pool the
    # least-clashing conformer has loosened coordination (see clash.py); disabled -> the
    # pool keeps its pre-A3 energy order (byte-identical).
    if clash.VDW_ACCEPTANCE_ENABLED:
        successful_mols.sort(key=clash.mol_clash_count)

    if optimizer:
        from .ml_optimizer import ASEOptimizer

        # We explicitly do NOT catch initialization exceptions here.
        # If the user asks for an optimizer and it fails to load,
        # we want to fail loudly rather than silently falling back to FF.
        opt = ASEOptimizer(method=optimizer, timeout=timeout)

        if opt:
            mols_to_optimize = successful_mols[:num_conformers]
            # Optimize the pooled conformers concurrently. Each opt.optimize deep-copies
            # its input and runs xtb in its own TemporaryDirectory subprocess pinned to a
            # single OpenMP thread (see ml_optimizer), so the calls share no state and
            # each xtb result is load-independent/deterministic. ThreadPoolExecutor.map
            # preserves INPUT ORDER, so the stable sort-by-energy below is unchanged.
            if len(mols_to_optimize) > 1:
                from concurrent.futures import ThreadPoolExecutor  # noqa: PLC0415

                # Worker cap: default to a safe count BELOW the core total (leave the
                # machine headroom -- each worker runs a single-threaded xtb), overridable
                # via ff_params["optimize_num_workers"]. Bounded by the number of
                # conformers either way.
                n_cpu = os.cpu_count() or 2
                default_workers = max(1, n_cpu - 2)
                requested = (ff_params or {}).get("optimize_num_workers")
                workers = int(requested) if requested else default_workers
                max_workers = max(1, min(len(mols_to_optimize), workers))
                with ThreadPoolExecutor(max_workers=max_workers) as ex:
                    opt_results = list(ex.map(opt.optimize, mols_to_optimize))
            else:
                opt_results = [opt.optimize(mol) for mol in mols_to_optimize]

            optimized_mols = []
            for (success, energy, new_mol), mol in zip(opt_results, mols_to_optimize):
                if success:
                    optimized_mols.append((energy, new_mol))
                else:
                    # Keep the original if optimization fails but penalize its rank
                    optimized_mols.append((float("inf"), mol))

            # Sort by energy
            optimized_mols.sort(key=lambda x: x[0])
            successful_mols = [m[1] for m in optimized_mols]

    _record_attempts(len(successful_mols))
    return successful_mols[:num_conformers]


def get_xyz_string(molecule):
    """Returns the XYZ string format for a given generated molecule."""
    atom_list = molecule.atom_list
    lines = [
        str(len(atom_list)),
        f"Generated by MetalloGen-3D, Charge: {molecule.chg}, "
        f"Multiplicity: {molecule.multiplicity}",
    ]
    for atom in atom_list:
        element = atom.get_element()
        coord = atom.get_coordinate()
        lines.append(f"{element:<3} {coord[0]:>12.8f} {coord[1]:>12.8f} {coord[2]:>12.8f}")
    return "\n".join(lines)
