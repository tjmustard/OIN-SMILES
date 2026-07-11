from concurrent.futures import ProcessPoolExecutor
from concurrent.futures import TimeoutError as FuturesTimeoutError

import scine_molassembler as masm
from rdkit import Chem

from .molassembler_adapter import (
    GeneratedStructure,
    MolassemblerAdapter,
    MolassemblerTimeoutError,
)
from .oin_parser import (
    OINParser,
    _extract_oin_constraints,
    construct_molassembler_mol,
    tokenize_unsanitized_smiles,
)

__all__ = ["OIN3DGenerator", "MolassemblerTimeoutError", "GeneratedStructure"]


def _dg_worker(mol_masm: masm.Molecule, seed: int = 42) -> str:
    """Module-level DG worker for ProcessPoolExecutor.

    Generates a conformer via Molassembler distance geometry and returns XYZ block.
    Runs in a separate process to isolate GIL-holding C++ code.

    Returns:
    -------
    str
        XYZ block string

    Raises:
    ------
    RuntimeError
        If DG generation fails
    """
    try:
        confs = masm.dg.generate_conformation(mol_masm, seed=seed)
        # Write to temporary file and read back as string
        import os
        import tempfile

        with tempfile.NamedTemporaryFile(suffix=".xyz", delete=False) as f:
            tmp_path = f.name

        try:
            masm.io.write(tmp_path, mol_masm, confs)
            with open(tmp_path) as f:
                xyz_str = f.read()
            return xyz_str
        finally:
            if tmp_path:
                try:
                    os.unlink(tmp_path)
                except OSError:
                    pass
    except Exception as e:
        raise RuntimeError(f"DG generation failed: {e}") from e


def _translate_eta_vertex_to_atoms(
    vertex_index: int,
    fragment_rank: int,
    frag_to_atom: dict,
) -> list:
    """Translate OIN vertex index (slot assignment) to atom indices for eta bonding.

    Maps a fragment rank at a polyhedral slot to the actual atom indices in the
    connected SMILES that should form eta bonds to the metal.

    For most ligands (monodentate), this returns the first atom of the fragment.
    For eta ligands (Cp, arene), this returns all atoms of the fragment.

    Parameters
    ----------
    vertex_index : int
        The slot index from the OIN string (e.g., 0, 1, 2 for the ligand slots).
        Currently used only for future multi-eta extensibility.
    fragment_rank : int
        The fragment rank (OIN's fragment-rank space). Must be in frag_to_atom.
    frag_to_atom : dict
        Mapping from fragment rank to list of atom indices in connected SMILES.
        Example: {0: [0], 1: [1], 2: [2]} for cisplatin.

    Returns:
    -------
    list[int]
        List of atom indices in the connected SMILES that should bond to the metal
        at this slot. Length 1 for monodentate; length > 1 for eta ligands.

    Raises:
    ------
    ValueError
        If fragment_rank is not in frag_to_atom.
    """
    if fragment_rank not in frag_to_atom:
        raise ValueError(
            f"Fragment rank {fragment_rank} not found in fragment-to-atom mapping. "
            f"Valid ranks: {sorted(frag_to_atom.keys())}"
        )

    atom_indices = frag_to_atom[fragment_rank]
    return atom_indices


def parse_oin_direct(oin_smiles: str) -> GeneratedStructure:
    """Parse OIN-SMILES and generate 3D coordinates via direct pipeline.

    Pipeline: regex preprocessing → AST tokenization → Molassembler instantiation → DG

    Parameters
    ----------
    oin_smiles:
        OIN-SMILES string in v3.6 inline format (e.g., "[Pd_SQP].[Cl]{0}.[Cl]{1}")

    Returns:
    -------
    GeneratedStructure
        xyz: XYZ block string with 3D coordinates
        mol: RDKit Mol with bond topology and 3D conformer (may be None)

    Raises:
    ------
    ValueError
        If OIN parsing or Molassembler construction fails
    TimeoutError
        If DG conformer generation exceeds timeout
    """
    try:
        # Step 1: Regex preprocessing
        stripped_smiles, constraints, frag_to_atom = _extract_oin_constraints(oin_smiles)

        # Build connected SMILES for Molassembler (join metal to ligands)
        # Track atom offsets to update frag_to_atom mapping for eta bonds
        fragments = stripped_smiles.split(".")
        updated_frag_to_atom = {}

        if len(fragments) == 1:
            # Already connected
            connected_smiles = stripped_smiles
            updated_frag_to_atom = frag_to_atom
        else:
            # Join all fragments to the metal (atom 0)
            # Metal is fragments[0], ligands are fragments[1:]
            metal_frag = fragments[0]
            rw_mol = Chem.RWMol()

            # Add metal atom
            metal_mol = Chem.MolFromSmiles(metal_frag, sanitize=False)
            if metal_mol is None:
                raise ValueError(f"Cannot parse metal fragment: {metal_frag}")
            metal_atom = metal_mol.GetAtomWithIdx(0)
            rw_mol.AddAtom(Chem.Atom(metal_atom.GetAtomicNum()))
            rw_mol.GetAtomWithIdx(0).SetFormalCharge(metal_atom.GetFormalCharge())

            # Fragment 0 (metal) → atom 0
            updated_frag_to_atom[0] = [0]

            # Add ligand atoms and bonds, tracking atom offsets
            atom_offset = 1
            for frag_rank, lig_smiles in enumerate(fragments[1:], start=1):
                lig_mol = Chem.MolFromSmiles(lig_smiles, sanitize=False)
                if lig_mol is None:
                    raise ValueError(f"Cannot parse ligand fragment: {lig_smiles}")

                # Track atoms for this fragment in the new connected SMILES
                lig_atom_indices = []

                for atom in lig_mol.GetAtoms():
                    new_atom = Chem.Atom(atom.GetAtomicNum())
                    new_atom.SetFormalCharge(atom.GetFormalCharge())
                    rw_mol.AddAtom(new_atom)
                    lig_atom_indices.append(atom_offset + atom.GetIdx())

                for bond in lig_mol.GetBonds():
                    i = atom_offset + bond.GetBeginAtomIdx()
                    j = atom_offset + bond.GetEndAtomIdx()
                    rw_mol.AddBond(i, j, bond.GetBondType())

                # Add bond from metal to first binding atom in ligand (use SINGLE for Molassembler)
                rw_mol.AddBond(0, atom_offset, Chem.BondType.SINGLE)

                # Update mapping: fragment_rank (in OIN space) → atom indices
                # (in connected SMILES space)
                updated_frag_to_atom[frag_rank] = lig_atom_indices
                atom_offset += lig_mol.GetNumAtoms()

            connected_smiles = Chem.MolToSmiles(rw_mol.GetMol(), isomericSmiles=True)

        # Step 2: AST tokenization on connected SMILES
        atoms, bonds = tokenize_unsanitized_smiles(connected_smiles)

        # Step 3: Molassembler instantiation (pass updated fragment mapping for eta bonds)
        mol_masm = construct_molassembler_mol(
            atoms, bonds, constraints, mol_rdkit=None, frag_to_atom=updated_frag_to_atom
        )

        # Step 4: DG conformer generation with timeout (10s default)
        with ProcessPoolExecutor(max_workers=1) as executor:
            future = executor.submit(_dg_worker, mol_masm, seed=42)
            try:
                xyz = future.result(timeout=10.0)
            except FuturesTimeoutError:
                raise TimeoutError(
                    "Conformer generation exceeded 10s; complex geometry may not be supported"
                )

        # Step 5: Reconstruct RDKit mol for return (unsanitized, preserving @/@@ tags)
        mol_rdkit = Chem.MolFromSmiles(connected_smiles, sanitize=False)
        if mol_rdkit is not None:
            try:
                Chem.SanitizeMol(mol_rdkit)
            except Exception:
                pass  # Best effort; keep unsanitized if sanitization fails

        return GeneratedStructure(xyz=xyz, mol=mol_rdkit)

    except TimeoutError:
        raise  # Re-raise timeout as-is
    except Exception as e:
        raise ValueError(f"Failed to generate 3D structure from OIN-SMILES: {e}") from e


class OIN3DGenerator:
    """Generates 3D XYZ structures from OIN-SMILES strings.

    **Warning:** This class is an internal implementation detail of OIN-SMILES.
    Its API is subject to change without notice. Users should prefer
    SMILESToXYZ.convert() instead.

    Uses direct parser pipeline: regex → AST → Molassembler → DG.
    """

    def __init__(
        self,
        timeout: int = 300,
        dg_strategy: str = "single",
        ensemble_size: int = 10,
        engine: str = "metallogen",
        optimizer: str | None = "xtb",
        ff_preset: str | None = None,
        ff_params: dict | None = None,
        seed: int = 42,
    ) -> None:
        """Initialize the generator with a parser and a generation backend.

        ``engine="metallogen"`` (default) uses the MetalloGen dummy-metal + CoordMap
        embed backend (``oinsmiles.generator3d``), refined by the ``optimizer`` (default
        ``"xtb"`` standard g-xTB; pass ``"ff"``/``None`` for the FF-only path, or
        ``"mace-omol-0-extra-large-1024"`` / ``"mace-omol25"`` for higher accuracy).
        The MACE models require ``mace-torch`` + the model weights and fail loudly if
        unavailable — use ``optimizer="ff"`` in environments without them.
        ``engine="legacy"`` uses the Molassembler/stitch adapter (the
        reference for Zone-A P stereo enforcement; ``optimizer`` is ignored there). Both
        expose ``generate(parsed)``, so the ``generate()`` delegation below is identical.

        ``seed`` is the base ETKDG seed for the metallogen engine; a fixed value
        (default 42) makes generation reproducible. It is honored by
        ``engine="metallogen"`` only — the legacy adapter is already deterministic
        and ignores it.
        """
        self.parser = OINParser()
        self.engine = engine
        if engine == "metallogen":
            from .metallogen_adapter import MetalloGenAdapter

            self.adapter = MetalloGenAdapter(
                timeout=timeout,
                dg_strategy=dg_strategy,
                ensemble_size=ensemble_size,
                optimizer=optimizer,
                ff_preset=ff_preset,
                ff_params=ff_params,
                seed=seed,
            )
        elif engine == "legacy":
            self.adapter = MolassemblerAdapter(
                timeout=timeout,
                dg_strategy=dg_strategy,
                ensemble_size=ensemble_size,
            )
        else:
            raise ValueError(f"Unknown engine {engine!r}; expected 'legacy' or 'metallogen'")

    def generate(self, oin_string: str) -> GeneratedStructure:
        """Convert an OIN-SMILES string to a 3D structure.

        **Internal implementation detail.** Not part of stable public API.
        Use SMILESToXYZ.convert() instead.

        Parameters
        ----------
        oin_string:
            OIN-SMILES string (V3.0 inline or V2.4 sidecar format).

        Returns:
        -------
        GeneratedStructure
            ``.xyz`` contains the XYZ block string.
            ``.mol`` contains an RDKit Mol with bond connectivity and a 3D
            conformer (None if connectivity could not be determined).

        Raises:
        ------
        ValueError
            If OIN parsing or structure generation fails
        TimeoutError
            If conformer generation exceeds timeout
        """
        # Use legacy parser (which integrates direct parser components internally)
        # The OINParser.parse() method uses _extract_oin_constraints(),
        # tokenize_unsanitized_smiles(), and construct_molassembler_mol()
        # as part of its processing pipeline for v3.6 inline format.
        parsed = self.parser.parse(oin_string)
        return self.adapter.generate(parsed)
