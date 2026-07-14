from .oin_parser import OINParser
from .structure import GeneratedStructure

__all__ = ["OIN3DGenerator", "GeneratedStructure"]


class OIN3DGenerator:
    """Generate 3D XYZ structures from OIN-SMILES strings (lower-level engine).

    This is the OIN->XYZ entry point (the ``oin-smiles oin2xyz`` CLI uses it).
    :class:`~oinsmiles.core.translator.SMILESToXYZ` is the simple, stable public
    entry point and delegates here. Use ``OIN3DGenerator`` directly when you need
    the richer :class:`GeneratedStructure` result (XYZ **and** the bonded RDKit
    mol) or the full set of generation knobs (``optimizer``, ``ensemble_size``,
    ``ff_preset``, ``seed``, ...); those lower-level knobs may evolve between
    releases.

    Parses the OIN string with :class:`OINParser` and generates a conformer via
    the MetalloGen backend (dummy-metal + RDKit ``CoordMap`` embed + constrained
    MMFF/UFF cleanup, refined by the ``optimizer``).
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
        """Initialize the generator with a parser and the MetalloGen backend.

        The MetalloGen backend (``oinsmiles.generator3d``) uses a dummy-metal +
        CoordMap embed refined by the ``optimizer`` (default ``"xtb"`` standard
        g-xTB; pass ``"ff"``/``None`` for the FF-only path, or
        ``"mace-omol-0-extra-large-1024"`` / ``"mace-omol25"`` for higher
        accuracy). The MACE models require ``mace-torch`` + the model weights and
        fail loudly if unavailable -- use ``optimizer="ff"`` in environments
        without them.

        ``engine`` is retained for backward compatibility and must be
        ``"metallogen"`` (the only supported backend).

        ``seed`` is the base ETKDG seed for the metallogen engine; a fixed value
        (default 42) makes generation reproducible.
        """
        if engine != "metallogen":
            raise ValueError(f"Unknown engine {engine!r}; expected 'metallogen'")

        self.parser = OINParser()
        self.engine = engine

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

    def generate(self, oin_string: str) -> GeneratedStructure:
        """Convert an OIN-SMILES string to a 3D structure.

        For just the XYZ string, prefer :meth:`SMILESToXYZ.convert`, which wraps
        this. Call this method (or :meth:`SMILESToXYZ.generate`) when you also
        need the bonded RDKit mol (``.mol``).

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
        parsed = self.parser.parse(oin_string)
        return self.adapter.generate(parsed)
