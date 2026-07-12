from .oin_parser import OINParser
from .structure import GeneratedStructure

__all__ = ["OIN3DGenerator", "GeneratedStructure"]


class OIN3DGenerator:
    """Generates 3D XYZ structures from OIN-SMILES strings.

    **Warning:** This class is an internal implementation detail of OIN-SMILES.
    Its API is subject to change without notice. Users should prefer
    SMILESToXYZ.convert() instead.

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
        )

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
        parsed = self.parser.parse(oin_string)
        return self.adapter.generate(parsed)
