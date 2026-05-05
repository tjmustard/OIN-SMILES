from .oin_parser import OINParser
from .molassembler_adapter import MolassemblerAdapter, MolassemblerTimeoutError, GeneratedStructure

__all__ = ["OIN3DGenerator", "MolassemblerTimeoutError", "GeneratedStructure"]


class OIN3DGenerator:
    """Generates 3D XYZ structures from OIN-SMILES strings.

    Delegates parsing to OINParser and 3D generation to MolassemblerAdapter.
    """

    def __init__(
        self,
        timeout: int = 60,
        dg_strategy: str = "single",
        ensemble_size: int = 10,
    ) -> None:
        self.parser = OINParser()
        self.adapter = MolassemblerAdapter(
            timeout=timeout,
            dg_strategy=dg_strategy,
            ensemble_size=ensemble_size,
        )

    def generate(self, oin_string: str) -> GeneratedStructure:
        """Convert an OIN-SMILES string to a 3D structure.

        Parameters
        ----------
        oin_string:
            OIN-SMILES string (V3.0 inline or V2.4 sidecar format).

        Returns
        -------
        GeneratedStructure
            ``.xyz`` contains the XYZ block string.
            ``.mol`` contains an RDKit Mol with bond connectivity and a 3D
            conformer (None if connectivity could not be determined).

        Raises
        ------
        MolassemblerTimeoutError
            If conformer generation exceeds the configured timeout.
        RuntimeError
            If Molassembler fails to generate a conformer.
        """
        parsed = self.parser.parse(oin_string)
        return self.adapter.generate(parsed)
