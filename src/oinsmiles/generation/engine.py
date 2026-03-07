from .oin_parser import OINParser
from .molassembler_adapter import MolassemblerAdapter, MolassemblerTimeoutError

__all__ = ["OIN3DGenerator", "MolassemblerTimeoutError"]


class OIN3DGenerator:
    """Generates 3D XYZ structures from OIN-SMILES strings.

    Delegates parsing to OINParser and 3D generation to MolassemblerAdapter.
    """

    def __init__(self, timeout: int = 60) -> None:
        self.parser = OINParser()
        self.adapter = MolassemblerAdapter(timeout=timeout)

    def generate(self, oin_string: str) -> str:
        """Convert an OIN-SMILES string to a 3D XYZ block.

        Parameters
        ----------
        oin_string:
            OIN-SMILES string (V3.0 inline or V2.4 sidecar format).

        Returns
        -------
        str
            XYZ block string of the generated 3D conformer.

        Raises
        ------
        MolassemblerTimeoutError
            If conformer generation exceeds the configured timeout.
        RuntimeError
            If Molassembler fails to generate a conformer.
        """
        parsed = self.parser.parse(oin_string)
        return self.adapter.generate(parsed)
