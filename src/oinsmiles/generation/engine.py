from .oin_parser import OINParser
from .molassembler_adapter import MolassemblerAdapter, MolassemblerTimeoutError

class OIN3DGenerator:
    def __init__(self, scaling_factor: float = 1.0):
        self.parser = OINParser()
        self.adapter = MolassemblerAdapter()
        
    def generate(self, oin_string: str, extra_params: dict = None):
        """
        Orchestrates the conversion from OIN string to 3D structure using Molassembler.
        """
        parsed = self.parser.parse(oin_string)
        
        # Note: extra_params are for Architector (force fields etc).
        # Molassembler is a pure DG engine. We'll ignore them for now.
        structure = self.adapter.generate(parsed)
        return structure
