from .oin_parser import OINParser
from .architector_adapter import ArchitectorAdapter
from .wrapper import ArchitectorWrapper

class OIN3DGenerator:
    def __init__(self, scaling_factor: float = 1.0):
        self.parser = OINParser()
        self.adapter = ArchitectorAdapter(scaling_factor=scaling_factor)
        self.wrapper = ArchitectorWrapper()
        
    def generate(self, oin_string: str, extra_params: dict = None):
        parsed = self.parser.parse(oin_string)
        arch_args = self.adapter.convert(parsed)
        
        if extra_params:
            if "parameters" not in arch_args:
                arch_args["parameters"] = {}
            arch_args["parameters"].update(extra_params)
            
        structure = self.wrapper.run(**arch_args)
        return structure
