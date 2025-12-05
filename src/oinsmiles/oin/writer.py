from typing import List, Tuple, Dict

class OINWriter:
    def __init__(self):
        pass

    def write(self, smiles: str, coords: List[Tuple[int, float, float, float]], 
              dative_bonds: List[Tuple[int, int]] = None, 
              haptic_groups: List[str] = None, 
              geometry_label: str = None) -> str:
        """
        Constructs the OIN string.
        """
        tags = []

        # w tag (Coordinates)
        if coords:
            # Format: idx:x,y,z;idx:x,y,z
            coord_strings = [f"{idx}:{x:.4g},{y:.4g},{z:.4g}" for idx, x, y, z in coords]
            tags.append(f"w:{';'.join(coord_strings)}")

        # d tag (Dative Bonds)
        if dative_bonds:
            bond_strings = [f"{src}.{tgt}" for src, tgt in dative_bonds]
            tags.append(f"d:{';'.join(bond_strings)}")

        # m tag (Haptic Groups)
        if haptic_groups:
            tags.append(f"m:{';'.join(haptic_groups)}")

        # g tag (Geometry Label)
        if geometry_label:
            tags.append(f"g:{geometry_label}")

        if not tags:
            return smiles

        # Join tags with pipe | as per OIN v1.2
        oin_block = "|".join(tags)
        return f"{smiles} |{oin_block}|"
