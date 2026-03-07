
from rdkit import Chem

def test_inline_roundtrip():
    print("Test 5: Inline Handler Roundtrip (Sanitize=False)")
    
    # Input from OINSanitizer Debug Output
    # c1cc[c]c(-c2cccc[n]2)c1
    s_in = "c1cc[c]c(-c2cccc[n]2)c1"
    print(f"Input: {s_in}")
    
    # Simulate OINInlineHandler
    m = Chem.MolFromSmiles(s_in, sanitize=False)
    if not m:
        print("Failed to parse")
        return
        
    # Apply Map Num (Simulate tagging)
    # Atom 2 is [c].
    a = m.GetAtomWithIdx(3) # c1cc[c] -> 0,1,2,3? Brackets preserve atomic status?
    # Index check:
    # c:0, c:1, c:2, [c]:3, c:4, c:5 ?
    # Let's check indices
    for at in m.GetAtoms():
        # print(f"{at.GetIdx()}: {at.GetSymbol()} IsAro={at.GetIsAromatic()}")
        if at.GetSymbol() == 'C' and at.GetNumExplicitHs() == 0 and at.GetNoImplicit():
             print(f"Target Atom {at.GetIdx()} found.")
             at.SetAtomMapNum(1000)
             
    s_out = Chem.MolToSmiles(m, canonical=False)
    print(f"Output: {s_out}")
    
    if "[cH" in s_out:
        print("FAIL: Got [cH]")
    elif "[c" in s_out:
        print("PASS: Got [c]")
        
        # Check replace logic
        import re
        def replace_map(match):
            content = match.group(1)
            # is_pure_organic = re.fullmatch(r"^(N|P|S|c|n|o|p|s)$", content)
            # if is_pure_organic: return f"{content}{{0}}"
            return f"[{content}]{{0}}"
            
        tagged = re.sub(r"\[([^:\]]+):(\d+)\]", replace_map, s_out)
        print(f"Tagged: {tagged}")

if __name__ == "__main__":
    test_inline_roundtrip()
