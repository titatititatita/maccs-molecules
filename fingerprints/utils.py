from rdkit import Chem
from rdkit.Chem import MACCSkeys


def make_fingerprints(smiles_string: str) -> str:
    smiles = smiles_string.split('\n')
    print(smiles)
    fingerprint_flag = False
    result = []
    for s in smiles:
        if check_smile(s):
            fingerprint_flag = True
            molecule = Chem.MolFromSmiles(s)
            fingerprint = MACCSkeys._pyGenMACCSKeys(molecule)
            result.append((s + ": " + key_to_string(fingerprint)))
    if not fingerprint_flag:
        return "All smiles were written incorrectly! Check the input."
    return '\n'.join([line for line in result])


def key_to_string(fingerprint) -> str:
    fingerprint_tuple = tuple(fingerprint.GetOnBits())
    result = [0] * 166
    for bit in fingerprint_tuple:
        result[bit] = 1
    result = ''.join([str(bit) for bit in result])
    return result


def check_smile(smile: str) -> bool:
    molecule = Chem.MolFromSmiles(smile)
    if molecule is None:
        return False
    return True
