import sys
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from collections import Counter
import base64
import io
import json

# Girolami atomic volume mapping
volume_map = {
    "H": 1,
    "Li": 2, "Be": 2, "B": 2, "C": 2, "N": 2, "O": 2, "F": 2,
    "Na": 4, "Mg": 4, "Al": 4, "Si": 4, "P": 4, "S": 4, "Cl": 4,
    "K": 5, "Ca": 5, "Ga": 5, "Ge": 5, "As": 5, "Se": 5, "Br": 5,
    "Rb": 7.5, "Sr": 7.5, "In": 7.5, "Sn": 7.5, "Sb": 7.5, "Te": 7.5, "I": 7.5,
    "Cs": 9, "Ba": 9, "Tl": 9, "Pb": 9, "Bi": 9
}

def predict_density(smiles):
    mol_no_H = Chem.MolFromSmiles(smiles)
    if mol_no_H is None:
        raise ValueError(f"Invalid SMILES string: {smiles}")

    mol_H = Chem.AddHs(mol_no_H)
    atom_counts = Counter(atom.GetSymbol() for atom in mol_H.GetAtoms())

    molar_volume = 0.0
    for atom, count in atom_counts.items():
        if atom not in volume_map:
            raise ValueError(f"Element '{atom}' not supported by Girolami method")
        molar_volume += volume_map[atom] * count

    molecular_weight = Descriptors.MolWt(mol_H)
    density = molecular_weight / (3.85 * molar_volume)
    error_density = molecular_weight / (3.74 * molar_volume)
    volume = 3.85 * molar_volume  # Adding volume for consistency with requirements

    img = Draw.MolToImage(mol_no_H, size=(400, 400))
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

    return {
        "smiles": smiles,
        "density": density,
        "error_density": error_density,
        "volume": volume,
        "image": img_base64
    }

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python girolami.py 'SMILES_STRING'")
        sys.exit(1)

    smiles = sys.argv[1]
    try:
        result = predict_density(smiles)
        print(json.dumps([result]))
    except Exception as e:
        print(json.dumps([{"smiles": smiles, "error": str(e)}]))