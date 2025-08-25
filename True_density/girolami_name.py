import sys
import json
import requests
from rdkit import Chem
from rdkit.Chem import Descriptors
from collections import Counter
import base64
import io
from rdkit.Chem import Draw

volume_map = {
    "H": 1,
    "Li": 2, "Be": 2, "B": 2, "C": 2, "N": 2, "O": 2, "F": 2,
    "Na": 4, "Mg": 4, "Al": 4, "Si": 4, "P": 4, "S": 4, "Cl": 4,
    "K": 5, "Ca": 5, "Ga": 5, "Ge": 5, "As": 5, "Se": 5, "Br": 5,
    "Rb": 7.5, "Sr": 7.5, "In": 7.5, "Sn": 7.5, "Sb": 7.5, "Te": 7.5, "I": 7.5,
    "Cs": 9, "Ba": 9, "Tl": 9, "Pb": 9, "Bi": 9
}

def name_to_smiles(name):
    url = f"https://cactus.nci.nih.gov/chemical/structure/{name}/smiles"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text.strip()
    else:
        raise ValueError(f"No SMILES found for: {name}")

def predict_density(smiles, compound):
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
    volume = 3.85 * molar_volume

    img = Draw.MolToImage(mol_no_H, size=(400, 400))
    buffered = io.BytesIO()
    img.save(buffered, format="PNG")
    img_base64 = base64.b64encode(buffered.getvalue()).decode("utf-8")

    return {
        "compound": compound,
        "smiles": smiles,
        "density": density,
        "error_density": error_density,
        "volume": volume,
        "image": img_base64
    }

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python girolami_name.py 'compound_name'")
        sys.exit(1)

    compound_name = sys.argv[1]
    try:
        smiles = name_to_smiles(compound_name)
        result = predict_density(smiles, compound_name)
        print(json.dumps([result]))
    except Exception as e:
        print(json.dumps([{"compound": compound_name, "error": str(e)}]))