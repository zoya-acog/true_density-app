import sys
import json
import base64
import io
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw

def predict_density(smiles):
    mol_no_H = Chem.MolFromSmiles(smiles)
    if mol_no_H is None:
        raise ValueError("Invalid SMILES string")

    mol_H = Chem.AddHs(mol_no_H)
    total_volume = 0.0
    processed_atoms = set()

    benzene = Chem.MolFromSmarts("c1ccccc1")
    for match in mol_H.GetSubstructMatches(benzene):
        total_volume += 75.2
        processed_atoms.update(match)

    naphthalene = Chem.MolFromSmarts("c1ccc2ccccc2c1")
    for match in mol_H.GetSubstructMatches(naphthalene):
        total_volume += 123.7
        processed_atoms.update(match)

    cooh = Chem.MolFromSmarts('[C](=O)[OH]')
    total_volume += len(mol_H.GetSubstructMatches(cooh)) * -2.6

    conh = Chem.MolFromSmarts('[C](=O)[NH]')
    conh2 = Chem.MolFromSmarts('[C](=O)[NH2]')
    total_volume += (len(mol_H.GetSubstructMatches(conh)) + len(mol_H.GetSubstructMatches(conh2))) * -2.8

    nh_n = Chem.MolFromSmarts('[NH][N]')
    total_volume += len(mol_H.GetSubstructMatches(nh_n)) * -0.3

    ring_info = mol_H.GetRingInfo()
    for ring in ring_info.AtomRings():
        if len(ring) in [5, 6]:
            if not any(mol_H.GetAtomWithIdx(i).GetIsAromatic() for i in ring):
                if not any(i in processed_atoms for i in ring):
                    total_volume += -3.0
                    processed_atoms.update(ring)

    for atom in mol_H.GetAtoms():
        idx = atom.GetIdx()
        if idx in processed_atoms:
            continue

        symbol = atom.GetSymbol()
        bonds = atom.GetBonds()
        bond_types = [b.GetBondType().name for b in bonds]

        if symbol == 'H':
            total_volume += 6.9
        elif symbol == 'C':
            s = bond_types.count('SINGLE')
            d = bond_types.count('DOUBLE')
            t = bond_types.count('TRIPLE')
            if d == 2:
                total_volume += 15.3
            elif d == 1 and s == 2:
                total_volume += 13.7
            elif t == 1 and s == 1:
                total_volume += 15.3
            elif s == 4:
                total_volume += 11.0
            else:
                total_volume += 15.3
        elif symbol == 'O':
            if 'DOUBLE' in bond_types:
                total_volume += 14.0
            else:
                total_volume += 9.2
        elif symbol == 'N':
            if 'TRIPLE' in bond_types or len(bonds) == 1:
                total_volume += 16.0
            elif 'DOUBLE' in bond_types or atom.GetIsAromatic():
                total_volume += 12.8
            else:
                total_volume += 7.2
        elif symbol == 'S':
            total_volume += 23.8
        elif symbol == 'F':
            total_volume += 12.8
        elif symbol == 'Cl':
            total_volume += 26.7
        elif symbol == 'Br':
            total_volume += 33.0
        elif symbol == 'I':
            total_volume += 45.0
        else:
            raise ValueError(f"Unsupported element: {symbol}")

    mol_weight = Descriptors.MolWt(mol_H)
    density = (mol_weight * 1.660) / total_volume
    error_density = (mol_weight * 1.645) / total_volume
    volume = total_volume



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
        print("Usage: python immirzi_perini.py 'SMILES_STRING'")
        sys.exit(1)

    smiles = sys.argv[1]
    try:
        result = predict_density(smiles)
        print(json.dumps([result]))
    except Exception as e:
        print(json.dumps([{"smiles": smiles, "error": str(e)}]))