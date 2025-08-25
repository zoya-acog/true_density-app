import sys
import json
import argparse
import csv
import requests
import base64
import io
from rdkit import Chem
from rdkit.Chem import Descriptors, Draw
from girolami import predict_density as girolami_predict
from girolami_group import predict_density as girolami_group_predict
from girolami_name import name_to_smiles, predict_density as girolami_name_predict
from immirzi_perini import predict_density as immirzi_predict
from immirzi_group import predict_density as immirzi_group_predict

def main():
    parser = argparse.ArgumentParser(description="CLI for density prediction using Girolami or Immirzi methods.")
    parser.add_argument("command", help="Command in format method-inputtype, e.g., girolami-smiles")
    parser.add_argument("input", help="SMILES string, compound name, or file path")
    parser.add_argument("--json", action="store_true", help="Output in JSON format")
    args = parser.parse_args()

    try:
        method, input_type = args.command.split("-")
        if method not in ["girolami", "immirzi"]:
            raise ValueError(f"Unsupported method: {method}")
        if input_type not in ["smiles", "name", "file"]:
            raise ValueError(f"Unsupported input type: {input_type}")

        results = []
        if input_type == "file":
            with open(args.input, newline="") as csvfile:
                reader = csv.DictReader(csvfile)
                for row in reader:
                    compound = row.get("compound", "N/A")
                    smiles = row["smiles"]
                    try:
                        if method == "girolami":
                            result = girolami_group_predict(smiles, compound)
                        else:
                            result = immirzi_group_predict(smiles, compound)
                        results.append(result)
                    except Exception as e:
                        results.append({"compound": compound, "smiles": smiles, "error": str(e)})
        else:
            if input_type == "name":
                compound = args.input
                smiles = name_to_smiles(compound)
                if method == "girolami":
                    result = girolami_name_predict(smiles, compound)
                else:
                    result = immirzi_predict(smiles)
                    result["compound"] = compound
            else:
                smiles = args.input
                if method == "girolami":
                    result = girolami_predict(smiles)
                else:
                    result = immirzi_predict(smiles)
                result["compound"] = "N/A"
            results.append(result)

        if args.json:
            print(json.dumps(results))
        else:
            for res in results:
                if "error" in res:
                    print(f"Error for {res.get('compound', res['smiles'])}: {res['error']}")
                else:
                    print(f"Compound: {res.get('compound', 'N/A')}")
                    print(f"SMILES: {res['smiles']}")
                    print(f"Density: {res['density']:.3f} g/cm³")
                    print(f"Error Density: {res['error_density']:.3f} g/cm³")
                    print(f"Volume: {res['volume']:.2f} Å³")
                    print("---")

    except Exception as e:
        if args.json:
            print(json.dumps([{"error": str(e)}]))
        else:
            print(f"Error: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()