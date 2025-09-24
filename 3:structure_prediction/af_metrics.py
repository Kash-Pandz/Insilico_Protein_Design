import json
import numpy as np
import re
from pathlib import Path
import argparse
from typing import Dict, Any


def calculate_metrics(json_dir: str) -> Dict[str, Dict[str, Any]]:
    """
    Calulates AlphaFold prediction metrics.

    Args:
        json_dir (str): Path to the folder containing score JSON files.

    Returns:
        Dictionary of metrics per model.
    """
    json_dir_path = Path(json_dir)
    metrics_dict: Dict[str, Dict[str, Any]] = {}

    for filepath in json_dir_path.glob("*_scores_rank_*.json"):
        with open(filepath, "r") as f:
            data = json.load(f)

        # Clean model name (remove "_scores_rank_*")
        model_name = re.sub(r"_scores_rank_.*$", "", filepath.stem)

        # Extract metrics
        plddt = np.array(data.get("plddt", []), dtype=float)
        avg_plddt = float(np.mean(plddt)) if plddt.size > 0 else None
        iptm = data.get("iptm")
        ptm = data.get("ptm")

        pae = np.array(data.get("predicted_aligned_error", []), dtype=float)
        mean_pae = float(np.mean(pae)) if pae.size > 0 else None

        metrics_dict[model_name] = {
            "plddt": avg_plddt,
            "iptm": iptm,
            "ptm": ptm,
            "mean_pae": mean_pae,
        }
        
    return metrics_dict


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Get model metrics (mean pLDDT, ipTM, pTM, mean PAE) from JSON files."
    )
    parser.add_argument("folder", type=str, help="Path to folder with JSON files.")

    args = parser.parse_args()

    metrics_dict = calculate_metrics(args.folder)

    # Save metrics to JSON
    folder_path = Path(args.folder)
    output_file = folder_path / f"{folder_path.name}_metrics.json"
    with open(output_file, "w") as f:
        json.dump(metrics_dict, f, indent=2)
    print(f"Metrics saved to {output_file}")


if __name__ == "__main__":
    main()


    





