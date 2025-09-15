import argparse
import subprocess
from pathlib import Path
from typing import Optional
from loguru import logger
from itertools import product


def parse_args():
    parser = argparse.ArgumentParser(description="Run RFDiffusion")
    parser.add_argument("--base_path", required=True, help="Path to RFDiffusion")
    parser.add_argument("--name", required=True, help="Job Name")
    parser.add_argument("--outdir", required=True, help="Output Directory")
    parser.add_argument("--contigs", required=True, help="Contigs Map")
    parser.add_argument("--pdb_path", nargs="+", help="Path to input PDB files")
    parser.add_argument("--iterations", type=int, default=50)
    parser.add_argument("--num_designs", type=int, default=10)
    parser.add_argument("--uncond_design", action="store_true", help="Unconditional design (no PDB input required))
    parser.add_argument("--guide_scale", type=float, default=1.0)
    parser.add_argument("--guide_potentials")
    parser.add_argument("--deterministic", action="store_true")
    parser.add_argument(
        "--partial_diff", 
        nargs="+", 
        type=int, 
        help="List of denoising steps for partial diffusion (e.g. --partial_diff 5 10 12")
    )
    return parser.parse_args()
                        
    
def run_diff(args):
    # Atleast one input should be provided
    pdbs = args.pdbs or [None]
    if not pdbs and not arg.uncond_design:
        raise ValueError("Provide --pdbs or use --uncond_design")

    # Results directory 
    results_dir = Path(args.outdir) / args.name / "results"
    results_dir.mkdir(parents=True, exist_ok=True)

    steps = args.partial_diff or [None]

    for pdb, step in product(pdbs, steps):
        job_name = args.name
        if pdb:
            job_name += f"_{pdb.stem}"
        if step:
            job_name += f"_partial_{step}"
            
        logger.add(results_dir / f"{job_name}.log", level="INFO")
        logger.info(f"Starting job: {job_name})

        # Subdirectory for partial diffusion steps
        subdir = results_dir if step is None else results_dir / f"partial_{step}"
        subdir.mkdir(parents=True, exist_ok=True)

        prefix = subdir / job_name

        cmd_args = {
            "inference.output_prefix": str(prefix),
            "inference.num_designs": args.num_designs
        }

        # Iterations based on partial diffusion
        cmd_args["diffuser.partial_T" if step else "diffuser.T"] = step or args.iterations

        # Unconditional design settings
        if args.uncond_design:
            cmd_args["potential.guide_scale"] = args.guide_scale
            if args.guide_potentials:
                cmd_args["potentials.guide_potentials"] = f'["{args.guide_potentials}"]'

        # Deterministic mode
        if args.deterministic:
            cmd_args["inference.deterministic"] = True
            
        # Add PDB input if not unconditional
        if pdb and not args.uncond_design:
            cmd_args["inference.input_pdb"] = str(pdb.resolve())

        cmd = ["python", "run_inference.py"] + [f"{k}={v}" for k, v in cmd_args.items()]

        logger.info("Running: {}", " ".join(cmd))
        subprocess.run(cmd, cwd=args.base_path, check=True)
        logger.success(f"Step {step or 'full'} finished for {job_name}")


def main():
    args = parse_args()
    run_diff(args)


if __name__ == "__main__":
    main()

            
                



    

    



