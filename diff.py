import subprocess
from pathlib import Path
import argparse
from loguru import logger


def run_diff(
    name: str,
    outdir: str,
    contigs: str,
    base_path: str,
    pdb: str = None,
    iterations: int = 50,
    num_designs: int = 10,
    uncond_design: bool = False,
    guide_scale: float = 1.0,
    guide_potentials: str = None,
    deterministic: bool = False,
    partial_diff: list[int] = None,
):

  # Output
  results_dir = Path(outdir) / name / "results"
  results_dir.mkdir(parents=True, exist_ok=True)
  logger.add(results_dir / f"{name}.log", level="INFO")

  logger.info(f"Starting job: {name}")

  # partial diff
  steps = partial_diff if partial_diff else [None]

  for step in steps:
      subdir = results_dir if step is None else results_dir / f"partial_{step}"
      subdir.mkdir(parents=True, exist_ok=True)
      prefix = subdir / (f"{name}_partial_{step}" if step else name)

      # Add general diffusion parameters
      args = {
          "inference.output_prefix": str(prefix),
          "inference.num_designs": num_designs,
          "contigmap.contigs": f"[{contigs}]"
      }

      # Add diffusion or partial diffusion steps
      if step:
        args["diffuser.partial_T"] = step
      else:
        args["diffuser.T"] = iterations
        
      # Add potentials to unconditional design
      if uncond_design:
          args["potential.guide_scale"] = guide_scale
          if guide_potentials:
              args["potentials.guide_potentials"] = f'["{guide_potential}"]'
            
      # Deterministic run      
      if deterministic:
          args["inference.deterministic"] = True

      # Option of input pdb or unconditional design  
      if pdb and not uncond_design:
          args["inference.input_pdb"] = str(Path(pdb).resolve())

      # Add run_inference script
      cmd = f"cd {Path(base_path).expanduser()} && python run_inference.py " \
            + " ".join(f"{k}={v}" for k, v in args.items())

      logger.infor(f"Running: {cmd}")
      subprocess.run(cmd, shell=True, check=True)
      logger.success(f"Step {step or 'full'} finished for {name}")
    
if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--base_path", required=True, help="Path to rfDiffusion")
    p.add_argument("--name", required=True, help="Base job name")
    p.add_argument("--outdir", required=True, help="Output directory")
    p.add_argument("--contigs", required=True, help="Contigs definition")
    p.add_argument("--pdbs", nargs="*", help="One or more input PDB files")
    p.add_argument("--iterations", type=int, default=50)
    p.add_argument("--num_designs", type=int, default=10)
    p.add_argument("--uncond", action="store_true", help="Unconditional design (no PDB needed)")
    p.add_argument("--guide_scale", type=float, default=1.0)
    p.add_argument("--guide_potentials")
    p.add_argument("--deterministic", action="store_true")
    p.add_argument("--partial", help="Comma-separated steps, e.g. 25,50,75")
    args = p.parse_args()

    partial_list = [int(x) for x in args.partial.split(",")] if args.partial else None

    if not args.pdbs and not args.uncond:
        raise ValueError("Provide --pdbs or use --uncond")

    pdbs = args.pdbs or [None]  # run unconditional if no pdbs
    for pdb in pdbs:
        job_name = args.name if not pdb else f"{args.name}_{Path(pdb).stem}"
        run_diff(
            name=job_name,
            outdir=args.outdir,
            contigs=args.contigs,
            base_path=args.base_path,
            pdb=pdb,
            iterations=args.iterations,
            num_designs=args.num_designs,
            uncond=args.uncond,
            guide_scale=args.guide_scale,
            guide_potentials=args.guide_potentials,
            deterministic=args.deterministic,
            partial=partial_list,
        )


