import os
import glob
import subprocess
import argparse
from loguru import logger
from Bio import SeqIO
from pathlib import Path
import sys

def parse_arguments():
    parser = argparse.ArgumentParser(description="ProteinMPNN Sequence Design Script")
    parser.add_argument("--pdb_dir", type=str, required=True, help="Directory containing PDB files.")
    parser.add_argument("--temps", type=float, nargs='+', required=True, help="Sampling temperatures (e.g., 0.1 0.2 0.3).")
    parser.add_argument("--fixed_res", type=str, required=True, help="Fixed residues (e.g., 'C11 C30 C90').")
    parser.add_argument("--omit_AA", type=str, default=None, help="Residues to be omitted from designs (e.g., 'C').")
    parser.add_argument("--model_type", choices=["protein_mpnn", "soluble_mpnn"], default="protein_mpnn")
    parser.add_argument("--checkpoint", type=str, default=None, help="Path to model checkpoint.")
    parser.add_argument("--number_of_batches", type=int, default=2, help="Number of batches to run (default: 2).")
    parser.add_argument("--batch_size", type=int, default=10, help="Batch size for sequence design (default: 10).")
    return parser.parse_args()

def run_seq_design(pdb_dir, temps, fixed_residues, omit_aa, checkpoint_proteinmpnn, number_of_batches, batch_size):
    pdb_files = glob.glob(os.path.join(pdb_dir, "*.pdb"))
    if not pdb_files:
        logger.error("No PDB files found in the specified directory!")
        sys.exit(1)

    # Select model type and checkpoint
    if args.model_type == "protein_mpnn":
        checkpoint_flag = "--checkpoint_protein_mpnn"
        checkpoint = args.checkpoint or "./model_params/proteinmpnn_v_48_020.pt"
    else:
        checkpoint_flag = "--checkpoint_soluble_mpnn"
        checkpoint = args.checkpoint or "./model_params/solublempnn_v_48_020.pt"

    if not Path(checkpoint).exists():
        logger.error(f"Checkpoint file not found: {checkpoint}")
        sys.exit(1)
    
    for pdb_file in pdb_files:
        pdb_name = Path(pdb_file).stem

        for temp in temps:
            out_dir = Path(f"{pdb_name}_temp_{temp}")
            out_dir.mkdir(parents=True, exist_ok=True)
            
            mpnn_cmd = [
                "python", "run.py",
                "--seed", "111",
                "--pdb_path", pdb_file,
                "--out_folder", str(out_dir),
                "--model_type", args.model_type,
                checkpoint_flag, checkpoint,
                "--number_of_batches", str(number_of_batches),
                "--batch_size", str(batch_size),
                "--temperature", str(temp),
                "--fixed_residues", fixed_residues,
            ]

            if omit_aa:
            mpnn_cmd.extend(["--omit_AA", omit_aa])

            logger.info(f"Running {args.model_type} on {pdb_name} at T={temp}")
            subprocess.run(mpnn_cmd, check=True)
    
            # Process FASTA files
            seqs_dir = out_dir / "seqs"
            logger.info(f"Checking directory for fasta files: {seqs_dir}")

            if seqs_dir.exists():
                logger.info(f"Listing files in {seqs_dir}: {os.listdir(seqs_dir)}")

            fasta_files = list(seqs_dir.glob("*.fa"))
            if not fasta_files:
                logger.warning(f"No fasta files found in {seqs_dir}")
                continue

            input_fasta = fasta_files[0]
            temp_str = str(temp).replace('.', '')
            output_fasta = out_dir / f"{pdb_name}_seqs_{temp_str}.fa"

            records = list(SeqIO.parse(input_fasta, "fasta"))
            with open(output_fasta, "w") as out_fname:
                for i, record in enumerate(records[1:], 1): # skip the first sequence
                    new_header = f"{pdb_name}_temp_{temp_str}_seq_{i}"
                    record.id = new_header
                    record.description = ""
                    SeqIO.write(record, out_fname, "fasta")

            logger.info(f"Fasta file {input_fasta.name} headers renamed and saved to {output_fasta.name}")

    logger.info("Sequence design and FASTA processing completed!")

if __name__ == "__main__":
    logger.remove()
    logger.add(sys.stdout, level="INFO", format="{time} {level} {message}")
    args = parse_arguments()
    run_seq_design(args)
