#!/bin/bash
module unload python
module load cuda/12.4.0
module load anaconda

# Activate conda environment
conda activate colabfold-conda

# Input fasta directory and output directory where predicted structures will be stored
OUTPUT_DIR=gen_seq_preds/
FASTA_DIR=gen_seq_fastas/
mkdir -p $OUTPUT_DIR

# Run command for colabfold_batch with input fasta 
colabfold_batch $FASTA_DIR $OUTPUT_DIR \
  --num-models 1 \
  --model-type alphafold2_multimer_v3 \
  --num-recycle 3
