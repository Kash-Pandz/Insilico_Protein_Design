#!/bin/bash

OUTPUT_DIR=gen_seq_preds/
FASTA_DIR=gen_seq_fastas/
mkdir -p $OUTPUT_DIR

colabfold_batch $FASTA_DIR $OUTPUT_DIR --num-models 1 --model-type alphafold2_multimer_v3
