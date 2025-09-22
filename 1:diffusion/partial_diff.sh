#!/bin/bash

# Configs
PARTIAL_T_VALUES=(1 5 7 10 12 15 20)
INPUT_DIR="/PATH/TO/INPUT/DIR"
OUTPUT_BASE_DIR="/PATH/TO/OUTPUT/DIR" 
RFDIFFUSION_SCRIPT="/scripts/run_inference.py"
CONTIGS="[A1-1/93-93/A95-95/63-63/A159-159/28-28/A188-188/40-40]"
NUM_DESIGNS=1

# Loop over all PDB files in input directory
for INPUT_PDB in "$INPUT_DIR"/*.pdb; do
    PDB_NAME=$(basename "$INPUT_PDB" .pdb)
    for T in "${PARTIAL_T_VALUES[@]}"; do
        OUTPUT_DIR="$OUTPUT_BASE_DIR/$PDB_NAME/pdb_pT${T}"
        mkdir -p "$OUTPUT_DIR"
        echo "Running RFdiffusion for $PDB_NAME, partial_T=${T}..."
        
        python "$RFDIFFUSION_SCRIPT" \
            inference.output_prefix="$OUTPUT_DIR" \
            inference.input_pdb="$INPUT_PDB" \
            "contigmap.contigs=$CONTIGS" \
            inference.num_designs=$NUM_DESIGNS \
            diffuser.partial_T=$T
        
        echo "Completed $PDB_NAME, partial_T=${T}"
    done
done

echo "Partial diffusion runs completed!"
