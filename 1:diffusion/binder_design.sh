#!/bin/bash

# Configs
OUTPUT_BASE_DIR="PATH/TO/OUTPUT/DIR"
RFDIFFUSION_SCRIPT="/opt/RFdiffusion-old/scripts/run_inference.py"
INPUT_PDB="PATH/TO/INPUT/PDB" 
NUM_DESIGNS=10
NOISE_SCALE_CA=0
NOISE_SCALE_FRAME=0

# Binder design paramaters
CONTIGS="[A17-132/0 40-120]"
HOTSPOT_RES="[A56,A115,A123]"

# Output directory
OUTPUT_DIR="$OUTPUT_BASE_DIR/binder_design"
mkdir  -p "$OUTPUT_DIRE"

echo "Running binder design for target: $INPUT_PDB"

python "$RFDIFFUSION_SCRIPT" \
    inference.output_prefix="$OUTPUT_DIR" \
    inference.input_pdb="$INPUT_PDB" \
    "contigmap.contigs=$CONTIGS" \
    "ppi.hotspot_res=$HOTSPOT_RES" \
    inference.num_designs=$NUM_DESIGNS \
    denoiser.noise_scale_ca=$NOISE_SCALE_CA \
    denoiser.noise_scale_frame=$NOISE_SCALE_FRAME

echo "Completed binder design for $INPUT_PDB"
