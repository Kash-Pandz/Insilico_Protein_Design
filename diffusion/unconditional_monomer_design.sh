#!/bin/bash

# Configs
OUTPUT_BASE_DIR="PATH_TO_OUTPUT_DIR"
RFDIFFUSION_SCRIPT="/opt/RFdiffusion-old/scripts/run_inference.py"
NUM_DESIGNS=10
CONTIGS="[100-200]" 
GUIDE_POTENTIAL='["type:monomer_ROG,weight:1,min_dist:5"]'
GUIDE_SCALE=2
GUIDE_DECAY="quadratic"

# Output directory
OUTPUT_DIR="$OUTPUT_BASE_DIR/monomer_design"
mkdir -p "$OUTPUT_DIR"

echo "Running unconditional monomer generation with contigs=$CONTIGS..."

python "$RFDIFFUSION_SCRIPT" \
    inference.output_prefix="$OUTPUT_DIR" \
    "contigmap.contigs=$CONTIGS" \
    inference.num_designs=$NUM_DESIGNS
    "potentials.guiding_potentials=$GUIDE_POTENTIAL" \
    potentials.guide_scale=$GUIDE_SCALE \
    potentials.guide_decay="$GUIDE_DECAY"

echo "Completed unconditional monomer generation!"
