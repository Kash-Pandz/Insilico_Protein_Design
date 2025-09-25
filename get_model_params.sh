#!/bin/bash

e.g. bash get_model_params.sh "./model_params"

mkdir -p $1

# ProteinMPNN Weights
wget -q https://files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_020.pt -O $1"/proteinmpnn_v_48_020.pt"

# SolubleMPNN Weights
wget -q https://files.ipd.uw.edu/pub/ligandmpnn/solublempnn_v_48_020.pt -O $1"/solublempnn_v_48_020.pt"

# LigandMPNN Weights
wget -q https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_010_25.pt -O $1"/ligandmpnn_v_32_010_25.pt"
wget -q https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_sc_v_32_002_16.pt -O $1"/ligandmpnn_sc_v_32_002_16.pt"
