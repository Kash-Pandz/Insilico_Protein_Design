# Protein_Design

This repository showcases an in-silico protein design pipeline using RFDiffusion + MPNN Sequence Design + ColabFold. It also includes an evaluation pipeline of various design tasks.


 ProteinMPNN and LigandMPNN, which are graph neural network-based models for protein sequence design. LigandMPNN extends the capabilities of ProteinMPNN by incorporating ligand atoms into the message-passing framework. 

# Overview

- `ligand_mpnn.py`: Sequence Design (with ligand)
- `protein_mpnn.py`: Sequence Design (without ligand)
- `colabfold.sh`: Starting protein prediction/forward folding with localcolab fold
- `plddt.py`: predicted local distance difference test (pLDDT) for all structure predictions
- `self_consistency.py`: Self-consistency metrics TMScore and RMSD for designs
- `seq_metrics.py`: The global similarity, active site similarity, sequence properties
- `pseudo_log_likelihood.py`: log likelihood for designed sequences used for design filtering

# Setup

### LigandMPNN/ProteinMPNN
- Installation instructions and inference scripts: https://github.com/dauparas/LigandMPNN
- Model weights can be downloaded:
  
`wget -q https://files.ipd.uw.edu/pub/ligandmpnn/proteinmpnn_v_48_020.pt -O $1"/proteinmpnn_v_48_020.pt"`

`wget -q https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_010_25.pt -O $1"/ligandmpnn_v_32_010_25.pt"`

`wget -q https://files.ipd.uw.edu/pub/ligandmpnn/ligandmpnn_v_32_020_16.pt -O $1"/ligandmpnn_v_32_020_16.pt"`


### Colab Fold
- Follow installation instructions:
  -  https://github.com/sokrypton/ColabFold
  -  https://github.com/YoshitakaMo/localcolabfold

# Steps:

1) Starting apo or holo protein structure (crystal or alphafold)
2) Define design space - consider fixing functionally important residues (e.g. catalytic residues and evolutionary conserved residues). Catalytic residues can be determined from literature or using tools such as M-CSA. Evolutionary conserved residues can be determined using multiple sequence alignments (e.g. HHBlits)
3) Run ProteinMPNN (`protein_mpnn.py`) or LigandMPNN (`ligand_mpnn.py`) with desired parameters such as fixed residues, biasing towards/away from particular residues and omitting residues. A broad sampling temperature range can be used (τ ∈ {0.1, 0.2, ..., 1}) to generate sequence diversity. 
4) Selection criteria for laboratory testing:
   - Structure-based metrics (pLDDT, RMSD and TMScore) from  structure predictions made with `colabfold.sh` (alternatively use the colabfold notebook provided from https://github.com/sokrypton/ColabFold).
   - Sequence-based metrics:
      - sequence properties: global similarity, active site similarity, hydropathy and charged fraction
      - pseudo-log-likelihoods (ESM-2): How biologically plausible a sequence is? Does it align well to natural sequences? A higher PLL (closer to zero, since log-probs are negative) suggests the sequence is more natural-like. 


# N.B.  

- LigandMPNN can be used to redesign active/binding site residues alone using the --redesigned_residues flag. This will only redesign specified residues and fix everything else. The ligand_mpnn.py provided is for sequence redesign fixing a specified number of residues.
- LigandMPNN outputs protein structures of the generated sequences with packed side-chains. These structures can be optimised using Pyrosetta or OpenMM energy minimisation. Both will require parameterisation of the ligand. It is recommended to "refold" the generated sequences with an independent folding software e.g. AF2, ESMFold and ColabFold. 
