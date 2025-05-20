# mpnn_seq_design

This repository provides tools for applying ProteinMPNN and LigandMPNN, which are graph neural network-based models for protein sequence design. LigandMPNN extends the capabilities of ProteinMPNN by incorporating ligand atoms into the message-passing framework. 


# Steps:

1) Starting apo or holo protein structure (crystal or alphafold)
2) Define design space - consider fixing functionally important residues (e.g. catalytic residues and evolutionary conserved residues). Catalytic residues can be determined from literature or using tools such as m-CSA. Evolutionary conserved residues can be determined using multiple sequence alignments (e.g. HHBlits)
3) Run ProteinMPNN (protein_mpnn.py) or LigandMPNN (ligand_mpnn.py) with desired parameters such as fixed residues, biasing towards/away from particular residues and omitting residues. A broad sampling temperature range can be used (e.g. 0.1-1.0) to generate increased sequence diversity. 
4) Selection criteria for laboratory testing:
   - Structure-based metrics (pLDDT, RMSD and TMScore) from  structure predictions made with colabfold.sh (alternatively use the colabfold notebook provided from https://github.com/sokrypton/ColabFold).
   - Sequence-based metrics:
      - sequence properties: global similarity, active site similarity, hydropathy and charged fraction
      - pseudo-log-likelihoods (ESM-2): How biologically plausible a sequence is? Does it align well to natural sequences? A higher PLL (closer to zero, since log-probs are negative) suggests the sequence is more natural-like. 


# N.B.  

- LigandMPNN can be used to redesign active/binding site residues alone using the --redesigned_residues flag. This will only redesign specified residues and fix everything else. The ligand_mpnn.py provided is for sequence redesign fixing a specified number of residues.
- LigandMPNN outputs protein structures of the generated sequences with packed side-chains. These structures can be optimised using Pyrosetta or OpenMM energy minimisation. Both will require parameterisation of the ligand. It is recommended to "refold" the generated sequences with an independent folding software e.g. AF2, ESMFold and ColabFold. 


# Resources

ProteinMPNN and LigandMPNN weights and inference script can be found:

https://github.com/dauparas/LigandMPNN?tab=readme-ov-file

Local Colab installation instructions can be found:

https://github.com/sokrypton/ColabFold

https://github.com/YoshitakaMo/localcolabfold

HHBlits package can be found:

https://github.com/soedinglab/hh-suite

Alternatively jobs can be submitted to https://toolkit.tuebingen.mpg.de/tools/hhblits

Biopython: https://biopython.org/

MDAnalysis: https://www.mdanalysis.org/

TMScore: https://zhanggroup.org/TM-score/TMscore.cpp

transformers:  https://huggingface.co/docs/transformers/installation
