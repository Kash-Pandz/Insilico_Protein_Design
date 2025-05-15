# mpnn_seq_design

This repository provides tools for applying ProteinMPNN and LigandMPNN, which are graph neural network-based models for protein sequence design. LigandMPNN extends the capabilities of ProteinMPNN by incorporating ligand atoms into the message-passing framework. 


# Steps:

1) Starting apo or holo protein structure (crystal or alphafold)
2) Define design space - consider fixing functionally important residues (e.g. catalytic residues and evolutionary conserved residues). Catalytic residues can be determined from literature or using tools such as m-CSA. Evolutionary conserved residues can be determined using multiple sequence alignments (e.g. HHBlits)
3) Run ProteinMPNN (protein_mpnn.py) or LigandMPNN (ligand_mpnn.py) with desired parameters such as fixed residues, biasing towards/away from particular residues and omitting residues. A broad sampling temperature range can be used (e.g. 0.1-1.0) to generate increased sequence diversity.
4) Optional: In order to diversify the input protein backbone, partial diffusion protocol in RFDiffusion/RFDiffusion-All-Atom can be used. 
5) Selection criteria:
   - Structure-based metrics (pLDDT, RMSD and TMScore) from  structure predictions made with colabfold.sh
   - Sequence-based metrics (sequence properties and pseudo-log-likelihoods) 

