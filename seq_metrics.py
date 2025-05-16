import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Align import substitution_matrices


def calculate_sequence_properties(seq):
    """Function to calculate sequence properties."""
    analysis = ProteinAnalysis(seq)
    aa_content = analysis.get_amino_acids_percent()
    
    charge_fraction = sum(aa_content.get(aa, 0) 
                          for aa in ["K", "R", "D", "E"]) 
    
    return {
    "length": len(seq),
    "molecular_weight": analysis.molecular_weight(),
    "pI": analysis.isoelectric_point(),
    "hydrophobicity": analysis.gravy(),
    "charge_fraction": charge_fraction,
    }



def calculate_global_similarity(ref_seq, seq):
    """Calculates global sequence similarity relative to reference sequence."""

    assert len(ref_seq) == len(seq), "Sequences must be of equal length."
    
    # Exclude gaps in alignment
    matches = sum(1 for r, a in zip(ref_seq, seq) if r == a and r != "-" and a != "-")
    total_positions = sum(1 for r, a in zip(ref_seq, seq) if r != "-" and a != "-")
    
    if total_positions == 0:
        return 0.0  
    
    identity = matches / total_positions
    return identity



def calculate_active_site_similarity(ref_seq, target_seq, matrix, gap_open_penalty=-10, gap_extend_penalty=-1):
    """Calculate active site similarity score using BLOSUM62 matrix with gap penalties."""
    score = 0
    in_gap = False 

    for ref_aa, target_aa in zip(ref_seq, target_seq):
        if ref_aa == '-' or target_aa == '-': 
            if not in_gap:
                score += gap_open_penalty
                in_gap = True
            else: 
                score += gap_extend_penalty
        else:
            in_gap = False
            try:
                score += matrix[ref_aa, target_aa]
            except KeyError:
                score += -4
    return score



def map_active_site_positions(alignment, reference_positions):
    """Map the active site positions in the multiple sequence alignment."""
    # Reference sequence (first sequence in alignment)
    ref_seq = alignment[0].seq 
    mapped_positions = []
    ref_index = 0 
    for alignment_index, aa in enumerate(ref_seq):
        if aa != '-':
            ref_index += 1
        if ref_index in reference_positions:
            mapped_positions.append(alignment_index + 1)
    return mapped_positions



def get_seq_metrics(fasta_file, active_site_positions=None, gap_open_penalty=-10, gap_extend_penalty=-1):
    """Gather all sequence metrics."""

    # Read fasta file
    fasta = list(SeqIO.parse(fasta_file, "fasta"))
    matrix = substitution_matrices.load("BLOSUM62")
    ref_seq = fasta[0].seq

    # Map active site positions
    mapped_active_site_positions = []
    if active_site_positions:
        mapped_active_site_positions = map_active_site_positions(fasta, active_site_positions)
    
    results = []

    for seq_record in fasta:
        seq = str(seq_record.seq)
        seq_id = seq_record.id 

        # Sequence properties
        seq_props = calculate_sequence_properties(seq)

        # Global similarity relative to reference sequence
        global_similarity = calculate_global_similarity(ref_seq, seq)

        # Active site similarity given the active site positions
        if active_site_positions:
            ref_active_site_seq = ''.join(ref_seq[pos - 1] for pos in mapped_active_site_positions)
            target_active_site_seq = ''.join(seq[pos - 1] for pos in mapped_active_site_positions)

            active_site_score = calculate_active_site_similarity(ref_active_site_seq, target_active_site_seq, matrix,
                                                                            gap_open_penalty, gap_extend_penalty)
        else:
            active_site_score = None
            target_active_site_seq = None
        
    results.append({
        "identifier": seq_id,
        "sequence": seq,
        "global_similarity": global_similarity,
        "active_site_score": active_site_score,
        "active_site_seq": target_active_site_seq,
        **seq_props
    })
    
    return pd.DataFrame(results)
