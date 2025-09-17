import pandas as pd
from Bio import AlignIO, SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from Bio.Align import substitution_matrices


def calculate_sequence_properties(seq: str) -> dict:
    """Function to calculate sequence properties."""
    analysis = ProteinAnalysis(seq)
    aa_content = analysis.get_amino_acids_percent()
    
    charge_fraction = sum(
        aa_content.get(aa, 0) for aa in ["K", "R", "D", "E"]
    ) 
    
    return {
        "length": len(seq),
        "molecular_weight": analysis.molecular_weight(),
        "pI": analysis.isoelectric_point(),
        "hydrophobicity": analysis.gravy(),
        "charge_fraction": charge_fraction
    }


def calculate_global_similarity(ref_seq: str, seq: str) -> float:
    """Calculates global sequence similarity relative to reference sequence."""

    assert len(ref_seq) == len(seq), "Sequences must be of equal length."
    
    matches = sum(
        1 for r, a in zip(ref_seq, seq) if r == a and r != "-" and a != "-"
    )
    total_positions = sum(
        1 for r, a in zip(ref_seq, seq) if r != "-" and a != "-"
    )
    
    if total_positions == 0:
        return 0.0  
    
    return matches / total_positions


def calculate_active_site_similarity(
    ref_seq: str,
    target_seq: str,
    matrix,
    gap_open_penalty=-10,
    gap_extend_penalty=-1,
) -> int:
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


def map_active_site_positions(alignment, reference_positions: list) -> list:
    """Map the active site positions in the multiple sequence alignment."""
    ref_seq = alignment[0].seq 
    mapped_positions = []
    ref_index = 0 
    
    for alignment_index, aa in enumerate(ref_seq):
        if aa != '-':
            ref_index += 1
        if ref_index in reference_positions:
            mapped_positions.append(alignment_index + 1)
            
    return mapped_positions


def find_longest_repeat(seq: str, k: int) -> int:
    """Find the longest run of consecutive repeats of a motif of length k."""

    if k <= 0 or k > len(seq):
        return 0

    longest = 1
    current = 1
    prev = seq[0:k]

    for i in range(k, len(seq), k):
        motif = seq[i:i+k]
        if motif == prev:
            current += 1
            longest = max(longest, current)
        else:
            current = 1
            prev = motif

    return longest
        

def get_seq_metrics(
    fasta_file: str,
    active_site_positions: list = None, 
    gap_open_penalty: int = -10, 
    gap_extend_penalty: int = -1,
    max_repeat_k: int = 5
) -> pd.DataFrame:
    """Gather all sequence metrics."""

    alignment = AlignIO.read(fasta_file, "fasta")
    matrix = substitution_matrices.load("BLOSUM62")
    ref_seq = alignment[0].seq

    mapped_active_site_positions = []
    if active_site_positions:
        mapped_active_site_positions = map_active_site_positions(
            alignment, active_site_positions
        )
    
    results = []

    for seq_record in alignment:
        seq = str(seq_record.seq)
        seq_id = seq_record.id 

        # Remove gaps
        clean_seq = seq.replace("-", "")
        seqs_props = calculate_sequence_properties(clean_seq)

        # Global similarity relative to reference sequence
        global_similarity = calculate_global_similarity(ref_seq, seq)

        # Active site similarity given the active site positions
        if active_site_positions:
            ref_active_site_seq = "".join(
                ref_seq[pos - 1] for pos in mapped_active_site_positions
            )
            target_active_site_seq = "".join(
                seq[pos - 1] for pos in mapped_active_site_positions
            )

            active_site_score = calculate_active_site_similarity(
                ref_active_site_seq, 
                target_active_site_seq, 
                matrix,
                gap_open_penalty, 
                gap_extend_penalty
            )
        else:
            active_site_score = None
            target_active_site_seq = None

        # Longest repeats for k=1 to k=max_repeat
        repeats = {
            f"longest_repeat_{k}": find_longest_repeat(clean_seq, k)
            for k in range(1, max_repeat_k + 1)
        }
        
        results.append(
            {
                "identifier": seq_id,
                "sequence": seq,
                "global_similarity": global_similarity,
                "active_site_score": active_site_score,
                "active_site_seq": target_active_site_seq,
                **seq_props,
                **repeats
            }
        )
    
    return pd.DataFrame(results)
