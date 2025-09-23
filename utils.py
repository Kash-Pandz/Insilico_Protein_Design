import os
import MDAnalysis as mda

def get_res_string(res_input):
    """ 
    Returns a space-separated string of residues from a string or txt file.

    """
    if not res_input:
        return ""

    if os.path.isfile(res_input):
        with open(res_input, 'r') as f:
            lines = [line.strip() for line in f if line.strip()]
        res_str = ' '.join(lines)
    else:
        res_str = res_input

    # Create residue ranges
    expanded = []
    for part in res_str.split():
        if '-' in part:
            chain = part[0]
            start, end = map(int, part[1:].split('-'))
            expanded.extend([f"{chain}{i}" for i in range(start, end + 1)])
        else:
            expanded.append(part)
          
    return ' ' .join(expanded)


def get_contig_map(pdb_file: str, fixed_residues: list[int], chain_id="A": str) -> str:
    """
    Creates a contig map for RFDiffusion from PDB file.

    Args:
        pdb_file: Path to input PDB file.
        res_sel: Residue IDs (PDB numbering).
        chain_id: Chain ID to use for contig string.

    Returns:
        Contig map

    """
    # Load PDB file in MDAnalyis universe object
    universe = mda.Universe(pdb_file)

    # Extract sorted residue IDs
    resids = sorted[res.resid for res in universe.residues])

    contig_map = []

    # Preceding block before first selected residue
    first_res = resid[0]
    first_fixed = fixed_residues[0]
    if first_fixed > first_res:
        seg_len = first_fixed - first_res
        contig_map.append(f"{seg_len}-{seg_len}")

    # In-between fixed residue and gaps between them
    for i, fixed_res in enumerate(fixed_residues):
        if i > 0:
            prev_res = fixed_residues[i-1]
            gap_len = fixed_res - prev_res - 1
            if gap_len > 0:
                contig_map.append(f"{gap_len}-{gap_len}")
        contig_map.append(f"{chain_id}(fixed_res}-{fixed_res}")

   # Block after last fixed residue
   last_fixed = fixed_residues[-1]
   last_res = resids[-1]
   if last_res > last_fixed:
       seg_len = last_res - last_fixed
       contig_map.append(f"{seg_len}-{seg_len}")

   return "[" + "/".join(contig_map) + "]"
