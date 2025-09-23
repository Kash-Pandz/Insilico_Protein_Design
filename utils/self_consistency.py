import argparse
import subprocess
from pathlib import Path
import pandas as pd
import MDAnalysis as mda
from MDAnalysis.core.groups import AtomGroup
from MDAnalysis.analysis import rms


# Atom selection presets
ATOM_SELECTIONS = {
    'calpha': 'name CA',
    'backbone': 'backbone',
    'heavy': 'not name H*',
    'all': 'protein'
}


def select_atom_group (
    universe: mda.Universe, 
    atom_type: str,
    chain_id: str = None, 
    custom: str = None
) -> AtomGroup:
    """
    Select atom type for RMSD calculation.

    Args:
        universe: MDAnalysis Universe object.
        atom_type: 'calpha', 'backbone', 'heavy', 'all'.
        chain_id: Optional chain ID(s).
        custom: Optional custom selection string.

    Returns:
        AtomGroup object of selected atom type.
    """
    base_sel = ATOM_SELECTIONS.get(atom_type)
    if not base_sel:
        raise ValueError(f"Invalid atom atype: {atom_type}")
      
    selection = f"{base_sel} and ({custom})" if custom else base_sel
    if chain_id:
        selection += f" and segid {chain_id}"
      
    return universe.select_atoms(selection)


def get_rmsd(
    ref_atoms: AtomGroup, 
    sel_atoms: AtomGroup
) -> float:
    """
    Calculate rmsd between two AtomGroups.

    Args:
        ref_atoms: Reference AtomGroup.
        sel_atoms: Selection AtomGroup.

    Returns:
        RMSD as float (Å).
    """
    if len(ref_atoms) != len(sel_atoms):
        raise ValueError("Atom count mismatch between ref and sel.")
  
    return rms.rmsd(
        ref_atoms.positions, 
        sel_atoms.positions, 
        center=True, 
        superposition=True
    )


def get_tmscore(ref_pdb: Path, sel_pdb: Path, tmscore_bin: Path) -> float:
    """
    Calculate the TM-Score beteen two PDB files.

    Args:
        ref_pdb: Path to reference PDB.
        sel_pdb: Path to selection PDB.
        tmscore_bin: Path to TM-Score executable.

    Returns:
        TM-Score as float.
    """
    
    tmscore_cmd = subprocess.run(
        [str(tmscore_bin), str(ref_pdb), str(sel_pdb)],
        capture_output=True, 
        text=True
    )
    
    for line in tmscore_cmd.stdout.splitlines():
        if line.startswith("TM-score"):
            return float(line.split("=")[1].strip().split()[0])
    return 0.0


def main():
    parser = argparse.ArgumentParser(
        description="Protein design self-consistency analysis."
    )
    parser.add_argument("ref_pdb", type=Path, help="Reference PDB file (Path)")
    parser.add_argument("sel_dir", type=Path, help="Directory with selection PDBs (Path)")
    parser.add_argument("tmscore_bin", type=Path, help="TM-score executable (Path)")
    parser.add_argument(
        "--atom_type",
        choices=ATOM_SELECTIONS.keys(),
        default="calpha",
        help="Atom type for RMSD (str)"
    )
    parser.add_argument(
        "--chain",
        type=str,
        help="Required chain ID(s) for RMSD, e.g. 'A' or 'B_C'"
    )
    parser.add_argument(
        "--custom",
        type=str,
        help="Optional custom residue selection string"
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Output CSV file (Path or None)"
    )

    args = parser.parse_args()

    ref_universe = mda.Universe(args.ref_pdb)
    sel_pdbs = sorted(args.sel_dir.glob("*.pdb"))
    if not sel_pdbs:
        raise ValueError(f"No PDB files found in {args.sel_dir}")

    values = []

    for sel_pdb in sel_pdbs:
        sel_universe = mda.Universe(sel_pdb)
        row = {"ref_pdb": args.ref_pdb.stem, "sel_pdb": sel_pdb.stem)

        # Overall RMSD
        if args.chain:
            chains = args.chain.split("_")
            for c in chains:
                
          

        # Overall RMSD
        overall_rmsd = get_rmsd(
            select_atom_group(ref_universe, args.atom_type),
            select_atom_group(sel_universe, args.atom_type)
        )

        # Motif RMSD
    




  
