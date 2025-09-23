import MDAnalysis as mda
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import capped_distance
from MDAnalysis.core.groups import AtomGroup, ResidueGroup
import numpy as np
from typing import Dict, Tuple



def get_residues(pdb_file: str, selection: str) -> ResidueGroup:
    """ Obtain residues of the two binding partners in a pdb file."""

    universe = mda.Universe(pdb_file)
    residues = universe.select_atoms(selection).residues
    if len(residues) == 0:
        raise ValueError(f"Residue selection '{selection}' returned 0 atoms")
    return residues
  

def get_interface_residues(g1: ResidueGroup, g2: ResidueGroup, cutoff: float = 8.0) -> Dict[str, Tuple[Residue]]:
    """ Identify interface residues between two selections."""
    contacts, _ = capped_distance(
        g1.atoms.positions,
        g2.atoms.positions,
        max_cutoff=cutoff,
        return_distances=True
    )
    left, right = contacts.T

    interface_1 = tuple(
        sorted({g1.atoms[l].residue for l in left}, key=lambda g: g.resid
    )
    interface_2 = tuple(
        sorted({g2.atoms[r].residue for r in right}, key=lambda g: g.resid
    )

    return {"group1": interface_1, "group2": interface_2}
    
