"""
struct2block.py: Detects functional block from aligned PDB files.
Authors: Zhang Yujian
"""

import typer
from typing_extensions import Annotated
from rich import print
import numpy as np
from biotite.structure.io import pdb
from biotite import structure as struc
from biotite import sequence as seq

def struct2block(complex: Annotated[str, typer.Argument(help="The PDB file contains 1 model: Antigen-Ligand.")], \
                 anti: Annotated[str, typer.Argument(help="The PDB file contains 1 model: Antigen-Antibody.")], \
                 resolution: Annotated[float, typer.Argument(help="Resolution selected for calculation")] = 1.) -> float:
    """Calculate the block rate of antibody. block rate = V(ligand occupied by Antibody) / V(ligand)
    
    Args:
        complex(str): PDB file containing Antigen-Ligand model.
        anti(str): PDB file containing Antigen-Antibody model.
        resolution(float): Voxel size.

    Returns:
        blockRate(float): block rate.
    """

    # Load Antigen-Ligand complex
    ligComplex_file: pdb.PDBFile = pdb.PDBFile.read(complex)
    ligComplex: structure.AtomArray = ligComplex_file.get_structure(model=1)
    
    # Load Antigen-Antibody complex
    antiComplex_file: pdb.PDBFile = pdb.PDBFile.read(anti)
    antiComplex: struc.AtomArray = antiComplex_file.get_structure(model=1)
    
    # Find the same chains
    
    # Align and find antigen

    # Create voxel of ligand
    # Mark voxel occupied by antibody
    # Calculate
    
    print(type(superimposed))
    print(type(transformation))


if __name__ == "__main__":
    typer.run(struct2block)
    
