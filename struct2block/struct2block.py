"""
struct2block.py: Detects functional block from aligned PDB files.
Authors: Zhang Yujian
Date: Jul, 2025
"""

import typer
from typing_extensions import Annotated
from rich import print
import numpy as np
from biotite.structure.io import pdb
from biotite import structure as struc
from biotite import sequence as seq


def alignStruct(structA: struc.AtomArray, structB: struc.AtomArray) -> tuple[struc.AtomArray, struc.AffineTransformation]:
    """Align (superimpose) 2 structures with a shared peptide.

    Arg:
        structA (struc.AtomArray): Reference structure.
        structB (struc.AtomArray): Mobile structure.

    Returns:
        (superimposed structB, transformation)
    """

    # Find the same chains
    structAChains: list[str] = np.unique(structA.chain_id)
    structBChains: list[str] = np.unique(structB.chain_id)
    # Chains valid.
    if (len(structAChains) != 2):
        raise (Exception(f"{len(structAChains)} detected in Antigen-Ligand complex, 2 expected."))
    if (len(structBChains) != 2):
        raise (Exception(f"{len(structBChains)} detected in Antigen-Antibody complex, 2 expected."))
    # Get each chain sequence and find anchor chains
    anchorA: str = "mark"
    anchorB: str = "mark"
    for candiAnchorA in structAChains:
        for candiAnchorB in structBChains:
            if (struc.to_sequence(structA[structA.chain_id == candiAnchorA])[0][0] == \
                struc.to_sequence(structB[structB.chain_id == candiAnchorB])[0][0]):
                anchorA = candiAnchorA
                anchorB = candiAnchorB
                break
    if ((anchorA == "mark") or (anchorB == "mark")):
        raise (Exception("No shared chain in two complex."))
    # Calculate transformation and superimposition.
    _, transformation = struc.superimpose(structA[structA.chain_id == anchorA], structB[structB.chain_id == anchorB])
    structC: struc.AtomArray = transformation.apply(structB)
    return (structC, transformation)



def struct2block(complex: Annotated[str, typer.Argument(help="The PDB file contains 1 model: Antigen-Ligand.")], \
                 anti: Annotated[str, typer.Argument(help="The PDB file contains 1 model: Antigen-Antibody.")], \
                 resolution: Annotated[float, typer.Argument(help="Resolution selected for calculation")] = 1.) -> float:
    """Calculate the block rate of antibody. block rate = V(ligand occupied by Antibody) / V(ligand)
    
    Args:
        complex (str): PDB file containing Antigen-Ligand model.
        anti (str): PDB file containing Antigen-Antibody model.
        resolution (float): Voxel size.

    Returns:
        blockRate (float): block rate.
    """

    # Load Antigen-Ligand complex
    ligComplex_file: pdb.PDBFile = pdb.PDBFile.read(complex)
    ligComplex: structure.AtomArray = ligComplex_file.get_structure(model=1)
    # Load Antigen-Antibody complex
    antiComplex_file: pdb.PDBFile = pdb.PDBFile.read(anti)
    antiComplex: struc.AtomArray = antiComplex_file.get_structure(model=1)
    # Alignment
    superimposedAntiComplex, transformation = alignStruct(ligComplex, antiComplex)
    # Find ligand and antibody
    notLigandId: int = -1
    notAntibodyId: int = -1
    complexChains: list[str] = np.unique(ligComplex.chain_id)
    antiChains: list[str] = np.unique(antiComplex.chain_id)
    for candiLigandId in range(len(complexChains)):
        for candiAntiId in range(len(antiChains)):
            if (struc.to_sequence(ligComplex[ligComplex.chain_id == complexChains[candiLigandId]])[0][0] == \
                struc.to_sequence(antiComplex[antiComplex.chain_id == antiChains[candiAntiId]])[0][0]):
                notLigandId = candiLigandId
                notAntibodyId = candiAntiId
                break
    ligandId: str = np.delete(complexChains, notLigandId)[0]
    antibodyId: str = np.delete(antiChains, notAntibodyId)[0]
    print(ligandId, antibodyId)
    # Create voxel of ligand
    # Mark voxel occupied by antibody
    # Calculate    
    # print(type(superimposed))
    # print(type(transformation))


if __name__ == "__main__":
    typer.run(struct2block)
    
