"""
struct2block.py: Detects functional block from aligned PDB files.
Authors: Zhang Yujian
Date: Jul, 2025
"""

import typer
from typing import Optional
from typing_extensions import Annotated
from rich import print
import numpy as np
from biotite.structure.io import pdb
from biotite import structure as struc
from biotite import sequence as seq
from biotite.sequence import align
import os
from rich.progress import Progress


def alignStruct(structA: struc.AtomArray, structB: struc.AtomArray) -> tuple[struc.AtomArray, struc.AffineTransformation, str, str]:
    """Align (superimpose) 2 structures with a shared peptide.

    Arg:
        structA (struc.AtomArray): Reference structure.
        structB (struc.AtomArray): Mobile structure.

    Returns:
        (superimposed structB, transformation, anchorCom, anchorAnti)
    """

    # Find the same chains
    structAChains: list[str] = np.unique(structA.chain_id)
    structBChains: list[str] = np.unique(structB.chain_id)
    # Get each chain sequence and find anchor chains
    anchorA: str = "mark"
    anchorB: str = "mark"
    bestIdent: float = 0.
    alph: seq.LetterAlphabet = seq.ProteinSequence.alphabet
    scores: np.ndarray = np.identity(len(alph), dtype=int)
    matrix: align.SubstitutionMatrix = align.SubstitutionMatrix(alph, alph, scores)
    for candiAnchorA in structAChains:
        for candiAnchorB in structBChains:
            alignment: align.Alignment = align.align_optimal(
                struc.to_sequence(structA[structA.chain_id == candiAnchorA])[0][0],
                struc.to_sequence(structB[structB.chain_id == candiAnchorB])[0][0],
                matrix
            )[0]
            ident: float = align.get_sequence_identity(alignment)
            if (ident > bestIdent):
                anchorA = candiAnchorA
                anchorB = candiAnchorB
                bestIdent = ident
                bestAlignment = alignment
    if (bestIdent < 0.6):
        print(f"Warning. Low identity of shared antigen: {np.round(bestIdent, 4)*100}%.")
    alignCode: np.ndarray = align.get_codes(bestAlignment)
    anchorMask: np.ndarray = ((matrix.score_matrix()[alignCode[0], alignCode[1]]) & (alignCode != -1).all(axis=0))
    anchorMask = np.array(anchorMask, dtype=bool)
    superimpositionAnchor: np.ndarray = bestAlignment.trace[anchorMask]
    # Calculate transformation and superimposition.
    sharedA: struc.AtomArray = structA[structA.chain_id == anchorA]
    sharedB: struc.AtomArray = structB[structB.chain_id == anchorB]
    _, transformation = struc.superimpose(sharedA[sharedA.atom_name == "CA"][superimpositionAnchor[:, 0]], \
                                          sharedB[sharedB.atom_name == "CA"][superimpositionAnchor[:, 1]])
    structC: struc.AtomArray = transformation.apply(structB)
    return (structC, transformation, anchorA, anchorB)



def struct2block(complex: Annotated[str, typer.Argument(help="The PDB file contains 1 model: Antigen-Ligand.")], \
                 anti: Annotated[str, typer.Argument(help="The PDB file contains 1 model: Antigen-Antibody.")], \
                 prefix: Annotated[Optional[str], typer.Argument(help="The file you want to store the superimposed structures in. Antigen-Ligand complex will be in {prefix}_ligand.pdb. Antigen-Antibody will be in {prefix}_antibody.pdb")] = None, \
                 quiet: Annotated[bool, typer.Option("--quiet", "-q", help="Suppress the output.")] = False) -> float:
    """Calculate the steric clash volume (block rate) of antibody.  = V(ligand occupied by Antibody) / V(ligand)
    
    Args:
        complex (str): PDB file containing Antigen-Ligand model.
        anti (str): PDB file containing Antigen-Antibody model.
        prefix (str): The file prefix you want to store the superimposed complex structures.
        quiet (bool): If true, suppress the output.

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
    superimposedAntiComplex, transformation, notLigandId, notAntibodyId = alignStruct(ligComplex, antiComplex)
    # Find ligand and antibody
    complexChains: list[str] = np.unique(ligComplex.chain_id)
    antiChains: list[str] = np.unique(antiComplex.chain_id)
    ligandId: str = np.delete(complexChains, np.where(complexChains == notLigandId))
    antibodyId: str = np.delete(antiChains, np.where(antiChains == notAntibodyId))
    if (not quiet):
        print(f"Found Ligand: Chain {ligandId} in '{complex}'")
        print(f"Found Antibody: Chain {antibodyId} in '{anti}'")
    # Create voxel of ligand
    ligandStruct: struc.AtomArray = ligComplex[ligComplex.chain_id != notLigandId]
    antibodyStruct: struc.AtomArray = superimposedAntiComplex[superimposedAntiComplex.chain_id != "C"]
    # Remvoe hydrogen
    ligandStruct = ligandStruct[ligandStruct.element != "H"]
    antibodyStruct = antibodyStruct[antibodyStruct.element != "H"]
    xMin: np.float32 = np.floor(np.min(ligandStruct.coord[:,0])) - 3.
    xMax: np.float32 = np.ceil(np.max(ligandStruct.coord[:,0])) + 3.
    yMin: np.float32 = np.floor(np.min(ligandStruct.coord[:,1])) - 3.
    yMax: np.float32 = np.ceil(np.max(ligandStruct.coord[:,1])) + 3.
    zMin: np.float32 = np.floor(np.min(ligandStruct.coord[:,2])) - 3.
    zMax: np.float32 = np.ceil(np.max(ligandStruct.coord[:,2])) + 3.
    xSize: np.int64 = np.int64(xMax - xMin + 1)
    ySize: np.int64 = np.int64(yMax - yMin + 1)
    zSize: np.int64 = np.int64(zMax - zMin + 1)
    shape: tuple[np.int64, np.int64, np.int64] = (xSize, ySize, zSize)
    ligVoxels: np.ndarray = np.zeros(shape, dtype=bool)
    atomRadii: dict[str, np.float32] = {}
    if (not quiet):
        print(f"Created {ligVoxels.size} voxels totally.")
    with Progress() as progress:
        task1 = progress.add_task("Marking atoms of ligand...", total = len(ligandStruct))
        with open(os.path.join(os.path.dirname(__file__), "data", "van_der_Waals_Radii.csv"), mode="r") as f:
            lines: list[str] = list(map(lambda x:x.strip("\n"), f.readlines()[1:]))
            for line in lines:
                atomRadii[line.split(",")[0]] = np.float32(line.split(",")[1])
            for atom in ligandStruct:
                zoneCenter: np.ndarray = np.array([np.int64(atom.coord[0]-xMin), np.int64(atom.coord[1]-yMin), np.int64(atom.coord[2]-zMin)])
                for zoneX in range(zoneCenter[0]-3, zoneCenter[0]+4):
                    for zoneY in range(zoneCenter[1]-3, zoneCenter[1]+4):
                        for zoneZ in range(zoneCenter[2]-3, zoneCenter[2]+4):
                            dist: np.float32 = struc.distance(np.array([zoneX, zoneY, zoneZ]), atom.coord - np.array([xMin, yMin, zMin]))
                            if (dist <= atomRadii[atom.element]):
                                ligVoxels[zoneX, zoneY, zoneZ] = True
                progress.update(task1, advance=1)
    
        ligandVoxelNum: np.int64 = np.count_nonzero(ligVoxels)
    
        antiVoxels: np.ndarray = np.zeros(shape, dtype=bool)
        task2 = progress.add_task("Marking atoms of antibody...", total = len(antibodyStruct))
        with open(os.path.join(os.path.dirname(__file__), "data", "van_der_Waals_Radii.csv"), mode="r") as f:
            lines: list[str] = list(map(lambda x:x.strip("\n"), f.readlines()[1:]))
            for line in lines:
                atomRadii[line.split(",")[0]] = np.float32(line.split(",")[1])
            for atom in antibodyStruct:
                zoneCenter: np.ndarray = np.array([np.int64(atom.coord[0]-xMin), np.int64(atom.coord[1]-yMin), np.int64(atom.coord[2]-zMin)])
                for zoneX in range(zoneCenter[0]-3, zoneCenter[0]+4):
                    for zoneY in range(zoneCenter[1]-3, zoneCenter[1]+4):
                        for zoneZ in range(zoneCenter[2]-3, zoneCenter[2]+4):
                            dist: np.float32 = struc.distance(np.array([zoneX, zoneY, zoneZ]), atom.coord - np.array([xMin, yMin, zMin]))
                            if (dist > atomRadii[atom.element]):
                                try:
                                    antiVoxels[zoneX, zoneY, zoneZ] = True
                                except IndexError:
                                    pass
                                else:
                                    pass
                progress.update(task2, advance=1)
    ligVoxels = ligVoxels & antiVoxels    
    # Calculate
    blockRate: np.float32 = np.count_nonzero(ligVoxels) / ligandVoxelNum
    if (not quiet):
        print(f"Block rate: {np.round(blockRate * 100, 2)}%. (Ligand volume {ligandVoxelNum} Å^3, antibody occupies {np.count_nonzero(ligVoxels)} Å^3)")
    if (prefix):
        file1 = pdb.PDBFile()
        file2 = pdb.PDBFile()
        file1.set_structure(ligComplex)
        file2.set_structure(superimposedAntiComplex)
        file1.write(prefix + "_ligand.pdb")
        file2.write(prefix + "_antibody.pdb")
    return (blockRate)

def main():
    typer.run(struct2block)


if __name__ == "__main__":
    main()
    
