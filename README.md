# Struct2Block
**Detect functional blocking.**

If you have two different proteins that bind to the same target protein. Struct2Block is an easy tool designed for quickly evaluating the following situations:

* An antibody binds to a receptor. To what extent does it sterically hinder the ligand.
* 2 different ligands bind to the same receptor. To what degree does they compete.

## 1 Install

### 1.1 pip
```console
pip install struct2block
```

### 1.2 conda
```conda
conda install chou_uken::struct2block
mamba install chou_uken::struct2block
```


## 2 How to use

### 2.1 As a Command Line Application

To get usage instructions.
```console
struct2block --help
```

For example, you can use this command line:
```console
struct2block A.pdb B.pdb
```
It will use A.pdb as antigen-ligand complex and B.pdb as Antibody-antigen complex. Then some details and results will be printed in the console. If you don't need any output, use option `-q` or `--quiet`.
You can also output superimposed structures with following command line.
```console
struct2block A.pdb B.pdb prefix
```
Then your structures will be output as prefix_antibody.pdb and prefix_ligand.pdb.

### 2.2 As a Python API
```python
def struct2block(complex: str, anti: str, prefix: str=None, quiet: bool=False) -> float:
    """Calculate the steric clash volume (block rate) of antibody.
    Block rate = Volume of ligand occupied by antibody / Total volume of ligand
    
    Args:
        complex (str): PDB file containing Antigen-Ligand model.
        anti (str): PDB file containing Antigen-Antibody model.
        prefix (str): The file prefix you want to store the superimposed complex structures.
        quiet (bool): If true, suppress the output.

    Returns:
        blockRate (float): block rate.
    """
```

For example:
```python
from struct2block.struct2block import struct2block
br: float = struct2block(complex="antigen-ligand.pdb", anti="antigen-antibody.pdb", quiet=False)
```

## 3 How does it work
First, Struct2Block find the most similar shared chain in two PDB files as antigen. Then it calculates the space (**V_ligand**) occupied by the ligand in the antigen-ligand complex. Next, it calculates the volume of that same space that is now occupied by the antibody in the antigen-antibody complex (**V_antibody ∩ ligand**). Then the steric clash volume of ligand (called 'block rate')

block rate = V(antibody ∩ ligand) / V(ligand)
