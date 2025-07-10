# Struct2Block
**Detect functional block.**
</hr>
If you have two different proteins that bind to the same target protein. Struct2Block is an easy tool to use for you who want to check these things quickly:

* An antibody binds to a receptor. To what extent does it sterically hinder the ligand.
* 2 different ligand binds to the same receptor. To what degree does they compete.

## Install


## How to use




## How does it work
It calculate the space (A) occupied by the ligand (V(ligand)). Then calculate the space occupied by antibody in space A (V(antibody | ligand)). Then the steric clash volume of ligand (called 'block rate')

block rate = V(antibody | ligand) / V(ligand)
