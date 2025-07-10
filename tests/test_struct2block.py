import pytest
from src.struct2block.struct2block import struct2block
import os


def test_struct2block():
    ligComplexPath: str = os.path.join(os.path.dirname(__file__), "data", "spike_s1-ace2_ecd.pdb")
    antiComplexPath: str = os.path.join(os.path.dirname(__file__), "data", "9cff-antibody.pdb")
    br: float = struct2block(ligComplexPath, antiComplexPath)
    assert br >= 0., "Error"

