#!/usr/bin/env python3
from pymatgen import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

structure = Structure.from_file("POSCAR")
analyze = SpacegroupAnalyzer(structure)
analyze.get_space_group_symbol()