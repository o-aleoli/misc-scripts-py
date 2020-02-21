#!/usr/bin/env python3
from pymatgen import Structure
from pymatgen.io.vasp.inputs import Kpoints

density = 150

structure = Structure.from_file( "POSCAR" )
kmesh = Kpoints.automatic_density_by_vol( structure, density )
kmesh.write_file("KPOINTS")
