#!/usr/bin/env python3
from pymatgen.io.vasp.outputs import BSVasprun

raw = BSVasprun("vasprun.xml")
bandstructure = raw.get_band_structure("KPOINTS")
bandstructure.get_band_gap()
