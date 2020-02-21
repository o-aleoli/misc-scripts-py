#!/usr/bin/env python3
import pymatgen as mg
from pymatgen.io.vasp.outputs import BSVasprun, Vasprun
from pymatgen import Spin
from pymatgen.electronic_structure.plotter import BSPlotter, BSDOSPlotter, DosPlotter
import matplotlib.pyplot as plt

run = BSVasprun("vasprun.xml", parse_projected_eigen=True)
bs = run.get_band_structure("KPOINTS")
BSPlotter(bs).plot_brillouin()
