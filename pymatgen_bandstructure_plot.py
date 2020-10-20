#!/bin/env python3
"""Alexandre Olivieri (olivieri.alexandre0@gmail.com)

Simple plot generation from the vasprun.xml and KPOINTS files."""

from pymatgen.io.vasp import BSVasprun
from pymatgen.electronic_structure.plotter import BSPlotterProjected
from matplotlib import pyplot as plt

VRUN = BSVasprun('vasprun.xml', parse_projected_eigen=True)
BS = VRUN.get_band_structure(kpoints_filename='KPOINTS', line_mode=True)
BANDPLOTTER = BSPlotterProjected(BS)
plot = BANDPLOTTER.get_projected_plots_dots_patom_pmorb(
        {'O':['p'], 'Nb':['dyz', 'dx2', 'dz2']},
        {'O':['all'], 'Nb':['all']},
        sum_atoms={'O':['all'], 'Nb':['all']},
        sum_morbs={'O':['p']},
        num_column=3
)
plot.savefig('bandstructure.pdf', dpi=1000)
