#!/usr/bin/env python3
'''Alexandre Olivieri (olivieri.alexandre0@gmail.com)

Post-processing script for VASP output

Reads plain text .dat files in INPUT_FILE and generates a projected density 
of states (pDOS) plot with desired dimensions and DPI.
It uses the LaTeX backend of matplotlib to renderize the eulervm fonts'''

import matplotlib as mpl
import numpy as np
from matplotlib.backends.backend_pgf import PdfPages
from matplotlib import rc
import matplotlib.pyplot as plt
import matplotlib.style

INPUT_FILE = [
    'NbSurf_eg', 'NbSurf_t2g',
    'NbSubsurf_eg', 'NbSubsurf_t2g',
    'NbBulk_eg', 'NbBulk_t2g',
    'total'
]
OUTPUT_FILE = 'pdosNblayer'
XLIM = (-7, 10)
YLIM = (-5, 5)
DATA = []
for i in INPUT_FILE:
    DATA.append(np.load(i+'.npz'))

MM2IN = 25.4
WIDTH = (0.5)*(210.0 - 20.0 - 30.0)/MM2IN
HEIGHT = (0.5)*(297.0 - 20.0 - 30.0)/MM2IN

"""
INPUT_DATA is a pickled, compressed file with the following numpy arrays:
    * pdos, an np.array, with dimensions (ispin, numpoints), where ispin=2 and numpoints=301 (default from VASP)
    * eigen_energy, an np.array, with dimensions (numpoints)
    * eigen_adjust, a float
    * e_fermi, a float
"""
with PdfPages(OUTPUT_FILE+'.pdf') as OUTPUT:
    plt.rcParams.update({
        'font.family': 'serif',
        'figure.figsize': [WIDTH, HEIGHT],
        'savefig.dpi': 300,
        'text.usetex': True,
        'pgf.rcfonts': False,
        'pgf.texsystem': 'pdflatex',
        'pgf.preamble': [r'\usepackage{mathpazo,eulervm}\usepackage[utf8x]{inputenc}'],
        'legend.loc': 'upper right',
        'legend.frameon': False,
        'legend.fancybox': False,
    })
    FIG, (AX0, AX1, AX2) = plt.subplots(
        nrows=3,
        ncols=1,
        sharex=True,
        sharey=True,
        constrained_layout=True
    )
    x_axis = DATA[0]['eigen_energy'] + DATA[0]['eigen_adjust']
# Surface
    AX0.plot(
        x_axis,  DATA[0]['pdos'][0, ...],
        x_axis, -DATA[0]['pdos'][1, ...],
        c='0.3333',
        label=r'$e_g$'
    )

    AX0.fill_between(
        x_axis, DATA[0]['pdos'][0, ...], -DATA[0]['pdos'][1, ...],
        facecolor='none',
        edgecolor='0.3333',
        hatch='---'
    )

    AX0.plot(
        x_axis,  DATA[1]['pdos'][0, ...],
        x_axis, -DATA[1]['pdos'][1, ...],
        c='0.6667',
        label=r'$t_{2g}$'
    )

    AX0.fill_between(
        x_axis, DATA[1]['pdos'][0, ...], -DATA[1]['pdos'][1, ...],
        facecolor='none',
        edgecolor='0.6667',
        hatch='|||'
    )

    AX0.plot(
        x_axis,  DATA[-1]['pdos'][0, ...],
        x_axis, -DATA[-1]['pdos'][1, ...],
        c='black',
        label='Total'
    )

    AX0.plot(
        [DATA[0]['e_fermi'], DATA[0]['e_fermi']],
        [np.amax(DATA[-1]['pdos'][0, ...]), -np.amax(DATA[-1]['pdos'][1, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='red',
        label=r'$\varepsilon_F$'
    )
    
    AX0.set_title('Surface')
# Subsurface
    AX1.plot(
        x_axis,  DATA[2]['pdos'][0, ...],
        x_axis, -DATA[2]['pdos'][1, ...],
        c='0.3333',
        label=r'$e_g$'
    )

    AX1.fill_between(
        x_axis, DATA[2]['pdos'][0, ...], -DATA[2]['pdos'][1, ...],
        facecolor='none',
        edgecolor='0.3333',
        hatch='---'
    )

    AX1.plot(
        x_axis,  DATA[3]['pdos'][0, ...],
        x_axis, -DATA[3]['pdos'][1, ...],
        c='0.6667',
        label=r'$t_{2g}$'
    )

    AX1.fill_between(
        x_axis, DATA[3]['pdos'][0, ...], -DATA[3]['pdos'][1, ...],
        facecolor='none',
        edgecolor='0.6667',
        hatch='|||'
    )

    AX1.plot(
        x_axis,  DATA[-1]['pdos'][0, ...],
        x_axis, -DATA[-1]['pdos'][1, ...],
        c='black',
        label='Total'
    )
    
    AX1.plot(
        [DATA[0]['e_fermi'], DATA[0]['e_fermi']],
        [np.amax(DATA[-1]['pdos'][0, ...]), -np.amax(DATA[-1]['pdos'][1, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='red',
        label=r'$\varepsilon_F$'
    )

    AX1.set_title('Subsurface')
# Bulk
    AX2.plot(
        x_axis,  DATA[4]['pdos'][0, ...],
        c='0.3333',
        label=r'$e_g$'
    )

    AX2.plot(
        x_axis, -DATA[4]['pdos'][1, ...],
        c='0.3333'
    )

    AX2.fill_between(
        x_axis, DATA[4]['pdos'][0, ...], -DATA[4]['pdos'][1, ...],
        facecolor='none',
        edgecolor='0.3333',
        hatch='---'
    )

    AX2.plot(
        x_axis,  DATA[5]['pdos'][0, ...],
        c='0.6667',
        label=r'$t_{2g}$'
    )

    AX2.plot(
        x_axis, -DATA[5]['pdos'][1, ...],
        c='0.6667'
    )

    AX2.fill_between(
        x_axis, DATA[5]['pdos'][0, ...], -DATA[5]['pdos'][1, ...],
        facecolor='none',
        edgecolor='0.6667',
        hatch='|||'
    )

    AX2.plot(
        x_axis,  DATA[-1]['pdos'][0, ...],
        c='black',
        label='Total'
    )

    AX2.plot(
        x_axis, -DATA[-1]['pdos'][1, ...],
        c='black'
    )

    AX2.plot(
        [DATA[0]['e_fermi'], DATA[0]['e_fermi']],
        [np.amax(DATA[-1]['pdos'][0, ...]), -np.amax(DATA[-1]['pdos'][1, ...])],
        linestyle='dashed',
        linewidth=2.0,
        c='red',
        label=r'$\varepsilon_F$'
    )

    AX2.set_title('Bulk-like')

    for i in DATA:
        i.close()

    AX2.set(
        xlabel='Energy (eV)',
        xlim=XLIM,
        ylim=YLIM,
        yticks=[-5, 5]
    )

    AX1.set(ylabel='pDOS (a.u.)')

    AX2.legend(
        loc='lower center',
        bbox_to_anchor=(0.5, -1),
        ncol=3
    )

    plt.figtext(0.01, 0.95, '(a)')

    OUTPUT.savefig()
