#!/usr/bin/env python3
'''Alexandre Olivieri (olivieri.alexandre0@gmail.com)

Post-processing script for VASP output

It extracts the electronic state's projection on each atom and orbital.
It sums over the declared number of atoms and the declared orbitals.
At the end, the script changes the Dirac-delta-like states into Gaussian-like
and builds the Y (or X) array.
The output is a pickled .npz file.'''

import xml.etree.ElementTree as et
import numpy as np

def Gaussian ( x, mu, sigma, height ):
    return height*np.exp( -0.5*np.power( (x - mu)/sigma, 2 ) )

def OrbSum ( dos, orb ):
    acc = np.zeros( dos.shape[:-1] )
    array_iterator = np.nditer( orb, flags=['f_index'] )
    while not array_iterator.finished:
        if array_iterator.value == 0:
            array_iterator.iternext()
        else:
            for spin in range( dos.shape[0] ):
                acc[spin, ...] += dos[spin, ..., array_iterator.iterindex]
            array_iterator.iternext()
    return acc

vasprun_root = et.parse( 'vasprun.xml' ).getroot()
# A sequence of atoms can be evaluated by declaring the first and last atom
# A list of atoms in any order can be evaluated by explicitly declaring each one
ion_list = np.array( [69, 86] )
orbitals = np.array( [0, 0, 0, 0, 1,  1,  1,  1,  1] )
#                    [s  py pz px dxy dyz dz2 dxz dx2]
eigen_adjust = 2.184686
pdos_on_y = True
output_path = "./Nb_total_d"
sigma = 0.1
ispin = int( vasprun_root.find( './/*[@name="ISPIN"]' ).text )
nedos = int( vasprun_root.find( './/*[@name="NEDOS"]' ).text )
e_fermi = float( vasprun_root.find( './/*[@name="efermi"]' ).text )
numpoints = nedos

if len( ion_list ) == 2:
    ion_list = np.arange( ion_list[0], ion_list[1] + 1 )

pdos_raw = np.empty( (len( ion_list ), ispin, nedos, len( orbitals ) + 1) )
ion_iterator = np.nditer( ion_list, flags=['f_index'] )

while not ion_iterator.finished:
    print( "extracting data from ion {} ({}/{})".format( np.array2string( ion_iterator.value ), ion_iterator.iterindex + 1, ion_iterator.itersize ) )
    for spin in range( ispin ):
        for numdos in range( nedos ):
            pdos_raw[ion_iterator.iterindex, spin, numdos] = np.fromstring( vasprun_root.find( './/partial/array/set' )[ion_iterator[0] - 1][spin][numdos].text, sep=' ' )
    ion_iterator.iternext()

vasprun_root.clear()

orb_iterator = np.nditer( orbitals, flags=['f_index'] )
acc = np.zeros( (ispin, nedos, len( orbitals )) )
eigen_energy = np.linspace( pdos_raw[0, 0, 0, 0], pdos_raw[0, 0, -1, 0], numpoints )

while not orb_iterator.finished:
    if orb_iterator.value == 0:
        orb_iterator.iternext()
    else:
        for ion in range( len( ion_list ) ):
            print( "summing orbital #{} of ion {}/{}".format( orb_iterator.iterindex + 1, ion + 1, len( ion_list ) ) )
            for spin in range( ispin ):
                acc[spin, ..., orb_iterator.iterindex] += pdos_raw[ion, spin, :, orb_iterator.iterindex + 1]
        orb_iterator.iternext()

pdos_sum = OrbSum( acc, orbitals )
pdos_gaussian = np.zeros( (ispin, numpoints) )
for spin in range( ispin ):
    for numdos in range( pdos_sum.shape[1] ):
        pdos_gaussian[spin, ...] += Gaussian( eigen_energy, pdos_raw[0, spin, numdos, 0], sigma, pdos_sum[spin, numdos] )

np.savez_compressed(
    output_path+'.npz',
    pdos=pdos_gaussian,
    eigen_energy=eigen_energy,
    eigen_adjust=eigen_adjust,
    e_fermi=e_fermi
    )