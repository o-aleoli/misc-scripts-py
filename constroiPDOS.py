#!/usr/bin/env python3

import xml.etree.ElementTree as et
import numpy as np

def Gaussian ( x, mu, sigma, height ):
    return height*np.exp( -0.5*np.power( (x - mu)/sigma, 2 ) )

def OrbSum ( dos, orb ):
    acc = np.zeros( dos.shape[:-1] )
    arrayIterator = np.nditer( orb, flags=['f_index'] )
    while not arrayIterator.finished:
        if arrayIterator.value == 0:
            arrayIterator.iternext()
        else:
            for spin in range( dos.shape[0] ):
                acc[spin, ...] += dos[spin, ..., arrayIterator.iterindex]
            arrayIterator.iternext()
    return acc

vasprun_root = et.parse( 'vasprun.xml' ).getroot()
# A sequence of atoms can be evaluated by declaring the first and last atom
# A list of atoms in any order can be evaluated by explicitly declaring each one
ion_list = np.array( [1, 86] )
orbitals = np.array( [1, 1, 1, 1, 1,  1,  1,  1,  1] )
#                    [s  py pz px dxy dyz dz2 dxz dx2]
sigma = 0.1
ispin = int( vasprun_root.find( './/*[@name="ISPIN"]' ).text )
nedos = int( vasprun_root.find( './/*[@name="NEDOS"]' ).text )
efermi = float( vasprun_root.find( './/*[@name="efermi"]' ).text )
numpoints = nedos
eigen_adjust = 0.0
outputPath = "./total.dat"

if len( ion_list ) == 2:
    ion_list = np.arange( ion_list[0], ion_list[1] + 1 )

pdos_raw = np.empty( (len( ion_list ), ispin, nedos, len( orbitals ) + 1) )
ionIterator = np.nditer( ion_list, flags=['f_index'] )

while not ionIterator.finished:
    print( "extracting data from ion {} ({}/{})".format( np.array2string( ionIterator.value ), ionIterator.iterindex + 1, ionIterator.itersize ) )
    for spin in range( ispin ):
        for numdos in range( nedos ):
            pdos_raw[ionIterator.iterindex, spin, numdos] = np.fromstring( vasprun_root.find( './/partial/array/set' )[ionIterator[0] - 1][spin][numdos].text, sep=' ' )
    ionIterator.iternext()

vasprun_root.clear()

orbIterator = np.nditer( orbitals, flags=['f_index'] )
pdos = np.zeros( (ispin, nedos, len( orbitals )) )
eigen_energy = np.linspace( pdos_raw[0, 0, 0, 0], pdos_raw[0, 0, -1, 0], numpoints )

while not orbIterator.finished:
    if orbIterator.value == 0:
        orbIterator.iternext()
    else:
        for ion in range( len( ion_list ) ):
            print( "summing orbital #{} of ion {}/{}".format( orbIterator.iterindex + 1, ion + 1, len( ion_list ) ) )
            for spin in range( ispin ):
                pdos[spin, ..., orbIterator.iterindex] += pdos_raw[ion, spin, :, orbIterator.iterindex + 1]
        orbIterator.iternext()

pdos_total = OrbSum( pdos, np.array( [1, 1, 1, 1, 1, 1, 1, 1, 1] ) )
pdos_total_gaussian = np.zeros( (ispin, numpoints) )
for spin in range( ispin ):
    for numdos in range( pdos_total.shape[1] ):
        pdos_total_gaussian[spin, ...] += Gaussian( eigen_energy, pdos_raw[0, spin, numdos, 0] + eigen_adjust, sigma, pdos_total[spin, numdos] )

print( "writing PDOS on file "+outputPath )
with open( outputPath, 'w' ) as output:
    for numdos in range( pdos_total_gaussian.shape[1] ):
        output.write( "{:f} {:f}\n".format( pdos_total_gaussian[0, numdos], eigen_energy[numdos] ) )
    output.write( "\n" )

    if pdos_total_gaussian.shape[0] == 1:
        for numdos in range( pdos_total_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( -pdos_total_gaussian[0, numdos], eigen_energy[numdos] ) )
    else:
        for numdos in range( pdos_total_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( -pdos_total_gaussian[1, numdos], eigen_energy[numdos] ) )
