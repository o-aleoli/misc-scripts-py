#!/usr/bin/env python3
import xml.etree.ElementTree as et
import numpy as np

def Gaussian ( x, mu, sigma, height ):
"""
Each evaluated point at mu is associated to a gaussian with FWHM sigma.
Works better when the vector x have more elements than mu.
"""
    return height*np.exp( -0.5*np.power( (x - mu)/sigma, 2 ) )

def OrbSum ( dos, orb ):
"""
Simply sums all the elements at the 2 index when the element in orb is one.
"""
    acc = np.zeros( dos.shape[:-1] )
    arrayIterator = np.nditer( orb, flags=['f_index'] )
    while not arrayIterator.finished:
        if arrayIterator.value == 0:
            arrayIterator.iternext()
        else:
            for spin in range( dos.shape[0] ):
                acc += dos[spin, :, arrayIterator.iterindex]
            arrayIterator.iternext()
    return acc

vasprun_root = et.parse( 'vasprun.xml' ).getroot()
# A sequence of atoms can be evaluated by declaring the first and last atom
# A list of atoms in any order can be evaluated by explicitly declaring each one
ion_list = np.array( [1, 86] )
orbitals = np.array( [1, 1, 1, 1, 1, 1, 1, 1, 1] )# [s py pz px dxy dyz dz2 dxz dx2]
sigma = 0.1
# Search at xml if the simulation solved for spin pseudo wavefunctions
ispin = int( vasprun_root.find( './/*[@name="ISPIN"]' ).text )
# Search at xml the number of states calculed
nedos = int( vasprun_root.find( './/*[@name="NEDOS"]' ).text )
efermi = float( vasprun_root.find( './/*[@name="efermi"]' ).text )
# Name of the plain text file that will save occupation data in pairs of columns
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
eigen_energy = np.linspace( pdos_raw[0, 0, 0, 0], pdos_raw[0, 0, -1, 0], 2*nedos )

while not orbIterator.finished:
"""
A quick and dirty way to sum the contribution from each angular momenta into a significant figure of merit.
Generally you will want to use it to, e.g., sum px, py and px into P but it can also be used to sum in
hybridizations such as dxy + dxz + dyz = e_g.
"""
    if orbIterator.value == 0:
        orbIterator.iternext()
    else:
        for ion in range( len( ion_list ) ):
            print( "summing orbital #{} of ion {}/{}".format( orbIterator.iterindex + 1, ion + 1, len( ion_list ) ) )
            for spin in range( ispin ):
                pdos[spin, :, orbIterator.iterindex] += pdos_raw[ion, spin, :, orbIterator.iterindex + 1]
        orbIterator.iternext()

pdos_s_gaussian = np.zeros( (ispin, 2*nedos) )
for spin in range( ispin ):
"""
Evokes the Gaussian function for each energy level and sum it according to the spin.
The spherical symmetry of S does not decompose into angular momenta values, so there is no sum over it.
"""
    for numdos in range( pdos.shape[1] ):
        pdos_s_gaussian[spin, :] += Gaussian( eigen_energy, pdos_raw[0, spin, numdos, 0], sigma, pdos[spin, numdos, 0] )

pdos_p = OrbSum( pdos, np.array( [0, 1, 1, 1, 0, 0, 0, 0, 0] ) )
pdos_p_gaussian = np.zeros( (ispin, 2*nedos) )
for spin in range( ispin ):
"""
Same logic for the last loop but for total P.
"""
    for numdos in range( pdos_p.shape[1] ):
        pdos_p_gaussian[spin, :] += Gaussian( eigen_energy, pdos_raw[0, spin, numdos, 0], sigma, pdos_p[spin, numdos] )

pdos_d = OrbSum( pdos, np.array( [0, 0, 0, 0, 1, 1, 1, 1, 1] ) )
pdos_d_gaussian = np.zeros( (ispin, 2*nedos) )
for spin in range( ispin ):
    for numdos in range( pdos_d.shape[1] ):
        pdos_d_gaussian[spin, :] += Gaussian( eigen_energy, pdos_raw[0, spin, numdos, 0], sigma, pdos_d[spin, numdos] )

pdos_total = OrbSum( pdos, np.array( [1, 1, 1, 1, 1, 1, 1, 1, 1] ) )
pdos_total_gaussian = np.zeros( (ispin, 2*nedos) )
for spin in range( ispin ):
    for numdos in range( pdos_d.shape[1] ):
        pdos_total_gaussian[spin, :] += Gaussian( eigen_energy, pdos_raw[0, spin, numdos, 0], sigma, pdos_total[spin, numdos] )

print( "writing PDOS on file "+outputPath )
with open( outputPath, 'w' ) as output:
    for numdos in range( pdos_s_gaussian.shape[1] ):
        output.write( "{:f} {:f}\n".format( eigen_energy[numdos], pdos_s_gaussian[0, numdos] ) )
    
    output.write( "\n" )
    if pdos_s_gaussian.shape[0] == 1:
"""
Writing of the second spin component at all cases for aesthetic reasons.
If there is no spin polarization, then the negative values are identical to teh positive ones.
If else, the true values are written.
"""
        for numdos in range( pdos_s_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( eigen_energy[numdos], -pdos_s_gaussian[0, numdos] ) )
    else:
        for numdos in range( pdos_s_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( eigen_energy[numdos], -pdos_s_gaussian[1, numdos] ) )
    
    output.write( "\n" )
    for numdos in range( pdos_p_gaussian.shape[1] ):
        output.write( "{:f} {:f}\n".format( eigen_energy[numdos], pdos_p_gaussian[0, numdos] ) )

    output.write( "\n" )    
    if pdos_p_gaussian.shape[0] == 1:
        for numdos in range( pdos_p_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( eigen_energy[numdos], -pdos_p_gaussian[0, numdos] ) )
    else:
        for numdos in range( pdos_p_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( eigen_energy[numdos], -pdos_p_gaussian[1, numdos] ) )
    
    output.write( "\n" )
    for numdos in range( pdos_d_gaussian.shape[1] ):
        output.write( "{:f} {:f}\n".format( eigen_energy[numdos], pdos_d_gaussian[0, numdos] ) )
    
    output.write( "\n" )
    if pdos_d_gaussian.shape[0] == 1:
        for numdos in range( pdos_d_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( eigen_energy[numdos], -pdos_d_gaussian[0, numdos] ) )
    else:
        for numdos in range( pdos_d_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( eigen_energy[numdos], -pdos_d_gaussian[1, numdos] ) )
    
    output.write( "\n" )
    for numdos in range( pdos_total_gaussian.shape[1] ):
        output.write( "{:f} {:f}\n".format( eigen_energy[numdos], pdos_total_gaussian[0, numdos] ) )
    
    output.write( "\n" )
    if pdos_total_gaussian.shape[0] == 1:
        for numdos in range( pdos_total_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( eigen_energy[numdos], -pdos_total_gaussian[0, numdos] ) )
    else:
        for numdos in range( pdos_total_gaussian.shape[1] ):
            output.write( "{:f} {:f}\n".format( eigen_energy[numdos], -pdos_total_gaussian[1, numdos] ) )
