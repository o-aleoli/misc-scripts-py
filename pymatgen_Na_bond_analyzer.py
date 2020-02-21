#!/usr/bin/env python3
from pymatgen import Structure
from pymatgen.analysis.bond_valence import BVAnalyzer
from pymatgen.analysis.local_env import CrystalNN
import matplotlib.pyplot as plt
import numpy as np

sites_Na = [[55, 59, 64, 68], [51, 54, 60, 63], [56, 58, 65, 67], [52, 53, 62, 64], [57, 57, 67, 67]]
directories = ["c5", "c4", "c3", "c2", "c1", "relax", "t1", "t2", "t3", "t4", "t5"]
o_fig = o_file = 'NaObonds'
all_dist_NaO = []# (strain,layers,atom,bond)

for element in directories:
"""
This populates the distances array with all the Na-O bonds detected by the CrystalNN routine.
The innermost index varies in size as each atom in sites_Na can have more or less neigbours.
The second innermost index size is always 4 due to array symmetry reasons. It registers the number of Na atoms at a layer
The third innermost index is 5. It is the number of layers in the slab.
The first index registers the strain steps taken.
"""
    input_struct = Structure.from_file( element+"/CONTCAR" )
    ox_struct = BVAnalyzer().get_oxi_state_decorated_structure( input_struct )
    dist_NaO = []
    for layer, atom_list in enumerate( sites_Na ):
        layer_atoms = []
        for atom in atom_list:
            atom_neigbors = CrystalNN( cation_anion=True ).get_nn( ox_struct, atom-1 )
            distances = []
            for neighbor in range( len( atom_neigbors ) ):
                distances.append( atom_neigbors[neighbor].distance( input_struct[atom-1] ) )
            layer_atoms.append( distances )
        dist_NaO.append( layer_atoms )
    all_dist_NaO.append( dist_NaO )

avrg_bond = np.zeros( [len(all_dist_NaO), len(all_dist_NaO[0]), len(all_dist_NaO[0][0])] )
avrg_layer = np.zeros( [len(all_dist_NaO), len(all_dist_NaO[0])] )
stdev_bond = np.zeros( [len(all_dist_NaO), len(all_dist_NaO[0]), len(all_dist_NaO[0][0])] )
stdev_layer = np.zeros( [len(all_dist_NaO), len(all_dist_NaO[0])] )

iterator = np.nditer( avrg_bond, flags=['multi_index'] )

for bond in iterator:
"""
Average bond length and standard deviation of each bond around each Na atom.
Here st dev is taken to measure how the polyhedra is deformed. A st dev close to 0 points to a regular polyhedra.
"""
    avrg_bond[iterator.multi_index[0]][iterator.multi_index[1]][iterator.multi_index[2]] = np.average( all_dist_NaO[iterator.multi_index[0]][iterator.multi_index[1]][iterator.multi_index[2]] )
    stdev_bond[iterator.multi_index[0]][iterator.multi_index[1]][iterator.multi_index[2]] = np.std( all_dist_NaO[iterator.multi_index[0]][iterator.multi_index[1]][iterator.multi_index[2]] )

iterator = np.nditer( avrg_bond, flags=['multi_index'] )

for atoms in iterator:
"""
Same logic of last loop but taking the average and st dev of the 4 atoms in each layer.
"""
    avrg_layer[iterator.multi_index[0]][iterator.multi_index[1]] = np.average( avrg_bond[iterator.multi_index[0]][iterator.multi_index[1]] )
    stdev_layer[iterator.multi_index[0]][iterator.multi_index[1]] = np.std( stdev_bond[iterator.multi_index[0]][iterator.multi_index[1]] )

iterator = np.nditer( avrg_layer, flags=['multi_index'] )
xaxis = np.array( range( -5,6 ) )

for layers in iterator:
"""
Creates a plt object with the NaO bond distance in y, strain in x and st dev as symmetric error lines
"""
    plt.errorbar( xaxis,
                  avrg_layer[..., iterator.multi_index[1]],
                  yerr=stdev_layer[..., iterator.multi_index[1]],
                  label='layer '+str( iterator.multi_index[1]+1 ),
                  marker='o',
                  linestyle='dashed',
                  elinewidth=1,
                  capsize=2,
                  capthick=1 )

plt.xlabel( 'Strain (%)' )
plt.xticks( xaxis )
plt.ylabel( 'Na-O bond length (Å)' )
plt.legend()
plt.savefig( o_fig+'.pdf', dpi=300, orientation='portrait', papertype='a4', format='pdf', transparent=True )

np.savetxt( o_file+'.csv', np.hstack( (avrg_layer, stdev_layer) ), fmt='%-3.4f' )# (strain, layer(0 to 4) and stdev(5 to 9))

