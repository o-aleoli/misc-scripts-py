#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np

# List of atoms to search at the CONTCAR file

listNeighbourNa = np.array([
    [55, 60],
    [64, 51],
    [60, 56],
    [51, 65],
    [56, 61],
    [65, 52],
    [61, 57],
    [52, 66],
    [57, 62],
    [66, 53],
    [62, 58],
    [53, 67],
    [58, 63],
    [67, 54],
    [63, 59],
    [54, 68]
])
listMiddleNa = np.array([
    [55, 64],
    [60, 51],
    [56, 65],
    [61, 52],
    [59, 68],
    [63, 54],
    [58, 67],
    [62, 53]
])
middleNa = np.array([57, 66])

posList = np.genfromtxt( 'CONTCAR', skip_header=8, skip_footer=86, delimiter='  ' )
simBox = np.genfromtxt( 'CONTCAR', skip_header=2, skip_footer=175, delimiter='    ' )

distFirstNa = np.zeros( (len( listNeighbourNa ), 3) )

for i in range( len( listNeighbourNa ) ):
"""
CONTCAR is in fractional cartesian coordinates. The atom distance is the matrix product of dist*simBox
More at: https://www.vasp.at/wiki/wiki/index.php/POSCAR
"""
    distFirstNa[i] = np.matmul(np.absolute(posList[listNeighbourNa[i][1] - 1] - posList[listNeighbourNa[i][0] - 1]), simBox)

distMiddleNa = np.zeros( (2, len( listMiddleNa ), 3) )
averagedDistMiddleNa = np.zeros( (int(0.5*listMiddleNa.shape[0]), 3) )

for i in range( len( listMiddleNa ) ):
    distMiddleNa[0][i] = np.matmul(np.absolute(posList[middleNa[0] - 1] - posList[listMiddleNa[i][0] - 1]), simBox)
    distMiddleNa[1][i] = np.matmul(np.absolute(posList[middleNa[1] - 1] - posList[listMiddleNa[i][1] - 1]), simBox)

for i in range( len( averagedDistMiddleNa ) ):
    averagedDistMiddleNa[i] = 0.25*(distMiddleNa[0][i] + distMiddleNa[1][i] + distMiddleNa[0][4+i] + distMiddleNa[1][4+i])

with open('distNa.txt', 'w') as output:
    output.write( '--------------------\nNeighbours Na\n--------------------\nát. át. ΔX (Å)\n' )
    for i in range( len( distFirstNa ) ):
        output.write( '{} {} {:.3E}\n'.format( listNeighbourNa[i][0], listNeighbourNa[i][1], distFirstNa[i][0] ) )
    output.write( '--------------------\nInterlayer distance\n--------------------\n')
    for i in range(  len( averagedDistMiddleNa ) ):
        output.write( '{} 🠆 5 {:.3E}\n'.format( i+1, averagedDistMiddleNa[i][0] ) )
    output.write( '--------------------\nFrom the middle Na\n--------------------\nát. át. ΔX (Å)\n' )
    for i in range( int( 0.5*distMiddleNa.shape[1] ) ):
        output.write( '{} {} {:.3E}\n'.format( listMiddleNa[i][0], middleNa[0], distMiddleNa[0][i][0] ) )
    for i in range( int( 0.5*distMiddleNa.shape[1] ) ):
        output.write( '{} {} {:.3E}\n'.format( listMiddleNa[i][1], middleNa[1], distMiddleNa[1][i][0] ) )

