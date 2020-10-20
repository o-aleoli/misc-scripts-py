#!/usr/bin/env python3
# Alexandre Olivieri (olivieri.alexandre0@gmail.com)
#
# Post-processing script for VASP output
#
# It extracts the pickled .npz in a plain text .dat file with formatted columns

import sys, getopt
import numpy as np

def main( argv ):
    try:
        opts, args = getopt.getopt( argv, "hi:o:p:", ["in=", "out=", "pdos_on_y="] )
    except getopt.GetoptError:
        print( "pickled_extractor.py -i <input> -o <output> -p <pDOS on Y axis?>" )
        sys.exit(2)
    
    output_path = None
    input_file = None
    pdos_on_y = True
    for opt, arg in opts:
        if opt == "-h":
            print( "pickled_extractor.py -i <input> -o <output>" )
            sys.exit()
        elif opt in ( "-i", "--in" ):
            input_file = arg
        elif opt in ( "-o", "--out" ):
            output_path = arg
        elif opt in ( "-p", "--pdos_on_y" ):
            if arg == "True": 
                pdos_on_y = True
            elif arg == "False":
                pdos_on_y = False
            else:
                print( "-p should be either <True> or <False>" )
                sys.exit()

    pdos_gaussian = np.load( input_file + ".npz" )['pdos']
    eigen_energy = np.load( input_file + ".npz" )['eigen_energy']
    eigen_adjust = np.load( input_file + ".npz" )['eigen_adjust']
#   eigen_adjust = 2.184686

    with open( output_path+'.dat', 'w' ) as output:
        if pdos_on_y == True:
            for numdos in range( pdos_gaussian.shape[1] ):
                output.write( "{:f} {:f}\n".format(
                    eigen_energy[numdos] + eigen_adjust,
                    pdos_gaussian[0, numdos]
                    ))
            output.write( "\n" )
    
            if pdos_gaussian.shape[0] == 1:
                for numdos in range( pdos_gaussian.shape[1] ):
                    output.write( "{:f} {:f}\n".format(
                        eigen_energy[numdos] + eigen_adjust,
                        -pdos_gaussian[0, numdos]
                        ))
            else:
                for numdos in range( pdos_gaussian.shape[1] ):
                    output.write( "{:f} {:f}\n".format(
                        eigen_energy[numdos] + eigen_adjust,
                        -pdos_gaussian[1, numdos]
                        ))
            output.write( "\n" )
        else:
            for numdos in range( pdos_gaussian.shape[1] ):
                output.write( "{:f} {:f}\n".format(
                    pdos_gaussian[0, numdos],
                    eigen_energy[numdos] + eigen_adjust
                    ))
            output.write( "\n" )
    
            if pdos_gaussian.shape[0] == 1:
                for numdos in range( pdos_gaussian.shape[1] ):
                    output.write( "{:f} {:f}\n".format(
                        -pdos_gaussian[0, numdos],
                        eigen_energy[numdos] + eigen_adjust
                        ))
            else:
                for numdos in range( pdos_gaussian.shape[1] ):
                    output.write( "{:f} {:f}\n".format(
                        -pdos_gaussian[1, numdos],
                        eigen_energy[numdos] + eigen_adjust
                        ))
            output.write( "\n" )

if __name__ == "__main__":
    main(sys.argv[1:])