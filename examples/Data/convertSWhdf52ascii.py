#!/usr/bin/env python

from avni.data.SW import SWhdf5toascii
import glob
import os.path
import pandas as pd
import argparse #parsing arguments

#########################  MAIN   ######################################################

def main():
    parser = argparse.ArgumentParser(description='convert list of files to hdf5')
    parser.add_argument("file", type=str,
        help='Input file with list of files to add')
    parser.add_argument('-i', '--overtone', type=int, default=0, help='Overtone branch')
    parser.add_argument('-p', '--period', type=float, default=25.0, help='Period in s')
    parser.add_argument('-w', '--wave', type=str, default='L1', help='Wave type')
    parser.add_argument('-g', '--group', type=str, default='REM3D', help='Research group')
    parser.add_argument('-o', '--output', type=str, default='SW.avni.data.h5',
        help='Output hdf5 file to write or append to'
        )
    parser.add_argument('-d', '--datatype', type=str, default='summary',
        help='Type of data so that the appropriate folder is updated in hdf5'
        )

    arg = parser.parse_args()
    query = str(arg.overtone)+'/{:.1f}'.format(arg.period)+'/'+arg.wave+'/'+arg.group
    SWhdf5toascii(query=query, hdffile = arg.file, datatype = arg.datatype)
    return

if __name__ == "__main__":
    main()
#######################################################################################


