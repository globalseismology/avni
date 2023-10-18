#!/usr/bin/env python

from avni.data.SW import SWasciitohdf5
import glob
import os.path
import pandas as pd
import argparse #parsing arguments

#########################  MAIN   ######################################################

def main():
    parser = argparse.ArgumentParser(description='convert list of files to hdf5')
    parser.add_argument("files", type=str,
        help='Input file with list of files to add')
    parser.add_argument('-o', '--output', type=str, default='SW.avni.data.h5',
        help='Output hdf5 file to write or append to'
        )
    parser.add_argument('-d', '--datatype', type=str, default='summary',
        help='Type of data so that the appropriate folder is updated in hdf5'
        )

    arg = parser.parse_args()
    files = pd.read_csv(arg.files,header=None)
    SWasciitohdf5(files, hdffile = arg.output, datatype = arg.datatype)
    return

if __name__ == "__main__":
    main()
#######################################################################################


