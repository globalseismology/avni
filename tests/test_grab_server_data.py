#!/usr/bin/env python
"""This module contains an example of grabbing a file from REM3D server"""

import argparse #parsing arguments

########################### IMPORT REM3D MODULES   #####################################
import rem3d
import pytest
#########################################################
@pytest.fixture
def test_download(file,folder):
    """download a file"""
    rem3d.data.update_file(file,folder)

def main():
    parser = argparse.ArgumentParser(description='grab a file from the REM3D server')
    parser.add_argument('-i', '--file', type=str, default='S362ANI+M.vs.5.epix',
        help='file from REM3D server')
    parser.add_argument('-f', '--folder', type=str,default='.',
        help='output folder')
    arg = parser.parse_args()

    # update the file from the server
    rem3d.data.update_file(arg.file,arg.folder)
    return

if __name__== "__main__":
    main()
