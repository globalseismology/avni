#!/usr/bin/env python
"""This module contains an example of grabbing a file from REM3D server"""

import argparse #parsing arguments

########################### IMPORT REM3D MODULES   #####################################
import rem3d
#########################################################
def func():
    file = 'S362ANI+M.vs.5.epix'
    folder = '.'
    rem3d.data.update_file(file,folder)

def test_server():
    func()
