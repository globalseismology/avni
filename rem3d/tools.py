from __future__ import absolute_import
import scipy.constants
import ConfigParser
import pdb    #for the debugger pdb.set_trace()
import pkgutil
import os
####################### IMPORT REM3D LIBRARIES  #######################################
from . import constants
#######################################################################################

def get_installdir(checkwrite=True):
    """
    Get the installation directory for the rem3d module
    """
    loader=pkgutil.find_loader('rem3d')
    if loader is None:
        print "Warning: installation directory not found for rem3d. Using current directory "
        installdir = '.'
    else:
        installdir = os.path.dirname(loader.get_filename())
        if checkwrite:
            if not os.access(installdir, os.W_OK):
                print "Warning: Cannot I/O to rem3d directory due to permissions. Use current directory. "
                installdir = '.'
    return installdir
    
def get_filedir(checkwrite=True,makedir=True):
    """
    Get the local files directory. Make a new directory if doesn't exist (makedir==True)
    """
    installdir = get_installdir(checkwrite=checkwrite)
    filedir = installdir+'/'+constants.localfilefolder
    if checkwrite and makedir: 
        if not os.path.exists(filedir):
            os.makedirs(filedir)        
    return filedir

