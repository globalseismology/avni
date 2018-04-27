from __future__ import absolute_import
import scipy.constants
import ConfigParser
import pdb    #for the debugger pdb.set_trace()
import pkgutil
import StringIO
####################### IMPORT REM3D LIBRARIES  #######################################
from . import constants
#######################################################################################

def get_installdir(infile='install.cfg'):
    """
    Get the installation directory
    """
    string = pkgutil.get_data('rem3d',constants.localfiles+'/'+infile)
    buf = StringIO.StringIO(string)
    Config = ConfigParser.ConfigParser()
    Config.readfp(buf)
    installdir = Config.get('metadata','installdir')
    return installdir
