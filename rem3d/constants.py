from __future__ import absolute_import
import scipy.constants
import ConfigParser
import pdb    #for the debugger pdb.set_trace()
import pkgutil
import StringIO

####################### IMPORT REM3D LIBRARIES  #######################################

from . import tools
#######################################################################################

"""
download files
"""
downloadpage = 'https://maurya.umd.edu/files'
localfiles = 'files'

"""
Installation directory
"""
installdir = tools.get_installdir('install.cfg')
