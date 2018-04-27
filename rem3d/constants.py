from __future__ import absolute_import
import scipy.constants
import ConfigParser
import pdb    #for the debugger pdb.set_trace()
import pkgutil
import StringIO

"""
download files
"""
downloadpage = 'https://maurya.umd.edu/files'
localfiles = 'files'

"""
Installation directory
"""
string = pkgutil.get_data('rem3d',localfiles+'/install.cfg')
buf = StringIO.StringIO(string)
Config = ConfigParser.ConfigParser()
Config.readfp(buf)
installdir = Config.get('metadata','installdir')