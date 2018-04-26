from __future__ import absolute_import
import scipy.constants
import ConfigParser
Config = ConfigParser.ConfigParser()
Config.read('files/install.cfg')

"""
download files
"""
downloadpage = 'https://maurya.umd.edu/files'
localfolder = 'files'
installdir = Config.get('metadata','installdir')
