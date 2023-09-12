#!/usr/bin/env python

#####################  IMPORT STANDARD MODULES   #########################

# python 3 compatibility
from __future__ import absolute_import, division, print_function, unicode_literals
import sys
if (sys.version_info[:2] < (3, 0)):
    from builtins import tuple

import numpy as np
import gc
import warnings
from scipy import sparse
import h5py

#######################################################################################

def close_h5py():
    """Close all h5py files

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    for obj in gc.get_objects():   # Browse through ALL objects
        if isinstance(obj, h5py.File):   # Just HDF5 files
            try:
                obj.close()
            except:
                warnings.warn('Warning: HDF5 files already closed')
                pass # Was already closed

def store_sparse_hdf(h5f,varname: str,mat,compression: str = "gzip"):
    """Store a `csr` matrix in HDF5

    Parameters
    ----------
    h5f
        HDF5 file handle
    varname : str
        node prefix in HDF5 hierarchy
    mat : scipy.sparse.csr.csr_matrix
        sparse matrix to be stored
    compression : str, optional
        Compression type in HDF5, by default "gzip"

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """

    # Check the vector type
    msg = "This code only works for csr matrices"
    if not mat.__class__ == sparse.csr.csr_matrix: raise AssertionError(msg)
    try:  # Try loading the sparse array if it exists
        mat_original = load_sparse_hdf(h5f,varname)
        mat_write = sparse.vstack([mat_original,mat])
        del(h5f[varname])
    except KeyError:
        mat_write = mat

    # Write to a file
    for par in ('data', 'indices', 'indptr', 'shape'):
        arr = np.array(getattr(mat_write, par))
        h5f.create_dataset(varname+'/'+par, data=arr, compression=compression)


def load_sparse_hdf(h5f,varname: str):
    """Load a `csr` matrix from HDF5 file

    Parameters
    ----------
    h5f
        HDF5 file handle
    varname : str
        node prefix in HDF5 hierarchy

    Returns
    -------
    scipy.sparse.csr.csr_matrix
        A sparse `csr` matrix

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    # Check the vector type
    pars = []
    for par in ('data', 'indices', 'indptr', 'shape'):
        pars.append(h5f[varname][par].value)
    m = sparse.csr_matrix(tuple(pars[:3]), shape=pars[3])
    return m

def store_numpy_hdf(h5f,varname: str,array: np.ndarray,compression: str = "gzip", compression_opts: int = 9):
    """Store a named numpy array in HDF5 file

    Parameters
    ----------
    h5f
        HDF5 file handle
    varname : str
        node prefix in HDF5 hierarchy
    array : np.ndarray
        Named numpy array
    compression : str, optional
        Compression type in HDF5, by default "gzip"
    compression_opts : int, optional
        Compression level opts, by default 9

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    # Check if it is a named numpy array
    if not isinstance(array, np.ndarray) : raise ValueError('Only numpy arrays can be stored with store_numpy_hdf')
    if array.dtype.names is None:
        raise ValueError('Only named numpy arrays are allowed')
    else:
        fields = np.array(array.dtype.names,dtype='a15')

    try:  # Try loading the sparse array if it exists
        arr_original = load_numpy_hdf(h5f,varname)
        arr_write = np.hstack([arr_original,array])
        del(h5f[varname])
        print('Warning: appending to existing field: '+varname)
    except:
        arr_write = array

    # Write the file
    h5f.create_dataset(varname+'/fields',data=fields,compression=compression, compression_opts=compression_opts)
    for field in fields:
        # if string, change to utf for python2/3 compatibility
        if arr_write[field].dtype.kind == 'S' or arr_write[field].dtype.kind == 'U':
            outarr=np.array(arr_write[field].tolist(),dtype='a'+str(arr_write[field].dtype.itemsize))
            h5f.create_dataset(varname+'/columns/'+field, data=outarr,compression=compression, compression_opts=compression_opts)
        else:
            h5f.create_dataset(varname+'/columns/'+field, data=arr_write[field], compression=compression, compression_opts=compression_opts)

def load_numpy_hdf(h5f,varname: str) -> np.ndarray:
    """Read a named numpy array from HDF5 file

    Parameters
    ----------
    h5f
        HDF5 file handle
    varname : str
        node prefix in HDF5 hierarchy

    Returns
    -------
    np.ndarray
        Named numpy array

    :Authors:
        Raj Moulik (moulik@caa.columbia.edu)
    :Last Modified:
        2023.02.16 5.00
    """
    if (sys.version_info[:2] > (3, 0)):
        names = h5f[varname]['fields'].value
    else:
        names = [name.decode('utf-8') for name in h5f[varname]['fields'].value]
    formats = [h5f[varname]['columns'][field].dtype.kind+ str(h5f[varname]['columns'][field].dtype.itemsize) for field in names]
    dt = {'names':names, 'formats':formats}
    output = np.zeros(h5f[varname]['columns'][names[0]].value.shape, dtype=dt)
    for field in names:
        output[field]=h5f[varname]['columns'][field].value
    return output
