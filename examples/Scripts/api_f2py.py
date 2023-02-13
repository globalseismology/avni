#!/usr/bin/env python
''' Example for calling an f2py function through the API '''
from avni.api.client import Client as r3d
from avni.api.f2py import f2pyWrapper
import argparse
import numpy as np


def main():
    # parse
    parser = argparse.ArgumentParser(description='f2py calls via API')
    parser.add_argument('-k', '--key', type=str,default='',help='api key')
    arg = parser.parse_args()

    # initialize API connection
    if arg.key=='':
        conn=r3d()
    else:
        conn=r3d(api_key=arg.key) # the connection object

    # load f2py instance
    f2py=f2pyWrapper(conn)

    # list available functions and arguments
    funlist=f2py.listf2py()
    print("\nf2py function list")
    for fun in funlist.keys():
        print("\nfunction: "+fun)
        print("description: "+funlist[fun]['description'])
        print("input arguments: ")
        print(funlist[fun]['arguments'])
        print("output expected:")
        print(funlist[fun]['output'])

    print("\nBegin function testing!")
    print("\nf2py.callf2py(function='FUNCTIONNAME',args={})")
    f2pyresult=f2py.callf2py(function='FUNCTIONNAME',args={})
    print(f2pyresult)

    # subroutine ddelazgc(eplat,eplong,stlat,stlong,delta,azep,azst)
    print("\nf2py.callf2py(function='ddelazgc',args={'lat1':-25.,'lon1':191.,'lat2':-22.,'lon2':160.})")
    f2pyresult=f2py.callf2py(function='ddelazgc',args={'lat1':-25.,'lon1':191.,'lat2':-22.,'lon2':160.})
    print("Output tuple, "+funlist['ddelazgc']['output']+":")    
    print(f2pyresult)


    splpts=np.linspace(0,3000,50) # will be order F'd on server
    f2py_args={'depth':125.,'Nsplts':len(splpts),'splpts':splpts}
    print("\nf2py.callf2py(function='vbspl',args=f2py_args)")
    f2pyresult=f2py.callf2py(function='vbspl',args=f2py_args)
    print("Output tuple, "+funlist['vbspl']['output']+":")
    print(f2pyresult)

    print("\n")

if __name__== "__main__":
    main()
