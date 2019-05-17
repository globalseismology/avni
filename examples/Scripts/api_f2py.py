#!/usr/bin/env python
''' Example for calling an f2py function through the API '''
from rem3d.api.client import Client as r3d
from rem3d.api.f2py import f2pyWrapper
import argparse


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
    print(funlist['function_list'])
    print("funlist['argument'][function_name] will return arguments")


    print("\n begin function testing")
    f2pyresult=f2py.callf2py(f2pyfun='FUNCTIONNAME',f2py_args=[])
    print("\nf2py.callf2py(f2pyfun='FUNCTIONNAME',f2py_args=[])")
    print(f2pyresult)

    # subroutine ddelazgc(eplat,eplong,stlat,stlong,delta,azep,azst)
    f2pyresult=f2py.callf2py(f2pyfun='ddelazgc',f2py_args=[-25.,191.,-22.,160.])
    print("\nf2py.callf2py(f2pyfun='ddelazgc',f2py_args=[-25.,191.,-22.,160.])")
    print(f2pyresult)

    print("\n")

if __name__== "__main__":
    main()
