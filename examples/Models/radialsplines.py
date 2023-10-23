#!/usr/bin/env python
"""This module contains an example of outputting radial splines using avni codes"""

import argparse #parsing arguments
import numpy as np
################################ IMPORT AVNI MODULES   #####################################

import avni

#########################################################
def main():
    parser = argparse.ArgumentParser(description='outputs the radial splines')
    parser.add_argument('-u', '--upper_bound', type=float, default=2891.,
        help='Upper bound for depth')
    parser.add_argument('-l', '--lower_bound', type=float, default=25.,
        help='Lower bound for depth')
    parser.add_argument('-i', '--interval', type=float, default=50.,
        help='Depth interval in km')
    parser.add_argument('-k','--knots', nargs='+', help='Knot locations in km. Can give multiple values e.g. 24.4 75. 150. 225. 300. 410. 530. 650. 650. 820. 1320. 1820. 2320. 2550. 2791. 2891.', required=True)

    arg = parser.parse_args()

    # Get spline knots
    knots = [float(knot) for knot in arg.knots]
    #Get depths
    start = float(arg.lower_bound);end=float(arg.upper_bound);inter=float(arg.interval)
    depths = np.arange(start,end,inter)
    vercof, _ = avni.tools.eval_vbspl(depths,knots)
    writearr = np.concatenate((depths.reshape((-1, 1)),vercof),1)
    np.savetxt('vercof.txt',writearr,fmt='%.5f')
    return

if __name__== "__main__":
    main()

