"""
This module has several trogonometric routines in degrees.
Typically used to interface with fortran codes.
"""

from math import cos, sin, tan, acos, asin, atan, atan2, degrees, radians
from numba import jit

@jit(nopython=True)
def cosd(x):
  return cos(radians(x))

@jit(nopython=True)
def sind(x):
  return sin(radians(x))

@jit(nopython=True)
def tand(x):
  return tan(radians(x))

@jit(nopython=True)
def acosd(x):
  return degrees(acos(x))

@jit(nopython=True)
def asind(x):
  return degrees(asin(x))

@jit(nopython=True)
def atand(x):
  return degrees(atan(x))

@jit(nopython=True)
def atan2d(y, x):
  return degrees(atan2(y, x))