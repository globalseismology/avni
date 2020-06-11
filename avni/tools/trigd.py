"""
This module has several trogonometric routines in degrees.
Typically used to interface with fortran codes.
"""

from math import cos, sin, tan, acos, asin, atan, atan2, degrees, radians
from numba import jit

@jit
def cosd(x):
  return cos(radians(x))

@jit
def sind(x):
  return sin(radians(x))

@jit
def tand(x):
  return tan(radians(x))

@jit
def acosd(x):
  return degrees(acos(x))

@jit
def asind(x):
  return degrees(asin(x))

@jit
def atand(x):
  return degrees(atan(x))

@jit
def atan2d(y, x):
  return degrees(atan2(y, x))