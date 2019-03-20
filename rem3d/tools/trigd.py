"""
This module has several trogonometric routines in degrees.
Typically used to interface with fortran codes.
"""

from math import cos, sin, tan, acos, asin, atan, atan2, degrees, radians

def cosd(x):
  return cos(radians(x))

def sind(x):
  return sin(radians(x))

def tand(x):
  return tan(radians(x))

def acosd(x):
  return degrees(acos(x))

def asind(x):
  return degrees(asin(x))

def atand(x):
  return degrees(atan(x))

def atan2d(y, x):
  return degrees(atan2(y, x))