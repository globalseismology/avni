"""
This module has several trogonometric routines in degrees.
Typically used to interface with fortran codes.
"""
#####################  IMPORT STANDARD MODULES   #########################

from numpy import cos, sin, tan, arccos, arcsin, arctan, arctan2, degrees, radians
import typing as tp
from numba import jit

##########################################################################

@jit
def cosd(x: tp.Union[float,int]) -> float:
  return cos(radians(x))

@jit
def sind(x: tp.Union[float,int]) -> float:
  return sin(radians(x))

@jit
def tand(x: tp.Union[float,int]) -> float:
  return tan(radians(x))

@jit
def acosd(x: tp.Union[float,int]) -> float:
  return degrees(arccos(x))

@jit
def asind(x: tp.Union[float,int]) -> float:
  return degrees(arcsin(x))

@jit
def atand(x: tp.Union[float,int]) -> float:
  return degrees(arctan(x))

@jit
def atan2d(y: tp.Union[float,int],
           x: tp.Union[float,int]) -> float:
  return degrees(arctan2(y, x))