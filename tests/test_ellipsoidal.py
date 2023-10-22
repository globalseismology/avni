#!/usr/bin/env python

from avni.mapping import ellipsoidal
import math

delta,azep,azst = ellipsoidal.get_distaz(0.,0.,0.,90.)
assert math.isclose(delta, 90., abs_tol=0.001),'Check ellipsodal get_distaz; incorrect delta'

delta,azep,azst= ellipsoidal.get_distaz(0.,0.,45.,0.)
geocentric_lat = ellipsoidal.geographic_to_geocentric(45.)
assert math.isclose(delta, geocentric_lat, abs_tol=0.001),'Check ellipsodal get_distaz and geographic_to_geocentric; incorrect delta'
assert math.isclose(azst, 180., abs_tol=0.001),'Check ellipsodal get_distaz; incorrect backazimuth azst'