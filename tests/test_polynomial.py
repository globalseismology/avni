#!/usr/bin/env python

import rem3d
import rem3d.constants as constants
import pdb
import numpy as np

rnorm = constants.R.to('km').magnitude
vercof,dvercof = rem3d.tools.bases.eval_polynomial(rnorm-30.,radius_range=[rnorm-35.,rnorm-20.] ,rnorm=rnorm,types=['CONSTANT','LINEAR'])
assert vercof[0] == 1.,'Check that the value of a CONSTANT basis is 1 within the domain'
vercof,dvercof = rem3d.tools.bases.eval_polynomial(rnorm-150.,radius_range=[rnorm-35.,rnorm-20.] ,rnorm=rnorm,types=['CONSTANT','LINEAR'])
assert np.all(vercof==0.), 'Check that the value of a CONSTANT basis is 0 outside the domain'

vercof,dvercof = rem3d.tools.bases.eval_polynomial([rnorm-150., rnorm-180.],radius_range=[rnorm-200.,rnorm-100.] ,rnorm=rnorm,types=['TOP','BOTTOM'])
top= bottom = 3.45 #m/s
coeff = np.array([[top,bottom]]).T
assert np.all(np.isclose(np.dot(vercof,coeff),top)), 'evaluation of TOP/BOTTOM should give constant value anywhere in the domain'
