	subroutine conv2geocen(xlatin, eplaout)
c
c convert a geographic coordinate to geocentric coordinate
c

	implicit double precision (a-h,o-z)

	real*4 xlatin, eplaout
	xlat=xlatin


c	data fac /1.d0/
	data fac /0.993305621334896/


	data rad /0.017453292519943/
	theta = (90.0-xlat)*rad
	if(theta.gt.1.0e-20) then 
		theta = 1.570796326794895-atan2(fac*cos(theta),sin(theta))
	else
		theta = 1.570796326794895-atan2(fac*cos(theta),1.0d-20)
	endif

	eplaout=90.d0-theta/rad

	return
	end
 
