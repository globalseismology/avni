      subroutine closfl(lufl,istat)
      include "giolib.h"
c
      if(lufl.gt.0.and.lufl.le.20) then
	if(iopen(lufl).gt.0) then
	  lup=lufl+1000
	  close(lup)
	  iopen(lufl)=-1
	  istat=2
	else
	  istat=3
	endif
      else
	write(6,"('lufl outside range:',i10)") lufl
	stop
      endif
      return
      end
