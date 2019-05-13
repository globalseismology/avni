      subroutine getflpos(lufl,iposition,ierror)
      include "giolib.h"
c
      ierror=1
      if(lufl.gt.0.and.lufl.le.20) then
	if(iopen(lufl).gt.0) then
	  iposition=ipos(lufl)
	  ierror=0
	endif
      else
	write(6,"('lufl outside range:',i10)") lufl
	stop
      endif
      return
      end
