      subroutine loadnbn2memory(igcmtchoice,allorder,qcmt,nread,ierror)
c
c---- this subroutine returns cmt date in the 2014-01-08 00:00:13 format for use in python
c     igmtchoice is 1 if both monthly+qcmt are used (allorder.nbn & qcmt.nbn)
c                   2 if only monthly are used (allorder.nbn)
c                   3 if qcmt is used
c
      logical exists
      character*4 catalog
      character*100 filein,allorder,qcmt
	  real*4 areal
c
c
Cf2py intent(inout) igcmtchoice,allorder,qcmt
Cf2py intent(out) nread,ierror

c
c---- open in the .nbn file
c
      do i=1,2
      	if(igcmtchoice.eq.1) then
      		if (i.eq.1) filein=allorder
      		if (i.eq.2) filein=qcmt
      	else if(igcmtchoice.eq.2) then
      	    if (i.eq.1) filein=allorder
            if (i.eq.2) goto 100
      	else if(igcmtchoice.eq.3) then
      	    if (i.eq.1) filein=qcmt
            if (i.eq.2) goto 100
		endif
c---- read in the nbn file
      	inquire(file=filein,exist=exists)
      	if(.not.exists) then
			write(6,"('catalog does not exist:',a)") filein(1:lnblnk(filein))
			ierror=1
			return
      	endif
c
c---- read in the nbn file
c
      	open(14,file=filein,access='direct',form='unformatted',
     #      recl=360)
      	call readnbn(14,nread,ierror)
      	close(14)
      	if(ierror.ne.0) then
			ierror=2
			return
     	endif
      enddo !allorder & qcmt
  100 continue
        write(6,"('read ',i6,' CMT events')") nread

      return
      end
