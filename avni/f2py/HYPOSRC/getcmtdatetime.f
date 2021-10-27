      subroutine getcmtdatetime(searchevent,igcmtchoice,datetime)
c
c---- this subroutine returns cmt date in the 2014-01-08 00:00:13 format for use in python
c     igmtchoice is 1 if both monthly+qcmt are used (allorder.nbn & qcmt.nbn)
c                   2 if only monthly are used (allorder.nbn)
c                   3 if qcmt is used
c
      logical exists
c---- for getcmtraj
      character*4 catalog
      character*16 event,searchevent
      character*16 stamp
      character*24 region
      character*19 datetime
      dimension xm(6), xerr(6), ev(3), ipl(3), iaz(3)
      dimension istr(2), idip(2), islp(2)
      character*80 seismdir,filein
	  real*4 areal
c      
c
      real*4 labo
      real*4 lato
      real*4 lole
      real*4 lori
      dimension timestart(6)
      dimension timeend(6)
      dimension time(6)
Cf2py intent(inout) searchevent,igcmtchoice
Cf2py intent(out) datetime

c
c---- open in the .nbn file
c
      call getenv('SEISM',seismdir)
      
      do i=1,2
      	if(igcmtchoice.eq.1) then
      		if (i.eq.1) filein=seismdir(1:lnblnk(seismdir))//'/allorder.nbn'
      		if (i.eq.2) filein=seismdir(1:lnblnk(seismdir))//'/qcmt.nbn'
      	else if(igcmtchoice.eq.2) then
      	    if (i.eq.1) filein=seismdir(1:lnblnk(seismdir))//'/allorder.nbn'
            if (i.eq.2) goto 100
      	else if(igcmtchoice.eq.3) then
      	    if (i.eq.1) filein=seismdir(1:lnblnk(seismdir))//'/qcmt.nbn'
            if (i.eq.2) goto 100
		endif	
c---- read in the nbn file			
      	inquire(file=filein,exist=exists)
      	if(.not.exists) then
			write(6,"('catalog does not exist:',a)") filein(1:lnblnk(filein))
			stop
      	endif
c
c---- read in the nbn file
c
      	open(14,file=filein,access='direct',form='unformatted',
     #      recl=360)
      	call readnbn(14,nread)
      	close(14)
C     	write(6,"('read ',i6,' CMT events')") nread
c
c---- open the file and parse it -- loop over all events in the catalog

	   call getcmtbyname(searchevent,ievt,
     # iyear,month,iday,ihour,minute,fsec,
     # elat,elon,edep,xmw,ierror)
	   if (ierror.eq.0) then
	   	  write(datetime,"(i4.4,'-',i2.2,'-',i2.2,' ',i2.2,':',i2.2,':',i2.2)") iyear,month,iday,ihour,minute,int(fsec)
C   	  write(6,'(i2.2)') int(fsec)
	      goto 100
	   endif 
      enddo !allorder & qcmt
  100 continue
c-----finished selection of events
      
      
      if (ierror.eq.1) write(6,"('Warning: No earthquakes found in Global CMT catalog in getcmtdate for ',a)") searchevent(1:lnblnk(searchevent))
        
      return  
      end
