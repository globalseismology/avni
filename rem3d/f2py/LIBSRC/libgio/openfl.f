      subroutine openfl(lufl,namein,iap,ifile,irec,istat,lrec)
c
c---- ifile is not used
c---- irec is the first record to read
c---- lrec is record length, if negative => a byte stream
c
      include 'giolib.h'
      character*(*) namein
      logical exists
      integer*4 ifarr(13)
      integer*4 lrec
      integer*4 irec
      integer fstat
c
      istat=0
      if(lufl.gt.0.and.lufl.le.20) then
      else
	stop 'lufl outside range'
      endif
c
      lup=1000+lufl
      inquire(file=namein,exist=exists)
      if(.not.exists) then
	stop 'openfl: file does not exist' 
      endif
c
      if(iopen(lufl).gt.0) then
	stop 'openfl: logical unit already in use'
      endif
c
      open(lup,file=namein,access='stream',form='unformatted')
      iopen(lufl)=1
      irw(lufl)=1
      if(iap.eq.4) then
	irw(lufl)=4
      endif
      irecl(lufl)=lrec
      ipos(lufl)=irec*abs(irecl(lufl))
c
      isize(lufl)=0
      ii=fstat(lup,ifarr)
      isize(lufl)=ifarr(8)
c      write(6,"('file size:',i10)") isize(lufl)
      istat=0
      return
      end
