      subroutine bffi(lufl,ifbin,ibuf,nbytes,istat,nread,irec)
      include "giolib.h"
      common/plevel/iprtlv
      byte ibuf(4)
c
      nread=0
      istat=9
c
      if(lufl.gt.0.and.lufl.le.20) then
      else
	stop 'lufl outside range'
      endif
c
      if(irec.eq.0) then
	iposition=ipos(lufl)+1
      else
	iposition=(irec-1)*abs(irecl(lufl))+1
      endif
      ntoread=nbytes
      if(irecl(lufl).gt.0) then
        ntoread=min(ntoread,irecl(lufl))
      endif
c
      if(iposition.gt.isize(lufl)) then
	if(iprtlv.gt.0) then
	write(6,"('trying to start read beyond EOF')") 
	endif
	istat=3
	return
      endif
c
      if(iposition+ntoread-1.gt.isize(lufl)) then
	write(6,"('trying to read across EOF')")
	ntoread=isize(lufl)-iposition+1
      endif
c
      lup=lufl+1000
      read(lup,pos=iposition,iostat=ios) (ibuf(i),i=1,ntoread)
      if(ios.eq.0) then
	nread=ntoread
	istat=2
	ipos(lufl)=iposition+nread-1
      else
	write(6,"('ios in bffi:',i8)") ios
      endif
      return
      end

