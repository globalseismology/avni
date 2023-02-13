      subroutine getcmtbyname(event,ievt,
     # iyear,month,iday,ihour,minute,fsec,
     # elat,elon,edep,xmw,ierror)

c
c---- this subroutine returns cmt parameters for a given
c---- event by copying them from common
c
Cf2py intent(inout) event
Cf2py intent(out) ievt,iyear,month,iday,ihour,minute,fsec,elat,elon,edep,xmw,ierror

       real fsec,elat,elon,edep,xmw
       integer ievt,iyear,month,iday,ihour,minute,ierror,iprtlv
      character*16 event
      character*16 ch16
      real*8 dadd
      common /plevel/ iprtlv
c
      include 'nbninfo.h'
c
      levent=lnblnk(event)
      ierror=1
      ievt=0
      do while(ierror.eq.1.and.ievt.lt.nreadnbn)
	ievt=ievt+1
	ch16=cmtnamenbn(ievt)
	lch16=lnblnk(ch16)
	if(lch16.eq.levent) then
	  if(event(2:lch16).eq.ch16(2:lch16)) then
	    ierror=0
            iyear=lyearnbn(ievt)
            month=monthnbn(ievt)
            iday=idaynbn(ievt)
	    julian=julday(iyear,month,iday)
            ihour=ihnbn(ievt)
            minute=minnbn(ievt)
            fsec=fsecnbn(ievt)
	    dadd=dble(torgnbn(ievt))
	    call tadder(iyear,julian,ihour,minute,fsec,i1,i2,i3,i4,f1,dadd)
	    call monday(i1,i2,month,iday)
	    iyear=i1
	    ihour=i3
	    minute=i4
	    fsec=f1
            elat=epanbn(ievt)
            elon=eponbn(ievt)
            edep=xdnbn(ievt)
            iexp=iexpnbn(ievt)
	    sc=scnbn(ievt)
	    xmw=0.6666667*alog10(sc*10.**iexp)-10.7333333
          endif
        endif
      enddo
      return
      end
