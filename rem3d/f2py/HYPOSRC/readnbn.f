c---- subroutine to read CMT catalog in binary nbn format
c---- modified 6/11/2009 to swap byte order depending on endianess
c
      subroutine readnbn(lu,nread)
      character*16 DEVENT
      character*16 timestamp
      character*24 region
      character*4 catalog
      DIMENSION XM(6),XERR(6),EV(3),IPL(3),IAZ(3),
     1   ISTR(2),IDIP(2),ISLP(2) 
c
c---- include file for storing the whole catalog
c
      include 'nbninfo.h'
c
      ios=0
      irec=1
      do while(ios.eq.0)
        read(lu,rec=irec,err=99) iversion,DEVENT,timestamp,catalog,
     1    LYEAR,MONTH,IDAY,IH,MIN,FSEC,EPLAT,
     1    EPLONG,DEPTH,XMB,XMS,region,
     #    itypcmt,itypstf,ISB,ICB,ICUTB,
     #    iss,ics,icuts,
     #    ISM,ICM,ICUTM,
     #    ires1,ires2,ires3,ires4,ires5,ires6,
     2    TORG,JH,JMIN,XSEC,ERRT,EPA,XLAT,INS,ERRA,EPO,XLON,IEW,ERRO,XD,
     3    ERRD,itypdep,DURT,IEXP,(XM(I),I=1,6),(XERR(I),I=1,6),(EV(I),I=1,3),
     4    (IPL(I),I=1,3),(IAZ(I),I=1,3),SC,(ISTR(I),I=1,2),(IDIP(I),I=1,2),
     5    (ISLP(I),I=1,2)
c
c---- check whether the bytes need to be swapped (use the month variable)
c
	if(irec.eq.1) then
	  if(month.gt.0.and.month.lt.13) then
	    iswapend=0
          else
	    idummy=month
	    call swap4of4(idummy)
	    if(idummy.gt.0.and.idummy.lt.13) then
              iswapend=1
c              write(6,"('I will swap bytes reading -.nbn file')")
            else
	      write(6,"('in readnbn -- file neither big or little endian')")
              stop
            endif
          endif
        endif
c
        if(iswapend.eq.1) then
	  call swap4of4(iversion)
	  call swap4of4(lyear)
	  call swap4of4(month)
	  call swap4of4(iday)
	  call swap4of4(ih)
	  call swap4of4(min)
	  call swap4of4(fsec)
	  call swap4of4(eplat)
	  call swap4of4(eplong)
	  call swap4of4(depth)
	  call swap4of4(xmb)
	  call swap4of4(xms)
	  call swap4of4(itypcmt)
	  call swap4of4(itypstf)
	  call swap4of4(isb)
	  call swap4of4(icb)
	  call swap4of4(icutb)
	  call swap4of4(iss)
	  call swap4of4(ics)
	  call swap4of4(icuts)
	  call swap4of4(ism)
	  call swap4of4(icm)
	  call swap4of4(icutm)
	  call swap4of4(ires1)
	  call swap4of4(ires2)
	  call swap4of4(ires3)
	  call swap4of4(ires4)
	  call swap4of4(ires5)
	  call swap4of4(ires6)
	  call swap4of4(torg)
	  call swap4of4(jh)
	  call swap4of4(jmin)
	  call swap4of4(xsec)
	  call swap4of4(errt)
	  call swap4of4(epa)
	  call swap4of4(xlat)
	  call swap4of4(ins)
	  call swap4of4(erra)
	  call swap4of4(epo)
	  call swap4of4(xlon)
	  call swap4of4(iew)
	  call swap4of4(erro)
	  call swap4of4(xd)
	  call swap4of4(errd)
	  call swap4of4(itypdep)
	  call swap4of4(durt)
	  call swap4of4(iexp)
	  do i=1,6
	    call swap4of4(xm(i))
	    call swap4of4(xerr(i))
          enddo
	  do i=1,3
	    call swap4of4(ev(i))
	    call swap4of4(ipl(i))
	    call swap4of4(iaz(i))
	  enddo
	  call swap4of4(sc)
	  do i=1,2
	    call swap4of4(istr(i))
	    call swap4of4(idip(i))
	    call swap4of4(islp(i))
	  enddo
	endif
c
        ivernbn(irec)=iversion
        cmtnamenbn(irec)=devent
        stampnbn(irec)=timestamp
        catalognbn(irec)=catalog
        regionnbn(irec)=region
        lyearnbn(irec)=lyear
        monthnbn(irec)=month
        idaynbn(irec)=iday
        ihnbn(irec)=ih
        minnbn(irec)=min
        fsecnbn(irec)=fsec
        eplatnbn(irec)=eplat
        eplongnbn(irec)=eplong
        depthnbn(irec)=depth
        xmbnbn(irec)=xmb
        xmsnbn(irec)=xms
        itypcmtnbn(irec)=itypcmt
        itypstfnbn(irec)=itypstf
        isbnbn(irec)=isb
        icbnbn(irec)=icb
        icutbnbn(irec)=icutb
        issnbn(irec)=iss
        icsnbn(irec)=ics
        icutsnbn(irec)=icuts
        ismnbn(irec)=ism
        icmnbn(irec)=icm
        icutmnbn(irec)=icutm
        ires1nbn(irec)=ires1
        ires2nbn(irec)=ires2
        ires3nbn(irec)=ires3
        ires4nbn(irec)=ires4
        ires5nbn(irec)=ires5
        ires6nbn(irec)=ires6
        torgnbn(irec)=torg
        jhnbn(irec)=jh
        jminnbn(irec)=jmin
c
        xsecnbn(irec)=xsec
        errtnbn(irec)=errt
        epanbn(irec)=epa
        xlatnbn(irec)=xlat
        insnbn(irec)=ins
        erranbn(irec)=erra
        eponbn(irec)=epo
        xlonnbn(irec)=xlon
        iewnbn(irec)=iew
        erronbn(irec)=erro
        xdnbn(irec)=xd
        errdnbn(irec)=errd
        itypdepnbn(irec)=itypdep
        durtnbn(irec)=durt
        iexpnbn(irec)=iexp
        do i=1,6
          xmnbn(i,irec)=xm(i)
          xerrnbn(i,irec)=xerr(i)
        enddo
        do i=1,3
          evnbn(i,irec)=ev(i)
          iplnbn(i,irec)=ipl(i)
          iaznbn(i,irec)=iaz(i)
        enddo
        scnbn(irec)=sc
        do i=1,2
          istrnbn(i,irec)=istr(i)
          idipnbn(i,irec)=idip(i)
          islpnbn(i,irec)=islp(i)
        enddo
        irec=irec+1
      enddo
   99 continue
      nread=irec-1
      nreadnbn=nread
      return
      end
