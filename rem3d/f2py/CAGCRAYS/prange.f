      subroutine prange(name,numlay,xbotlay,xtoplay,ianisoflag,ishflag,
     #           radso,isl,icm,ici,pmin,pmax)
c
c---- this subroutine calculates a range of slownesses for which
c---- a given phase may exist. dtm is a minimum bottoming depth
c---- given as input to the routine and is used to calculate a maximum 
c---- slowness -- a given ray has to exist at this depth
c
      common /plevel/ iprtlv
      character*16 name
      character*1 ca,cb
      character*3 c3a
c
      double precision xbotlay(numlay)
      double precision xtoplay(numlay)
      double precision r8
      double precision radso
      double precision vpv,vph,vsv,vsh,rho,qk,qm,eta
      logical isotropic
c
      parameter (radfac=3.1415926534/180.)
c
c----- evaluate wave slownesses at different interfaces
c
      isotropic=.true.
      if(ianisoflag.eq.1) then
        isotropic=.false.
      endif
c
      call evem(isl,radso,isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
      vptm=sngl(vph)
      if(ishflag.eq.1) then
        vstm=sngl(vsh)
      else
        vstm=sngl(vsv)
      endif
c
      r8=xbotlay(icm+1)
      call evem(icm+1,r8,isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
      vpbm=sngl(vph)
      if(ishflag.eq.1) then
        vsbm=sngl(vsh)
      else
        vsbm=sngl(vsv)
      endif
c
      call evem(icm,r8,isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
      vptc=sngl(vpv)
c
      r8=xbotlay(ici+1)
      call evem(ici+1,r8,isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
      vpbc=sngl(vpv)
c
      call evem(ici,r8,isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
      vpti=sngl(vph)
      if(ishflag.eq.1) then
        vsti=sngl(vsh)
      else
        vsti=sngl(vsv)
      endif
c
      pptm=radfac*sngl(radso)/vptm
      pstm=radfac*sngl(radso)/vstm
      ppbm=radfac*sngl(xtoplay(icm))/vpbm
      psbm=radfac*sngl(xtoplay(icm))/vsbm
      pptc=radfac*sngl(xtoplay(icm))/vptc
      ppbc=radfac*sngl(xtoplay(ici))/vpbc
      ppti=radfac*sngl(xtoplay(ici))/vpti
      psti=radfac*sngl(xtoplay(ici))/vsti
      if(iprtlv.gt.0) then
        write(6,"(8f10.3)") pptm,pstm,ppbm,psbm,pptc,ppbc,ppti,psti
      endif
      ioc=0
      iic=0
      inoc=0
c
c---- inserted 10 lines to deal with direct upgoing p and s.
c
      if(name(1:1).eq.'p'.and.name(2:2).eq.' ') then
        pmin=0.0
        pmax=pptm
        go to 50
      endif
      if(name(1:1).eq.'s'.and.name(2:2).eq.' ') then
        pmin=0.0
        pmax=pstm
        go to 50
      endif
c
c----------------------------------------------------------
c
      do 23 i=1,16
      ca=name(i:i)
      if(ca.eq.'K') inoc=inoc+1
      if(ca.eq.'c'.or.ca.eq.'K') ioc=ioc+1
      if(ca.eq.'I'.or.ca.eq.'J'.or.ca.eq.'i') iic=iic+1
   23 continue
c
c---- this line looks to be wrong; changed 7/17/04
ccc      if(iic+ioc.gt.0.and.ishflag.eq.1) then
      if(iic+inoc.gt.0.and.ishflag.eq.1) then
        write(6,"('SH phase cannot interact with core')")
        pmin=0.0
        pmax=0.0
        return
      endif
      if(iic.gt.0) go to 41
      if(ioc.gt.0) go to 31
c
c------ the ray does not interact with the core
c------ ray consists of 'S','s','P','p'.
c
      ip=0
      is=0
      pmax=pstm
      pmin=ppbm
      do 22 i=1,16
      if(i.gt.1)  ca=name(i-1:i-1)
      if(i.eq.1)  ca=' '
      cb=name(i:i)
      if(ca.eq.'p') ip=ip+1
      if(ca.eq.'P'.or.cb.eq.'P') ip=ip+1
cccc      if(ca.eq.'p'.and.cb.eq.'S') ip=ip-1
      if(ca.eq.'S'.or.cb.eq.'S') is=is+1
   22 continue
      if(ip.gt.0) then
          pmax=amin1(pmax,pptm)
          pmin=amax1(pmin,ppbm)
      endif
      if(is.gt.0) then
          pmin=amax1(pmin,psbm)
          pmax=amin1(pmax,pstm)
      endif
      go to 50
   31 continue
c
c------ the ray interacts with the outer core and the mantle
c------ ray consists of 'S','s','p','P','K','c'.
c
      pmax=psbm
      pmin=0.0
      is=0
      do 32 i=2,16
      ca=name(i-1:i-1)
      cb=name(i:i)
      if(i.gt.1)c3a=name(i-1:i+1)
      if(i.eq.1)c3a=' '//name(i:i+1)
      if(ca.eq.'K') pmin=amax1(pmin,ppbc)
      if((ca.eq.'K'.and.cb.eq.'S').or.(cb.eq.'K'.and.ca.eq.'S'))then
         pmax=amin1(pmax,pptc)
      endif
      if((ca.eq.'c'.and.cb.eq.'P').or.(cb.eq.'c'.and.ca.eq.'P'))then
         pmax=amin1(pmax,ppbm)
      endif
      if((ca.eq.'K'.and.cb.eq.'P').or.(cb.eq.'K'.and.ca.eq.'P'))then
         pmax=amin1(pmax,ppbm)
      endif
      if(i.eq.16) go to 32
      if((c3a.eq.' PS').or.(c3a.eq.' PP').or.(c3a.eq.'SPP').or.
     #   (c3a.eq.'PPP').or.(c3a.eq.'PPS').or.(c3a.eq.'SPP').or.
     #   (c3a.eq.'SP ').or.(c3a.eq.'pPS').or.(c3a.eq.'sPS').or.
     #   (c3a.eq.'pPP').or.(c3a.eq.'sPP').or.
     #   (c3a.eq.'SPS')) then
          pmin=amax1(pmin,ppbm)
      endif
      if(ca.eq.'S') is=is+1
      if(cb.eq.'S') is=is+1
   32 continue
      if(is.eq.0) pmax=amin1(pmax,ppbm)
      go to 50
   41 continue
      pmin=0.0
      pmax=ppbc
      do 42 i=1,16
      ca=name(i:i)
      if(ca.eq.'I') pmax=amin1(ppti,pmax)
   42 continue
      go to 50
   50 continue
      if(iprtlv.gt.0) then
        write(6,"('for the phase ',a16,' pmin=',f12.8,' pmax=',f12.8)")
     #        name,pmin,pmax
      endif
      return
      end
