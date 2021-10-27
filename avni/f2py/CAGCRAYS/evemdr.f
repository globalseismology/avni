      subroutine evemdr(ilay,rad,isotropic,rho,vpv,vph,vsv,vsh,
     #eta,qmu,qka)
c
c---- this returns the radial derivative of the different variables
c
      implicit double precision (a-h,o-z)
      logical isotropic
      include 'emcommon.h'
      ifilev=ibotlev(ilay)
      ilalev=itoplev(ilay)
      if(rad.gt.radlev(itoplev(ilay)).or.rad.lt.radlev(ibotlev(ilay))) then
        write(6,"('extrapolating:',i4,3f10.2)") ilay,rad,radlev(ibotlev(ilay)),
     #        radlev(itoplev(ilay))
      endif
      rho=drspledr(ifilev,ilalev,radlev,rholev,rhospl,rad)
      qka=drspledr(ifilev,ilalev,radlev,bigqklev,bigqkspl,rad)
      qmu=drspledr(ifilev,ilalev,radlev,bigqmlev,bigqmspl,rad)
      if(isotropic) then
        vpv=drspledr(ifilev,ilalev,radlev,vpiso,vpisospl,rad)
        vsv=drspledr(ifilev,ilalev,radlev,vsiso,vsisospl,rad)
        vph=vpv
        vsh=vsv
        eta=0.d0
      else
        vpv=drspledr(ifilev,ilalev,radlev,vpvlev,vpvspl,rad)
        vsv=drspledr(ifilev,ilalev,radlev,vsvlev,vsvspl,rad)
        vph=drspledr(ifilev,ilalev,radlev,vphlev,vphspl,rad)
        vsh=drspledr(ifilev,ilalev,radlev,vshlev,vshspl,rad)
        eta=drspledr(ifilev,ilalev,radlev,etalev,etaspl,rad)
      endif
      ierr=0
      return
      end
