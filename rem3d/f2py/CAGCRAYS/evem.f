      subroutine evem(ilay,rad,isotropic,rho,vpv,vph,vsv,vsh,
     #eta,qmu,qka,ierr)
      implicit double precision (a-h,o-z)
      logical isotropic
      include 'emcommon.h'
      ifilev=ibotlev(ilay)
      ilalev=itoplev(ilay)
      if(rad-radlev(itoplev(ilay)).gt.0.001d0.or.
     #   rad-radlev(ibotlev(ilay)).lt.-0.001d0) then
        write(6,"('extrapolating:',i4,3f10.2)") 
     #    ilay,rad,radlev(ibotlev(ilay)),radlev(itoplev(ilay))
	ierr=1
      endif
      rho=drsple(ifilev,ilalev,radlev,rholev,rhospl,rad)
      qka=drsple(ifilev,ilalev,radlev,bigqklev,bigqkspl,rad)
      qmu=drsple(ifilev,ilalev,radlev,bigqmlev,bigqmspl,rad)
      if(isotropic) then
        vpv=drsple(ifilev,ilalev,radlev,vpisolev,vpisospl,rad)
        vsv=drsple(ifilev,ilalev,radlev,vsisolev,vsisospl,rad)
        vph=vpv
        vsh=vsv
        eta=1.d0
      else
        vpv=drsple(ifilev,ilalev,radlev,vpvlev,vpvspl,rad)
        vsv=drsple(ifilev,ilalev,radlev,vsvlev,vsvspl,rad)
        vph=drsple(ifilev,ilalev,radlev,vphlev,vphspl,rad)
        vsh=drsple(ifilev,ilalev,radlev,vshlev,vshspl,rad)
        eta=drsple(ifilev,ilalev,radlev,etalev,etaspl,rad)
      endif
      ierr=0
      return
      end
