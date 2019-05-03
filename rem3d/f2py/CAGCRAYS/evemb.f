      subroutine evemb(ilay,rad,rho,grav,bigk,pres,vaisala,bullen,ierr)
      implicit double precision (a-h,o-z)
      logical isotropic
      include 'emcommonb.h'
      ifilev=ibotlev(ilay)
      ilalev=itoplev(ilay)
      if(rad-radlev(itoplev(ilay)).gt.0.001d0.or.
     #   rad-radlev(ibotlev(ilay)).lt.-0.001d0) then
        write(6,"('extrapolating:',i4,3f10.2)") 
     #    ilay,rad,radlev(ibotlev(ilay)),radlev(itoplev(ilay))
	ierr=1
      endif
      rho=drsple(ifilev,ilalev,radlev,rholev,rhospl,rad)
      bigk=drsple(ifilev,ilalev,radlev,bigklev,bigkspl,rad)
      grav=drsple(ifilev,ilalev,radlev,gravlev,gravspl,rad)
      pres=drsple(ifilev,ilalev,radlev,prelev,prespl,rad)
c
c---- calculate the Brunt-Vaisala frquency and Bullen parameter
c
      rhosi=rho*1000.
      rhoder=drspleder(maxlev,ifilev,ilalev,radlev,rholev,rhospl,rad)
      vaisala=-grav*rhoder/rhosi-rhosi*grav**2/bigk
c
      if(vaisala.ne.0..and.grav.ne.0.) then
        bullen=1.0+vaisala*bigk/(rhosi*grav**2)
      else
        bullen=1.0
      endif
      ierr=0
      return
      end
