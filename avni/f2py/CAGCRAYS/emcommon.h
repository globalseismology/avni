      character*80 emtitle
c
      common /empar1/ emtitle
c
      parameter (maxlev=1000)
      real*8 radlev(maxlev)
      real*8 rholev(maxlev)
      real*8 vpvlev(maxlev)
      real*8 vsvlev(maxlev)
      real*8 bigqklev(maxlev)
      real*8 bigqmlev(maxlev)
      real*8 vphlev(maxlev)
      real*8 vshlev(maxlev)
      real*8 etalev(maxlev)
      real*8 vsisolev(maxlev)
      real*8 vpisolev(maxlev)
      real*8 elllev(maxlev)
      real*8 eta2lev(maxlev)
      real*8 gravlev(maxlev)
c
      real*8 rhospl(3,maxlev)
      real*8 vpvspl(3,maxlev)
      real*8 vsvspl(3,maxlev)
      real*8 bigqkspl(3,maxlev)
      real*8 bigqmspl(3,maxlev)
      real*8 vphspl(3,maxlev)
      real*8 vshspl(3,maxlev)
      real*8 etaspl(3,maxlev)
      real*8 vpisospl(3,maxlev)
      real*8 vsisospl(3,maxlev)
      real*8 ellspl(3,maxlev)
      real*8 eta2spl(3,maxlev)
      real*8 gravspl(3,maxlev)
c
      real*8 work(3,maxlev)
c      
      common /empar2/ radlev,
     #       rholev,rhospl,
     #       vpvlev,vpvspl,
     #       vsvlev,vsvspl,
     #       bigqklev,bigqkspl,
     #       bigqmlev,bigqmspl,
     #       vphlev,vphspl,
     #       vshlev,vshspl,
     #       etalev,etaspl,
     #       vpisolev,vpisospl,
     #       vsisolev,vsisospl,
     #       elllev,ellspl,
     #       eta2lev,eta2spl,
     #       gravlev,gravspl,
     #       numemlev
c
      parameter (maxemreg=30)
      real*8 topemreg(maxemreg)
      real*8 botemreg(maxemreg)
      integer itoplev(maxemreg)
      integer ibotlev(maxemreg)
c
      common /empar3/ topemreg,
     #       botemreg,
     #       itoplev,
     #       ibotlev,
     #       numemreg
     
