      parameter (mxnbn=500000)
      common /nbncommon/ nreadnbn,ivernbn,cmtnamenbn,stampnbn,
     #  catalognbn,regionnbn,lyearnbn,monthnbn,idaynbn,ihnbn,minnbn,
     #  fsecnbn,eplatnbn,eplongnbn,depthnbn,xmbnbn,xmsnbn,
     #  itypcmtnbn,itypstfnbn,isbnbn,icbnbn,icutbnbn,
     #  issnbn,icsnbn,icutsnbn,ismnbn,icmnbn,icutmnbn,
     #  ires1nbn,ires2nbn,ires3nbn,ires4nbn,ires5nbn,ires6nbn,
     #  torgnbn,jhnbn,jminnbn,xsecnbn,errtnbn,epanbn,xlatnbn,
     #  insnbn,erranbn,eponbn,xlonnbn,iewnbn,erronbn,xdnbn,
     #  errdnbn,itypdepnbn,durtnbn,iexpnbn,xmnbn,xerrnbn,
     #  evnbn,iplnbn,iaznbn,scnbn,istrnbn,idipnbn,islpnbn
      integer ivernbn(mxnbn)
      character*16 cmtnamenbn(mxnbn)
      character*16 stampnbn(mxnbn)
      character*4 catalognbn(mxnbn)
      character*24 regionnbn(mxnbn)
      integer lyearnbn(mxnbn)
      integer monthnbn(mxnbn)
      integer idaynbn(mxnbn)
      integer ihnbn(mxnbn)
      integer minnbn(mxnbn)
      real fsecnbn(mxnbn)
      real eplatnbn(mxnbn)
      real eplongnbn(mxnbn)
      real depthnbn(mxnbn)
      real xmbnbn(mxnbn)
      real xmsnbn(mxnbn)
      integer itypcmtnbn(mxnbn)
      integer itypstfnbn(mxnbn)
      integer isbnbn(mxnbn)
      integer icbnbn(mxnbn)
      integer icutbnbn(mxnbn)
      integer issnbn(mxnbn)
      integer icsnbn(mxnbn)
      integer icutsnbn(mxnbn)
      integer ismnbn(mxnbn)
      integer icmnbn(mxnbn)
      integer icutmnbn(mxnbn)
      integer ires1nbn(mxnbn)
      integer ires2nbn(mxnbn)
      integer ires3nbn(mxnbn)
      integer ires4nbn(mxnbn)
      integer ires5nbn(mxnbn)
      integer ires6nbn(mxnbn)
      real torgnbn(mxnbn)
      integer jhnbn(mxnbn)
      integer jminnbn(mxnbn)
      real xsecnbn(mxnbn)
      real errtnbn(mxnbn)
      real epanbn(mxnbn)
      real xlatnbn(mxnbn)
      integer insnbn(mxnbn)
      real erranbn(mxnbn)
      real eponbn(mxnbn)
      real xlonnbn(mxnbn)
      integer iewnbn(mxnbn)
      real erronbn(mxnbn)
      real xdnbn(mxnbn)
      real errdnbn(mxnbn)
      integer itypdepnbn(mxnbn)
      real durtnbn(mxnbn)
      integer iexpnbn(mxnbn)
      real xmnbn(6,mxnbn)
      real xerrnbn(6,mxnbn)
      real evnbn(3,mxnbn)
      integer iplnbn(3,mxnbn)
      integer iaznbn(3,mxnbn)
      real scnbn(mxnbn)
      integer istrnbn(2,mxnbn)
      integer idipnbn(2,mxnbn)
      integer islpnbn(2,mxnbn)
      integer nreadnbn
      save
