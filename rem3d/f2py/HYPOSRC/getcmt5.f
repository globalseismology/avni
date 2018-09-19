      subroutine getcmt5(ievt,
     # catalog,iyear,month,iday,ihour,minute,fsec,
     # catlat,catlon,catdep,xmb,xms,region,
     # event,stamp,elat,elon,edep,xm,sc,iexp,
     # ierror)
c
c---- this subroutine returns cmt parameters for a given
c---- event by copying them from common
c     getcmt5 returns several parameters for each earthquake in the catalog:
c
c     catalog (character*4) source for the hypcentral information (usually PDE)
c     iyear (integer*4) year of earthquake in catalog
c     month (integer*4) month of earthquake
c     iday  (integer*4) day
c     minute (integer*4) minute
c     fsec (real*4) second
c     catlat (real*4) reported latitude
c     catlon (real*4) reported longitude
c     catdep (real*4) reported depth
c     xmb (real*4) reported body-wave magnitude
c     xms (real*4) reported surface-wave magnitude
c     region (character*24) reported region
c     event (character*16) CMT catalog name of event
c     stamp (character*16) timestamp for CMT analysis
c     elat (real*4) centroid latitude
c     elon (real*4) centroid longitude
c     edep (real*4) centroid depth
c     xm(6) (real*4) array with moment tensor elements
c     sc (real*4) scalar seismic moment
c     iexp (integer*4) exponent to be used to make full scalar moment, m0=sc*10**iexp 
c     ierror (integer*4) error flag (0 ==> no error) 
c     xmw-MOMENT MAGNITUDE; xmwmin-Minimum xmw in selected catalog

c
      character*4 catalog
      character*16 event
      character*16 stamp
      character*24 region
      dimension xm(6)
      common /plevel/ iprtlv
c
      include 'nbninfo.h'
c
      ierror=0
      if(ievt.gt.nreadnbn) then
        ierror=1
      else
        catalog=catalognbn(ievt)
        iyear=lyearnbn(ievt)
        month=monthnbn(ievt)
        iday=idaynbn(ievt)
        ihour=ihnbn(ievt)
        minute=minnbn(ievt)
        fsec=fsecnbn(ievt)
        catlat=eplatnbn(ievt)
        catlon=eplongnbn(ievt)
        catdep=depthnbn(ievt)
        xmb=xmbnbn(ievt)
        xms=xmsnbn(ievt)
        region=regionnbn(ievt)
        event=cmtnamenbn(ievt)
        stamp=stampnbn(ievt)
        elat=epanbn(ievt)
        elon=eponbn(ievt)
        edep=xdnbn(ievt)
        do i=1,6
          xm(i)=xmnbn(i,ievt)
        enddo
        sc=scnbn(ievt)
        iexp=iexpnbn(ievt)
      endif
      return
      end
