      subroutine sw_hitcount(npath,eplat,eplon,stlat,stlon,xsize,
     #      xminla,xmaxla,xminlo,xmaxlo,hitlat,hitlon,hitcount)
Cf2py intent(inout) npath,eplat,eplon,stlat,stlon,xsize,xminla,xmaxla,xminlo,xmaxlo
Cf2py intent(out) hitlat,hitlon,hitcount
Cf2py depend(xminla,xmaxla,xminlo,xmaxlo) hitlat,hitlon,hitcount
c
      integer npath
      real*4 eplat(npath)
      real*4 eplon(npath)
      real*4 stlat(npath)
      real*4 stlon(npath)
c
c --- segments of the ray paths
c
      real*4 elat,elon,slat,slon,delta,delta2,az2,azatep
      real*4 seglat,seglon,segweight,dseg
      integer iseg,nseg
      real*4 geoco
      parameter (geoco=0.993277)
c
c --- cells
c
      integer icell
      real*4 xsize
      integer mxcell,ncell
      parameter(mxcell=100000)
      integer icount(mxcell)
      integer ihit(mxcell)
      real*4 areaarr(mxcell)
      real*4 xlalo(mxcell),xlahi(mxcell),xlami(mxcell)
      real*4 xlolo(mxcell),xlohi(mxcell),xlomi(mxcell)
      real*4 arealatzone,totarea,area,twopi,areamax
      real*4 xminla,xmaxla,xminlo,xmaxlo
      real*4 xlat,xlon
      parameter (twopi=6.2831853)
      integer ila,ilo,mxla,mxlo
      integer xlen,ylen
      real,dimension(:,:)::
     #      hitcount(int((xmaxlo-xminlo)/xsize),
     #      int((xmaxla-xminla)/xsize)),   !<-  c is allocatable, rank  
     #      hitlat(int((xmaxlo-xminlo)/xsize),
     #      int((xmaxla-xminla)/xsize)),   !<-  c is allocatable, rank  
     #      hitlon(int((xmaxlo-xminlo)/xsize),
     #      int((xmaxla-xminla)/xsize))   !<-  c is allocatable, rank  
      integer, dimension(:,:), allocatable :: icellarr
      integer ipath
c ----------------------------------------
c
c --- define blocks(cells) 
c
      mxla=int(180./xsize)
      mxlo=int(360./xsize)
      allocate(icellarr(mxla,mxlo)) 
      ncell=0
c      nlazone=180/isize
      totarea=0
c      xsize=float(isize)

c --- find centers of the cells

      icell=0
      xlat=xminla-xsize/2.0
      do while(xlat.lt.(xmaxla-xsize))
        xlat=xlat+xsize       
        xlon=xminlo-xsize/2.0
        do while(xlon.lt.(xmaxlo-xsize))
          xlon=xlon+xsize      
          icell=icell+1
          if(icell.gt.mxcell) stop 'icell.gt.mxcell'
          xlami(icell)=xlat
          xlomi(icell)=xlon      
        enddo
      enddo
      ncell=icell
      write(6,"('number of cells= ',i8)")ncell

c--- allocate arrays
c      xlen=int((xmaxlo-xminlo)/float(isize))
c      ylen=int((xmaxla-xminla)/float(isize))
      xlen=int((xmaxlo-xminlo)/xsize)
      ylen=int((xmaxla-xminla)/xsize)

c --- areas      

      totarea=0.0       
      areamax=0.0
      do icell=1,ncell
        xlat=xlami(icell)
        xlon=xlomi(icell)
c----   sind does not work with gfortran
c        arealatzone=sind(xlat+xsize/2.0)-sind(xlat-xsize/2.0)
        arealatzone=sin((xlat+xsize/2.0)*twopi/360.)-
     #      sin((xlat-xsize/2.0)*twopi/360.)
           area=twopi*arealatzone*xsize/360.0
        if(area.gt.areamax) areamax=area
        areaarr(icell)=area
        totarea=totarea+area
      enddo
c      write(6,"('total area= ',f)")totarea

c --- find indices for all cells

      do icell=1,ncell
        xlat=xlami(icell)
        xlon=xlomi(icell)
        xlalo(icell)=xlat-xsize/2.0
        xlahi(icell)=xlat+xsize/2.0
        xlolo(icell)=xlon-xsize/2.0
        xlohi(icell)=xlon+xsize/2.0
      enddo

C---- does not work with 1 degree pixel-ignores equator and mixes up IDT and greenwich meridien
        do icell=1,ncell
             seglat=xlami(icell)
             seglon=xlomi(icell) 
          ila=int(1+((90.-seglat)/xsize))
          !if (ila.eq.91) ila=90
          if (ila.eq.mxla+1) ila=mxla
          ilo=int(1+((180.+seglon)/xsize))
          if (ilo.eq.mxlo+1) ilo=mxlo
          !if (ilo.eq.181) ilo=1
          icellarr(ila,ilo)=icell
          hitlat(ila,ilo)=seglat
          hitlon(ila,ilo)=seglon
        enddo

c --- check if cells are sampled by the rays

      do icell=1,ncell
      icount(icell)=0
      enddo  

      do ipath=1,npath
c
        elat=atan(geoco*tan(eplat(ipath)*twopi/360.))
        elon=eplon(ipath)
        slat=atan(geoco*tan(stlat(ipath)*twopi/360.))
        slon=stlon(ipath)
        call delazgc(elat,elon,slat,slon,delta,azatep,az2)
c
c---- minor arc paths
c
      do icell=1,ncell
        ihit(icell)=0
       enddo
        nseg=int(2.*delta)+1
        dseg=delta/float(nseg)
        segweight=dseg/delta
        do iseg=1,nseg
          delta2=dseg*0.5+float(iseg-1)*dseg
          call pdaz(elat,elon,azatep,delta2,seglat,seglon)
          if(seglon.gt.180.0) seglon=seglon-360.
          if(seglon.lt.-180.0) seglon=seglon+360.
C---- does not work with 1 degree pixel-ignores equator and mixes up IDT and greenwich meridien
          ila=int(1+((90.-seglat)/xsize))
          !if (ila.eq.91) ila=90
          if (ila.eq.mxla+1) ila=mxla
          ilo=int(1+((180.+seglon)/xsize))
          if (ilo.eq.mxlo+1) ilo=mxlo
          !if (ilo.eq.181) ilo=1

        if(ila.lt.1.or.ila.gt.mxla) then
           ierror = 2
           return
        endif
        if(ilo.lt.1.or.ilo.gt.mxlo) then
           ierror=3
           return
        endif
          icell=icellarr(ila,ilo)
        ihit(icell)=ihit(icell)+1
        enddo ! --- iseg
      do icell=1,ncell
        if(ihit(icell).ge.1) then
          icount(icell)=icount(icell)+1            
        endif
       enddo
c
        if(mod(ipath,5000).eq.0) print*,'.... Processed ipath: ',ipath
      enddo ! --- ipath
c
c --- write out hitcounts for the minor arc paths
c
      do icell=1,ncell
          xlat=xlami(icell)
          xlon=xlomi(icell)
          if(xlon.lt.0.0) xlon=xlon+360.0
          ila=int(1+((90.-xlat)/xsize))
          !if (ila.eq.91) ila=90
          if (ila.eq.mxla+1) ila=mxla
          ilo=int(1+((180.+xlon)/xsize))
          if (ilo.eq.mxlo+1) ilo=mxlo
          !if (ilo.eq.181) ilo=1
          icellarr(ila,ilo)=icell
          hitlat(ila,ilo)=xlat
          hitlon(ila,ilo)=xlon      
          hitcount(ila,ilo)=float(icount(icell))*areamax/areaarr(icell)
      enddo
c      close(luo)
       return
      end
