      subroutine sw_hitcount(npath,eplat,eplon,stlat,stlon,isize,xminla,xmaxla,xminlo,xmaxlo,hitlat,hitlon,hitcount)
Cf2py intent(inout) eplat,eplon,stlat,stlon,mxpath,mxcell,npath,xminla,xmaxla,xminlo,xmaxlo
Cf2py intent(out) hitlat,hitlon,hitcount
c
      integer npath
      real, dimension(:), allocatable :: eplat,eplon,stlat,stlon
c
      character*80 pathfile
      character*80 hitfil
      character*80 getunx
      character*256 c256
      logical header
c
      logical exists
      logical eof
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
      integer isize,icell
      real*4 xsize
      integer mxcell,ncell,nlazone
      parameter(mxcell=100000)
      integer icount(mxcell)
      integer ihit(mxcell)
      real*4 areaarr(mxcell)
      real*4 xlalo(mxcell),xlahi(mxcell),xlami(mxcell)
      real*4 xlolo(mxcell),xlohi(mxcell),xlomi(mxcell)
      real*4 arealatzone,totarea,area,twopi,areamax
      real*4 xminla,xmaxla,xminlo,xmaxlo,xlat,xlon,value
      parameter (twopi=6.2831853)
      integer ila,ilo,mxla,mxlo
      parameter(mxla=180)
      parameter(mxlo=360)
      integer xlen,ylen
      integer icellarr(mxla,mxlo)      
      real, dimension(:,:), allocatable :: hitcount,hitlat,hitlon   !<-  c is allocatable, rank  
c
      integer ipath,luo,luin
c ----------------------------------------
c
c --- define blocks(cells) 
c
      ncell=0
      nlazone=180/isize
      totarea=0
      xsize=float(isize)

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
      xlen=int((xmaxlo-xminlo)/float(isize))
      ylen=int((xmaxla-xminla)/float(isize))
      allocate(hitcount(xlen,ylen))
      allocate(hitcount(xlen,ylen))
      allocate(hitcount(xlen,ylen))
      allocate(eplat(npath))
      allocate(eplon(npath))
      allocate(stlat(npath))
      allocate(stlon(npath))
c --- areas      

      totarea=0.0       
      areamax=0.0
      do icell=1,ncell
        xlat=xlami(icell)
        xlon=xlomi(icell)
        arealatzone=sind(xlat+xsize/2.0)-sind(xlat-xsize/2.0)
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
          ila=int(1+((90.-seglat)/isize))
          !if (ila.eq.91) ila=90
          if (ila.eq.mxla+1) ila=mxla
          ilo=int(1+((180.+seglon)/isize))
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
        elat=atand(geoco*tand(eplat(ipath)))
        elon=eplon(ipath)
        slat=atand(geoco*tand(stlat(ipath)))
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
          ila=int(1+((90.-seglat)/isize))
          !if (ila.eq.91) ila=90
          if (ila.eq.mxla+1) ila=mxla
          ilo=int(1+((180.+seglon)/isize))
          if (ilo.eq.mxlo+1) ilo=mxlo
          !if (ilo.eq.181) ilo=1

        if(ila.lt.1.or.ila.gt.mxla) stop 'ila'
        if(ilo.lt.1.or.ilo.gt.mxlo) then
              print*,seglat,seglon,mxla,mxlo,ila,ilo,ipath,isize
              stop 'ilo'
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
        if(mod(ipath,5000).eq.0) write(6,"('.... Proceesed ipath: ',i12,' out of ',i12)") ipath,npath
      enddo ! --- ipath
c
c --- write out hitcounts for the minor arc paths
c
c      write(6,"('writing ...')")
c      lstr=lnblnk(pathfile)
c      hitfil=pathfile(1:lstr)//'.hit'
c      write(6,"('writing to: ',a)") hitfil(1:lnblnk(hitfil))
c      open(luo,file=hitfil(1:lnblnk(hitfil)))
      do icell=1,ncell
c        if(icount(icell).gt.0) then
          xlat=xlami(icell)
          xlon=xlomi(icell)
          if(xlon.lt.0.0) xlon=xlon+360.0
          ila=int(1+((90.-xlat)/isize))
          !if (ila.eq.91) ila=90
          if (ila.eq.mxla+1) ila=mxla
          ilo=int(1+((180.+xlon)/isize))
          if (ilo.eq.mxlo+1) ilo=mxlo
          !if (ilo.eq.181) ilo=1
          icellarr(ila,ilo)=icell
          hitlat(ila,ilo)=xlat
          hitlon(ila,ilo)=xlon      
          hitcount(ila,ilo)=float(icount(icell))*areamax/areaarr(icell)
      
c      write(luo,*)xlat,xlon,value
c      endif
      enddo
c      close(luo)
       return
      end
