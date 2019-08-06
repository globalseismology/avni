      subroutine cagcrays_pm(lu,filename,epla,eplo,depth,stla,stlo,
     #       phase,ishflag,ianisoflag,itypcalc,delta,azep,azst,
     #       maxarr,narr,arrivaltime,slowness,dtddepth,ddeltadp,tstar,
     #       ellcorr,maxraypts,nrayptsx,iraytpx,iraylayx,
     #       rayradx,raydeltax,raylatx,raylongx,raytimex,ireadagain)
c-----------------------------------------------------------------------
c---- This subroutine calculates rays and ray parameters in a geocentric
c---- coordinate system.
c---- lu                 logical unit used for table file
c---- filename           name of travel time table
c---- epla,eplo,depth    source location (geocentric lat)
c---- stla,stlo          receiver location (geocentric lat)
c---- phase              name of phase
c---- ishflag            flag for sh phase (pure S phases) (0-SV)
c---- ianisoflag         flag for anisotropic calculation (0-isotropic)
c---- itypcalc           flag to determine calculation
c                        1 - calculate travel times, delta, dddp, dtdh
c                        2 - calculate ellipticity correction
c                        3 - calculate and return ray coordinates
c---- delta,azep,azst    delta, azimuth at epicenter and at station
c---- maxarr             maximum number of arrivals (in calling program)
c---- narr               number of arrivals returned
c---- arrivaltime        array with arrival times
c---- slowness           array with slowness values
c---- dtddepth           derivative of travel time w.r.t. depth
c---- ddeltadp           derivative of delta w.r.t. ray parameter
c---- tstar              array with tstar
c---- ellcorr            total ellipticity correction (to be added)
c---- maxraypts          maximum number of points in ray description
c---- nrayptsx           number of points in ray description
c---- rayradx            ray coordinate, radius
c---- raylatx            ray coordinate, latitude
c---- raylongx           ray coordinate, longitude
c---- raytimex           ray coordinate, time
c---- iraytpx            ray index (wave type, up or down)
c---- iraylayx           ray layer
c---- ireadagain		 ignores the last time cagcrays was called and reads reads files again
c-----------------------------------------------------------------------
c
      common/plevel/iprtlv
c
      character*80 filename
      character*16 phase
c
      dimension arrivaltime(maxarr)
      dimension slowness(maxarr)
      dimension ddeltadp(maxarr)
      dimension dtddepth(maxarr)
      dimension tstar(maxarr)
      dimension turad(100)
      dimension ellcorr(maxarr)
c
      parameter (mxa=10)
      parameter (mxpts=5000)
      dimension nraypts(mxa)
      dimension rayrad(mxpts,mxa)
      dimension raydelta(mxpts,mxa)
      dimension raylat(mxpts,mxa)
      dimension raylong(mxpts,mxa)
      dimension raytime(mxpts,mxa)
      dimension iraytp(mxpts,mxa)
      dimension iraylay(mxpts,mxa)
c
      dimension nrayptsx(maxarr)
      dimension rayradx(maxraypts,maxarr)
      dimension raydeltax(maxraypts,maxarr)
      dimension raylatx(maxraypts,maxarr)
      dimension raylongx(maxraypts,maxarr)
      dimension raytimex(maxraypts,maxarr)
      dimension iraytpx(maxraypts,maxarr)
      dimension iraylayx(maxraypts,maxarr)
c
      double precision dtor8
      double precision rtod8
      double precision twopi
      parameter ( twopi=6.2831853072d0 )
c      parameter ( dtor8=6.2831853072d0/360.0d0)
c      parameter ( rtod8=360.0d0/6.2831853072d0)
c
      parameter (maxlay=25)
      double precision xbotlay(maxlay)
      double precision xtoplay(maxlay)
      dimension ianiso(maxlay)
      dimension numlev(maxlay)
      dimension ifluid(maxlay)
      dimension idepl(maxlay)
c
      double precision rad8
      integer iss,irr,nn
      dimension np(maxlay),ns(maxlay),id(10,maxlay)
      dimension icn(maxlay)
      dimension iseq(1000)
      dimension idirseg(1000)
c
c      parameter (maxps=14000)
c      parameter (maxps=20500)
      parameter (maxps=16000)
      dimension p(maxps)
      dimension t(maxps),d(maxps),dddp(maxps),q(maxps),turn(maxps)
c
      logical isotropic
      double precision pray8
      double precision xbot8
      double precision xtop8
      double precision delta8(maxlay,3)
      double precision ttime8(maxlay,3)
      double precision dddp8(maxlay,3)
      double precision tstar8(maxlay,3)
      double precision xturn8(maxlay,3)
      parameter (maxray=1000)
      double precision rayrad8(maxray,maxlay,3)
      double precision raydelta8(maxray,maxlay,3)
      double precision raytime8(maxray,maxlay,3)
      dimension nray(maxlay,3)
c
      double precision rho8,vpv8,vph8,vsv8,vsh8,eta8,qka8,qmu8,r8
      double precision eta28,ell8
      logical exists
c

      character*80 filenamelast
      character*16 phaselast
      data filenamelast /' '/
      data phaselast /' '/
      data depthlast /-1.0/
      data ianisoflaglast /-1/
      data ishflaglast /-1/
c
      save
c
c-----------------------------------------------------------------------
c
      ! degree to radian and vice versa
      dtor8=twopi/360.0d0
      rtod8=360.0d0/twopi

	  if(ireadagain.eq.1) then
	  	filenamelast=' '
		phaselast=' '
		depthlast=-1.0
        ianisoflaglast=-1
        ishflaglast=-1
        call reset_ptab
      endif

      narr=0
      if(maxraypts.lt.1) then
        write(6,"('maxraypts must be .ge. 1 calling cagcrays')")
        stop
      endif
      if(maxarr.lt.1) then
        write(6,"('maxarr must be .ge. 1 calling cagcrays')")
        stop
      endif
c
c-----------------------------------------------------------------------
c
      if(itypcalc.le.0.or.itypcalc.gt.3) then
        write(6,"('cagcrays_pm: itypcalc.le.0.or.itypcalc.gt.3')")
        stop
      endif
c
c-----------------------------------------------------------------------
c
      if(itypcalc.ge.0) then
        call delazgc(epla,eplo,stla,stlo,delta,azep,azst)
      endif
c-----------------------------------------------------------------------
c
      if(itypcalc.ge.1) then
        lstring=lnblnk(filename)
        lstringlast=lnblnk(filenamelast)
        if(filename(1:lstring).eq.filenamelast(1:lstringlast)) then
          if(iprtlv.gt.1) then
            write(6,"('Travel time tables are already initialized')")
          endif
        else
          if(filenamelast(1:lstringlast).eq.' ') then
            if(ireadagain.eq.1) call closfl(lu,ierr)
          else
            call closfl(lu,ierr)
          endif
          inquire(file=filename,exist=exists)
		  if(.not.exists) then
	  		write(6,"(a,' does not exists')")filename(1:lnblnk(filename))
	  	  	stop
		  endif

          call openfl(lu,filename,1,0,0,ierr,-1)
          call reademfl(lu,maxlay,numlay,xbotlay,xtoplay,icm,ici,
     #       numlev,ifluid,ianiso,ierr)
          call getflpos(lu,initpos,ierr)
         do i=1,numlay
            idepl(i)=numlay+1-i
          enddo
          if(iprtlv.gt.0) then
            do i=1,numlay
              write(6,"(i3,f10.1,f10.1,3i4)") i,xbotlay(i),xtoplay(i),
     #                                  ifluid(i),ianiso(i),numlev(i)
            enddo
            write(6,"('inner core:',2i5)") ici,idepl(ici)
            write(6,"('outer core:',2i5)") icm,idepl(icm)
          endif
          filenamelast=filename
          depthlast=-1.0
          ianisoflaglast=-1
          ishflaglast=-1
        endif
c
c---- determine the source layer
c
        if(depth.eq.depthlast) then
        else
          rad8=xtoplay(numlay)-dble(depth)
          isl=numlay
          do ilay=numlay,1,-1
            if(rad8.le.xtoplay(ilay)) then
              isl=ilay
            endif
          enddo
          if(iprtlv.gt.0) then
            write(6,"('source layer:',2i5)") isl,idepl(isl)
          endif
        endif
c
c---- find the 670, 400, and 220-km discontinuities
c
        do i=1,numlay
	  if(dabs(xtoplay(i)-(6371.d0-670.d0)).lt.25.d0) then
	    i670=numlay-i+1
	  else if(dabs(xtoplay(i)-(6371.d0-400.d0)).lt.25.d0) then
    	    i400=numlay-i+1
	  else if(dabs(xtoplay(i)-(6371.d0-220.d0)).lt.25.d0) then
	    i220=numlay-i+1
	  endif
        enddo
c
c---- call psvrayin to get ray definition
c
        if(depth.eq.depthlast.and.phase.eq.phaselast.and.
     #     ianisoflag.eq.ianisoflaglast.and.ishflag.eq.ishflaglast) then
        else

          call psvrayin(phase,idepl(isl),idepl(icm),idepl(ici),i670,i400,i220,
     #       numlay,np,ns,id,iss,irr,nn)
          if(iprtlv.gt.1) then
            write(6,"(1x,'ray definition - source layer ',i2,/,
     #       1x,'                 outer core   ',i2,/,
     #       1x,'                 inner core   ',i2,/,
     #       1x,'source index iss              ',i2,/,
     #       1x,'receiver index irr            ',i2,/,
     #       1x,'number of layers for ray      ',i2)") isl,icm,ici,iss,irr,nn
            write(6,"(1x,'layer p-legs s-legs p1p1 p1p2 p2p2 s1s1 s1s2 s2s2',
     *        ' p1s1 p1s2 p2s1 p2s2')")
            do i=1,nn
              write(6,"(1x,i4,i7,i7,10i5)") i,np(i),ns(i),(id(j,i),j=1,10)
            enddo
          endif
c
c---- call rayseq to get a sequencing of the ray
c
          do i=1,numlay
            icn(i)=i
          enddo
          call rayseq(np,ns,numlay,iss,irr,idepl(isl),1,icn,iseq,idirseg,nseg)
          if(iprtlv.gt.1) then
            do i=1,nseg
              write(6,"(i4,i3)") i,iseq(i)
            enddo
          endif
c
c---- call readptab to set up table components for this source depth
c
          rad4=sngl(rad8)
          call readptab(lu,initpos,rad4,isl,ianisoflag,ierr)
c
c---- call prange to find out the range of p values for a given phase
c
          call prange(phase,numlay,xbotlay,xtoplay,ianisoflag,ishflag,
     #       rad8,isl,icm,ici,pmin,pmax)
c
c---- call getptab to create the phase-specific tables
c
          call getptab(ianisoflag,ishflag,
     #       iss,isl,np,ns,id,pmin,pmax,npps,p,t,d,dddp,q,turn)
        endif
        depthlast=depth
        phaselast=phase
        ianisoflaglast=ianisoflag
        ishflaglast=ishflag
c
c---- call rdrays to find the arrival times for the distance
c
        call fdrays(delta,npps,p,t,d,dddp,q,turn,slowness,
     #      arrivaltime,ddeltadp,tstar,turad,maxarr,narr)
        if(iprtlv.gt.1) then
        do iarr=1,narr
          write(6,"('p,turning radius:',2f10.4)") slowness(iarr),turad(iarr)
        enddo
        endif
c
c---- create diffracted P and S phases
c
        idiffphase=0
        diffract=0.
        if(narr.eq.0) then
          if(delta.gt.80.d0) then
            if(phase(1:2).eq.'P '.or.
     #         phase(1:3).eq.'pP '.or.phase(1:3).eq.'sP ') then
               idiffphase=2
            endif
            if(phase(1:2).eq.'S '.or.
     #         phase(1:3).eq.'sS '.or.phase(1:3).eq.'pS ') then
               idiffphase=1
            endif
            if(idiffphase.ne.0) then
               reach=d(1)*sngl(rtod8)
               slowness(1)=p(1)
               ddeltadp(1)=dddp(1)*(rtod8**2)
               tstar(1)=q(1)
               diffract=delta-reach
               arrivaltime(1)=t(1)+diffract*slowness(1)
               if(iprtlv.gt.0) then
               write(6,"('delta,reach',3g15.5)") delta,reach,slowness(1)
               endif
               narr=1
            endif
          endif
        endif
        do iarr=1,narr
          nrayptsx(iarr)=0
          ellcorr(iarr)=0
        enddo
      endif
c
c-----------------------------------------------------------------------
c
      if(itypcalc.gt.1) then
        do iarr=1,narr
          pray8=dble(slowness(iarr))*rtod8
          do ilay=1,numlay
            do irt=1,3
              xtop8=xtoplay(ilay)
              xbot8=xbotlay(ilay)
              nlev=numlev(ilay)
              isotropic=.true.
              if(ianisoflag.eq.1.and.ianiso(ilay).eq.1) then
                isotropic=.false.
              endif
c
c---- here we do the intergration for a particular p
c
              call wgray(pray8,ilay,xbot8,xtop8,nlev,isotropic,irt,
     #           delta8(ilay,irt),ttime8(ilay,irt),
     #           dddp8(ilay,irt),tstar8(ilay,irt),xturn8(ilay,irt),
     #           maxray,nray(ilay,irt),rayrad8(1,ilay,irt),
     #           raydelta8(1,ilay,irt),raytime8(1,ilay,irt))
c
              delta8(ilay,irt)=delta8(ilay,irt)*rtod8
              dddp8(ilay,irt)=dddp8(ilay,irt)*rtod8*rtod8
              do i=1,nray(ilay,irt)
                raydelta8(i,ilay,irt)=raydelta8(i,ilay,irt)*rtod8
              enddo
c
              if(iprtlv.gt.0) then
                write(6,"(i2,i3,7f10.3)")irt,ilay,xtop8,xbot8,
     #           delta8(ilay,irt),ttime8(ilay,irt),
     #           dddp8(ilay,irt),tstar8(ilay,irt),xturn8(ilay,irt)
              endif
c
c---- calculate the ray integrals for the upper and lower parts of source layer
c---- store the results in layers numlay+1 (upper) and numlay+2 (lower)
c
              if(ilay.eq.isl) then
                xtop8=xtoplay(ilay)
                xbot8=rad8
                if(xtop8-xbot8.lt.0.001d0) then
                  iupper=0
                else
                  iupper=1
                  klay=numlay+1
                  if(xtop8-xbot8.lt.1.d0) then
                    klev=min0(2,nlev)
                  else
                    klev=nlev
                  endif
                  call wgray(pray8,ilay,xbot8,xtop8,klev,isotropic,irt,
     #             delta8(klay,irt),ttime8(klay,irt),
     #             dddp8(klay,irt),tstar8(klay,irt),xturn8(klay,irt),
     #             maxray,nray(klay,irt),rayrad8(1,klay,irt),
     #             raydelta8(1,klay,irt),raytime8(1,klay,irt))
                  delta8(klay,irt)=delta8(klay,irt)*rtod8
                  dddp8(klay,irt)=dddp8(klay,irt)*rtod8*rtod8
                  do i=1,nray(klay,irt)
                    raydelta8(i,klay,irt)=raydelta8(i,klay,irt)*rtod8
                  enddo
                endif
                xtop8=rad8
                xbot8=xbotlay(ilay)
                if(xtop8-xbot8.lt.0.001d0) then
                  ilower=0
                else
                  ilower=1
                  klay=numlay+2
                  if(xtop8-xbot8.lt.1.d0) then
                    klev=min0(2,nlev)
                  else
                    klev=nlev
                  endif
                  call wgray(pray8,ilay,xbot8,xtop8,klev,isotropic,irt,
     #             delta8(klay,irt),ttime8(klay,irt),
     #             dddp8(klay,irt),tstar8(klay,irt),xturn8(klay,irt),
     #             maxray,nray(klay,irt),rayrad8(1,klay,irt),
     #             raydelta8(1,klay,irt),raytime8(1,klay,irt))
                  delta8(klay,irt)=delta8(klay,irt)*rtod8
                  dddp8(klay,irt)=dddp8(klay,irt)*rtod8*rtod8
                  do i=1,nray(klay,irt)
                    raydelta8(i,klay,irt)=raydelta8(i,klay,irt)*rtod8
                  enddo
                endif
              endif
            enddo
          enddo
c
c---- sum up the contributions from each layer
c
          iseq(1)=0
          if(iss.eq.1.and.ilower.eq.1) then
            iseq(1)=numlay+2
          else if(iss.eq.5.and.iupper.eq.1) then
            iseq(1)=numlay+1
          else if(iss.eq.2.and.ilower.eq.1) then
            iseq(1)=-(numlay+2)
          else if(iss.eq.6.and.iupper.eq.1) then
            iseq(1)=-(numlay+1)
          else
          endif
c
          dsum=0.
          tsum=0.
          dddpsum=0.
          tstarsum=0.
          nraypts(iarr)=0
          xxx=0.
          yyy=0.
          do iseg=1,nseg
            if(iseq(iseg).ne.0) then
              if(iseq(iseg).gt.0) then
                irt=2
              else
                irt=3-(ishflag*2)
              endif
              ilay=iabs(iseq(iseg))
              if(ilay.le.numlay) ilay=numlay+1-ilay
              dsum=dsum+sngl(delta8(ilay,irt))
              tsum=tsum+sngl(ttime8(ilay,irt))
              dddpsum=dddpsum+sngl(dddp8(ilay,irt))
              tstarsum=tstarsum+sngl(tstar8(ilay,irt))
              if(iprtlv.gt.0) then
              write(6,"('iseg,idir,irt,ilay,nray,delta',5i4,2f10.3)")
     #        iseg,idirseg(iseg),irt,ilay,nray(ilay,irt),dsum,tsum
              endif
              do ii=1,nray(ilay,irt)
                if(idirseg(iseg).eq.1) then
                  jj=ii
                else
                  jj=nray(ilay,irt)+1-ii
                endif
                nraypts(iarr)=nraypts(iarr)+1
                if(nraypts(iarr).gt.mxpts) then
                  stop 'ray points exceed mxpts -- cagcrays'
                endif
                rayrad(nraypts(iarr),iarr)=sngl(rayrad8(jj,ilay,irt))
                if(idirseg(iseg).eq.1) then
                  raydelta(nraypts(iarr),iarr)=xxx+
     #              sngl(raydelta8(jj,ilay,irt))
                  raytime(nraypts(iarr),iarr)=yyy+
     #              sngl(raytime8(jj,ilay,irt))
                else
                  raydelta(nraypts(iarr),iarr)=xxx+
     #            sngl(raydelta8(nray(ilay,irt),ilay,irt)-raydelta8(jj,ilay,irt))
                  raytime(nraypts(iarr),iarr)=yyy+
     #            sngl(raytime8(nray(ilay,irt),ilay,irt)-raytime8(jj,ilay,irt))
                endif
                if(ilay.le.numlay) then
                  iraylay(nraypts(iarr),iarr)=ilay
                else
                  iraylay(nraypts(iarr),iarr)=isl
                endif
                iraytp(nraypts(iarr),iarr)=irt*idirseg(iseg)
              enddo
              xxx=raydelta(nraypts(iarr),iarr)
              yyy=raytime(nraypts(iarr),iarr)
            endif
          enddo
          if(iprtlv.gt.0) then
          write(6,"('dsum,tsum,dddpsum,tstarsum:',4f10.4)")
     #          dsum,tsum,dddpsum,tstarsum
          endif
        enddo
c
c---- modify path if this is a diffracted wave
c
        if(narr.eq.1.and.idiffphase.ne.0.and.diffract.gt.0.0) then
          iarr=1
          ndiffract=1+int(diffract)
          diffinc=diffract/float(ndiffract)
          ifirstup=0
          do ipt=nraypts(iarr),1,-1
            if(iraytp(ipt,iarr).lt.0.and.iraytp(ipt-1,iarr).gt.0) then
              ifirstup=ipt
            endif
            if(ipt.ge.ifirstup) then
              rayrad(ipt+ndiffract,iarr)=rayrad(ipt,iarr)
              raydelta(ipt+ndiffract,iarr)=raydelta(ipt,iarr)+diffract
              iraytp(ipt+ndiffract,iarr)=iraytp(ipt,iarr)
              iraylay(ipt+ndiffract,iarr)=iraylay(ipt,iarr)
              raytime(ipt+ndiffract,iarr)=raytime(ipt,iarr)+
     #                 diffract*slowness(1)
            endif
          enddo
          if(iprtlv.gt.0) then
          write(6,"('before filling:',i4)") ndiffract
          endif
c
c---- define the upgoing leg to start with the diffracted part
c
         do i=1,ndiffract
           ipt=ifirstup-1+i
           rayrad(ipt,iarr)=rayrad(ifirstup-1,iarr)
           raydelta(ipt,iarr)=raydelta(ifirstup-1,iarr)+float(i-1)*diffinc
           iraytp(ipt,iarr)=iraytp(ifirstup+ndiffract,iarr)
           iraylay(ipt,iarr)=iraylay(ifirstup-1,iarr)
           raytime(ipt,iarr)=raytime(ifirstup-1,iarr)+
     #              float(i-1)*diffinc*slowness(1)
         enddo
          nraypts(iarr)=nraypts(iarr)+ndiffract
        endif
c
c---- calculate latitudes and longitudes for ray points
c
        do iarr=1,narr
          do ipt=1,nraypts(iarr)
            if(raydelta(ipt,iarr).gt.0.001) then
              call pdaz(epla,eplo,azep,raydelta(ipt,iarr),
     #             raylat(ipt,iarr),raylong(ipt,iarr))
            else
              raylat(ipt,iarr)=epla
              raylong(ipt,iarr)=eplo
            endif
          enddo
        enddo
c
c---- calculate travel time along straight line segments
c
        do iarr=1,narr
          ttt=0.
          ttte=0.
          do ipt=2,nraypts(iarr)
            x1=rayrad(ipt-1,iarr)*
     #         cosd(raylong(ipt-1,iarr))*sind(90.0-raylat(ipt-1,iarr))
            y1=rayrad(ipt-1,iarr)*
     #         sind(raylong(ipt-1,iarr))*sind(90.0-raylat(ipt-1,iarr))
            z1=rayrad(ipt-1,iarr)*cosd(90.0-raylat(ipt-1,iarr))
            r8=dble(rayrad(ipt-1,iarr))
            call evemell(r8,ell8,eta28,ierr)
            c1=(1.+sngl(ell8)*(0.3333333-cosd(90.0-raylat(ipt-1,iarr))**2))
c
            x2=rayrad(ipt,iarr)*
     #         cosd(raylong(ipt,iarr))*sind(90.0-raylat(ipt,iarr))
            y2=rayrad(ipt,iarr)*
     #         sind(raylong(ipt,iarr))*sind(90.0-raylat(ipt,iarr))
            z2=rayrad(ipt,iarr)*cosd(90.0-raylat(ipt,iarr))
            r8=dble(rayrad(ipt,iarr))
            call evemell(r8,ell8,eta28,ierr)
            c2=(1.+sngl(ell8)*(0.3333333-cosd(90.0-raylat(ipt,iarr))**2))
c
            xe1=x1*c1
            ye1=y1*c1
            ze1=z1*c1
            xe2=x2*c2
            ye2=y2*c2
            ze2=z2*c2
c
            dist=sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
            diste=sqrt((xe1-xe2)**2+(ye1-ye2)**2+(ze1-ze2)**2)
c
            if(dist.lt.0.001) then
            else
              if(iraytp(ipt-1,iarr).ne.iraytp(ipt,iarr)) then
                write(6,"('error -- not the same type of ray',i5)")ipt
              endif
              if(iraylay(ipt,iarr).ne.iraylay(ipt,iarr)) then
                write(6,"('error -- not the same layer')")
              endif
              r8=dble(rayrad(ipt-1,iarr))
              il=iraylay(ipt-1,iarr)
c
              isotropic=.true.
              if(ianisoflag.eq.1.and.ianiso(il).eq.1) then
                isotropic=.false.
              endif
              if(isotropic) then
                call evem(il,r8,.true.,rho8,vpv8,vph8,vsv8,vsh8,
     #             eta8,qmu8,qka8,ierr)
                if(iraytp(ipt-1,iarr).eq.2.or.iraytp(ipt-1,iarr).eq.-2) then
                  vgrp1=sngl(vpv8)
                else
                  vgrp1=sngl(vsv8)
                endif
                r8=dble(rayrad(ipt,iarr))
c --- the following line was missing - a bug
                il=iraylay(ipt,iarr)
                call evem(il,r8,.true.,rho8,vpv8,vph8,vsv8,vsh8,
     #             eta8,qmu8,qka8,ierr)
                if(iraytp(ipt-1,iarr).eq.2.or.iraytp(ipt-1,iarr).eq.-2) then
                  vgrp2=sngl(vpv8)
                else
                  vgrp2=sngl(vsv8)
                endif
              else
c
c---- calculate angles
c
                ddd=raydelta(ipt,iarr)-raydelta(ipt-1,iarr)
                sth1=sind(ddd)*rayrad(ipt-1,iarr)/dist
                sth2=sind(ddd)*rayrad(ipt,iarr)/dist
                if(abs(sth1).gt.1.) then
                  th1=90.0
                else
                  th1=asind(sth1)
                endif
                if(abs(sth2).gt.1.) then
                  th2=90.0
                else
                  th2=asind(sth2)
                endif
c                write(6,"('thetas:',7f12.6)")
c     #            sth1,sth2,th1,th2,ddd,rayrad(ipt-1,iarr),rayrad(ipt,iarr)
                iii=iabs(iraytp(ipt,iarr))
                r8=dble(rayrad(ipt-1,iarr))
                call evem(il,r8,.false.,rho8,vpv8,vph8,vsv8,vsh8,
     #             eta8,qmu8,qka8,ierr)
                call anivel(th1,iii,vpv8,vph8,vsv8,vsh8,eta8,vpha1,vgrp1,phi1)
                r8=dble(rayrad(ipt,iarr))
c --- the following line was missing - a bug
                il=iraylay(ipt,iarr)
                call evem(il,r8,.false.,rho8,vpv8,vph8,vsv8,vsh8,
     #             eta8,qmu8,qka8,ierr)
                call anivel(th2,iii,vpv8,vph8,vsv8,vsh8,eta8,vpha2,vgrp2,phi2)
                if(iprtlv.gt.0) then
                write(6,"(6f10.5)") r8,th2,vpha2,vgrp2,phi2
                endif
              endif
c
c---- calculate slowness
c
              slow1=1./vgrp1
              slow2=1./vgrp2
              ttt=ttt+dist*0.5*(slow1+slow2)
              ttte=ttte+diste*0.5*(slow1+slow2)
            endif
            if(iprtlv.gt.1) then
              write(6,"('ttt:',i4,3f10.4,i4)")
     #           ipt,ttt,ttte,raytime(ipt,iarr),iraylay(ipt,iarr)
            endif
          enddo
          if(iprtlv.gt.1) then
          write(6,"('Travel time along ray ',i2,':',2f10.3)") iarr,ttt,ttte
          endif
          ellcorr(iarr)=ttte-ttt
        enddo
      endif
c
c-----------------------------------------------------------------------
c
      if(itypcalc.gt.2) then
        do iarr=1,narr
          if(iarr.gt.maxarr) then
            stop 'too many arrivals - cagcrays'
          endif
          nrayptsx(iarr)=nraypts(iarr)
          if(nraypts(iarr).gt.maxraypts) then
            stop 'too many points for ray - cagcrays'
          endif
          do ipt=1,nraypts(iarr)
            iraytpx(ipt,iarr)=iraytp(ipt,iarr)
            iraylayx(ipt,iarr)=iraylayx(ipt,iarr)
            rayradx(ipt,iarr)=rayrad(ipt,iarr)
            raydeltax(ipt,iarr)=raydelta(ipt,iarr)
            raylatx(ipt,iarr)=raylat(ipt,iarr)
            raylongx(ipt,iarr)=raylong(ipt,iarr)
            raytimex(ipt,iarr)=raytime(ipt,iarr)
          enddo
        enddo
      endif
c
      if(ierr.eq.1) then
	write(6,"('ierr= ',i5,' in cagcrays')")ierr
	narr=0
      endif
      return
      end
c
      subroutine anivel(th,itp,vpv8,vph8,vsv8,vsh8,eta8,vpha,vgrp,phi)
      double precision vpv8,vph8,vsv8,vsh8,eta8
      if(itp.eq.1) then
        vpha=sqrt((cosd(th)*sngl(vsv8))**2+(sind(th)*sngl(vsh8))**2)
        phi=atand(0.5*sind(2.*th)*(sngl(vsh8)**2-sngl(vsv8)**2)/vpha**2)
        th2=th-phi
        vpha=sqrt((cosd(th2)*sngl(vsv8))**2+(sind(th2)*sngl(vsh8))**2)
        vgrp=vpha/cosd(phi)
      else if(itp.eq.2) then
        vpha=sngl(vpv8)
        vgrp=vpha
      else if(itp.eq.3) then
        vpha=sngl(vsv8)
        vgrp=vpha
      endif
      return
      end


 	  subroutine reset_ptab
c     initializes all arrays in rays.h to zeros.
c     is required when cagcrays is called multiple times with ireadagain
 	  include 'rays.h'
	   nlayers=0
	   do i=1,maxp
	     	pfull(i)=0.
	   		do j=1,maxl
	   			do k=1,3
	   				tfull(i,j,k)=0.
	   				dfull(i,j,k)=0.
	   				dddpfull(i,j,k)=0.
	   				tstarfull(i,j,k)=0.
	   				xturnfull(i,j,k)=0.
	   			enddo
	   		enddo
	   	enddo

		do j=1,maxl
			topfull(j)=0.
			botfull(j)=0.
			do k=1,3
				npfull(j,k)=0.
			enddo
		enddo

		do i=1,maxp
			do k=1,3
				tabove(i,k)=0.
				tbelow(i,k)=0.
				dabove(i,k)=0.
				dbelow(i,k)=0.
				dddpabove(i,k)=0.
				dddpbelow(i,k)=0.
				tstarabove(i,k)=0.
				tstarbelow(i,k)=0.
			enddo
		enddo

		do k=1,3
			npabove(k)=0.
			npbelow(k)=0.
		enddo

	   return
	   end









