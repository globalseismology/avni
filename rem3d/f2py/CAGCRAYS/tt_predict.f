      subroutine tt_predict(delstep,filename,iallpsv,iellip,
     #epla,eplo,depth,phase_file)
      character*80 filename,string
      character*16 phase
      real*8 delstep,epla,eplo,depth
      integer iallpsv,iellip
c
      parameter (maxarr=20)
      dimension ttime(maxarr)
      dimension slowness(maxarr)
      dimension dtddepth(maxarr)
      dimension ddeltadp(maxarr)
      dimension tstar(maxarr)
      integer imodel(maxarr)
      dimension ellcorr(maxarr)
      character*3 branch(maxarr)
      integer iphasebranchopen(maxarr)
c
      parameter (maxraypts=5000)
      dimension nraypts(maxarr)
      dimension rayrad(maxraypts,maxarr)
      dimension raydelta(maxraypts,maxarr)
      dimension raylat(maxraypts,maxarr)
      dimension raylong(maxraypts,maxarr)
      dimension raytime(maxraypts,maxarr)
      dimension iraytp(maxraypts,maxarr)
      dimension iraylay(maxraypts,maxarr)
c
      character*80 systemstring

      parameter(nstations=5000)
      character*80 station_file,phase_file
      character*12 a
      character*3 d
      character*3 branchlist(nstations)
      character*12 station(nstations),phaselist(nstations),station_name
      dimension slat(nstations)
      dimension slon(nstations)
      dimension deltarr(nstations)
      dimension narrivals(nstations)
      character*80 pathfilename
      dimension distrange(nstations,2)
      character*1 str1
      character*2 str2
      character*3 str3
      character*4 str4
      logical exists
      character*80 getunx

c
      common/plevel/iprtlv
c-------------------------------------------------------------------------------
c---- this gives the branches in 1,2,3 for main phases (1 stands for large delta
c     when we have ne arrival for S  and P. more branches are added at smaller delta
c
c
c      call chekcl('|     :r:1:RAYS file'//
c     #            '|   -s:o:1:Assume all phases are PSV [0]'//
c     #            '|   -t:o:1:step in delta [1.]'//
c     #            '|   -e:o:1:print out ellipticity corrections [0]'//
c     #            '|   -v:o:1:Verbosity level [0]|')
      lu=1
c      string=getunx('-t',1,nbyts)
c      read(string,*) delstep
      write(6,"('step in delta, in degrees')")
c      read(5,*) delstep

      write(6,"('type name of the .RAYS file: ')")
c      read(5,"(a)") filename
c      filename=getunx('',1,lfile)

c      string=getunx('-s',1,nbyts)
c      read(string,*) iallpsv
      write(6,"('Assume all phases are PSV')")
c      read(5,*) iallpsv

c     string=getunx('-v',1,nbyts)
c      read(string,*) iprtlv
      iprtlv= 0
      write(6,"('type if ellipticity is needed')")
c      read(5,*) iellip
c      string=getunx('-e',1,nbyts)
c      read(string,*) iellip
c	  if only calculating tt and t* use itypcalc=1, otherwise use itypcalc=3
      	itypcalc=3
c
c------------
      write(6,"('type epla,eplo,dep')")
c      read(5,*) epla,eplo,depth

c
c-----read in phase file
c
      write(6,"('type phasefile')")
c      read(5,"(a)") phase_file
c      phase_file='phasefile'
      print *, 'READING PHASE FILE'
      open(1,file=phase_file)
      ios=0
      np=0
      do while (ios.eq.0)
	  	read(1,"(a)",iostat=ios) string
	 	 if(ios.ne.0) then
		 else
	    	lstr=lnblnk(string)
	    	if(string(1:1).eq.'#') then
	    	else
	    		if(string(13:15).eq.' | '.and.string(20:22).eq.' | '.and.string(27:29).eq.' | '.and.string(29:30).ne.'  ') then
	    		else
	    			stop'string(13:15)==" | ".and.string(20:22)==" | ".and.string(27:29)==" | ".and.string(30).ne." "'
	    		endif
		  		read(string(1:lstr),"(a12,' | ',f4.1,' | ',f4.1,' | ',a3)") a,b,c,d
				!this is done to store multiple entries of the same phase into different models in get_branch
				iphasebefore=0
				np=np+1
				do ii=1,np-1
					if(phaselist(ii)(1:lnblnk(phaselist(ii))).eq.a(1:lnblnk(a))) then
						if(branchlist(ii)(1:lnblnk(branchlist(ii))).eq.d(1:lnblnk(d))) then
							write(6,"('two entries in phaselist have duplicate phase and branch - ',2a)")a(1:lnblnk(a)),d(1:lnblnk(d))
							stop
						endif
						iphasebefore=iphasebefore+1
					endif
				enddo
				imodel(np)=iphasebefore+1
				phaselist(np)=a
				distrange(np,1)=b
				distrange(np,2)=c
				if(b>c) stop'provide distance range in ascending order in phasefile'
				branchlist(np)=d
			endif
		endif
      enddo
      close(1)
      print *, 'DONE READING PHASE FILE'
      if(np.eq.0) stop'np.eq.0'
	  print*,'phaselist     distrange(1)-distrange(2)    branchlist    imodel'
      do ii=1,np
		print*,phaselist(ii),distrange(ii,1),'-',distrange(ii,2),' ',branchlist(ii),imodel(ii)
      enddo
c
	  do iphase=1,np
	 	phase=phaselist(iphase)
		lphase=lnblnk(phase)
		if(iallpsv.eq.1) then
			ishflag=0 ! all phases are P-SV
		else ! find out
			call get_ishflag(phase,ishflag) !see whether SH or P-SV
		endif
		ianisoflag=1   ! use anisotropic tracing,
				               ! checked that it gives same answer as ianisoflag=0 for iso models

		inquire(file='TT_OUTPUT',exist=exists)
		if(.not.exists) then
			systemstring='mkdir TT_OUTPUT'
        	ierror=system(systemstring)
        endif

c       make a list fr this phase
	    slat(1)=epla
	    slon(1)=eplo+distrange(iphase,1)
	    deltarr(1)=distrange(iphase,1)
	    icount=1
	    write (str1,"(i1)") icount
	 	station(1)='SYNTH_000'//str1
	    do while (slon(icount).lt.distrange(iphase,2))
	  		icount=icount+1
	  		slat(icount)=epla
	  		slon(icount)=slon(icount-1)+delstep
	    	deltarr(icount)=deltarr(icount-1)+delstep
			if(icount.gt.999) then
	  			write (str4,"(i4)") icount
	  			station(icount)='SYNTH_'//str4
	  		else if (icount.gt.99.and.icount.le.999) then
	  			write (str3,"(i3)") icount
	  			station(icount)='SYNTH_0'//str3
	  		else if (icount.gt.9.and.icount.le.99) then
	  			write (str2,"(i2)") icount
	  			station(icount)='SYNTH_00'//str2
	  		else if (icount.gt.0.and.icount.le.9) then
	  			write (str1,"(i1)") icount
	  			station(icount)='SYNTH_000'//str1
	  		endif
	    enddo

		do iarr=1,maxarr
			iphasebranchopen(iarr)=0
		enddo
		iphaseopen=0
      	! loop over all stations
      	do numstation=icount,1,-1
      		stla=slat(numstation)
      		stlo=slon(numstation)
      		delta=deltarr(numstation)
      		if(delta.ge.distrange(iphase,1).and.delta.le.distrange(iphase,2)) then
c
c               convert to geocentric coordinates as those are needed by cagcrays-doesnt matter for long
	  			call conv2geocen(epla,eplageoc)
          		call conv2geocen(stla,stlageoc)

	  			ireadagain=0  !This is 1 if we read the RAYS file everytime cagcrays_pm is called
        		call cagcrays_pm(lu,filename,eplageoc,eplo,depth,stlageoc,stlo,
     #     phase,ishflag,ianisoflag,itypcalc,delta,azep,azst,
     #     maxarr,narr,ttime,slowness,dtddepth,ddeltadp,tstar,ellcorr,
     #     maxraypts,nraypts,iraytp,iraylay,
     #     rayrad,raydelta,raylat,raylong,raytime,ireadagain)
c

	      		if(itypcalc.ge.1) then
    	    		if(narr.eq.0.or.narr.gt.1) write(6,"(a16,2x,i3,' arrivals at',f9.3,' degrees')") phase,narr,delta
      				if(narr.gt.0) then !if number of rays traced is greater than 0
c              		  find the branch of the phase
      				  call get_branch(imodel(iphase),phase,narr,delta,slowness,ddeltadp,branch,ierror)
        			endif
C        			do ii=1,narr
C        				print*,phase,delta,ddeltadp(ii),slowness(ii),' ',branch(ii)
C        			enddo
      			endif
c
c
				if(itypcalc.ge.2) then
					if(narr.gt.0) then !if number of rays traced is greater than 0
							do iarr=1,narr
							    !first check if this branch has been excluded or not
       				    		if(branchlist(iphase)(1:lnblnk(branchlist(iphase))).eq.'all'.or.
     #                       	branchlist(iphase)(1:lnblnk(branchlist(iphase))).eq.'ALL') then
       				    			includebranch=1
       				    			inquire(file='TT_OUTPUT/'//phase(1:lphase),exist=exists)
									if(.not.exists) then
										systemstring='mkdir TT_OUTPUT/'//phase(1:lphase)
										ierror=system(systemstring)
									endif
       				    			if(iarr.eq.1) then
									    station_name=station(numstation)
									    pathfilename='TT_OUTPUT/'//phase(1:lphase)//'/path_'//station_name(1:lnblnk(station_name))
       				    				open(8,file=pathfilename)
										write(8,"('#PATH:',5f9.3)") eplageoc,eplo,depth,stlageoc,stlo
										lfilename=lnblnk(filename)
										write(8,"('#MODEL:',a)") filename(1:lfilename)
										write(8,"('#IFANI:',i2)") ianisoflag
										write(8,"('#IFSH:',i2)") ishflag
										write(8,"('#PHASE:',a)") phase(1:lphase)
										write(6,"('..... written ',a)") pathfilename(1:lnblnk(pathfilename))
									endif

									string='TT_OUTPUT/'//phase(1:lphase)//'/ttandtstar.txt'
									if(iphaseopen.ne.1) then
										open(11,file=string)
										if (iellip.gt.0) then
											write(11,"('delta,p,T,dDdp,tstar,ellcorr,branch:')")
										else
											write(11,"('delta,p,T,dDdp,tstar,branch:')")
										endif
										iphaseopen=1
									endif

       				    		else if(branchlist(iphase)(1:lnblnk(branchlist(iphase))).eq.
     #                       	branch(iarr)(1:lnblnk(branch(iarr)))) then
       				    			includebranch=2
									inquire(file='TT_OUTPUT/'//phase(1:lphase)//branch(iarr)(1:lnblnk(branch(iarr))),exist=exists)
									if(.not.exists) then
										systemstring='mkdir TT_OUTPUT/'//phase(1:lphase)//branch(iarr)(1:lnblnk(branch(iarr)))
										ierror=system(systemstring)
									endif
									string='TT_OUTPUT/'//phase(1:lphase)//branch(iarr)(1:lnblnk(branch(iarr)))//'/ttandtstar.txt'
       				    			inquire(file=string(1:lnblnk(string)),exist=exists)
									if(iphasebranchopen(imodel(iphase)).eq.0) then
										open(11+imodel(iphase),file=string)
										if (iellip.gt.0) then
											write(11+imodel(iphase),"('delta,p,T,dDdp,tstar,ellcorr,branch:')")
										else
											write(11+imodel(iphase),"('delta,p,T,dDdp,tstar,branch:')")
										endif
										iphasebranchopen(imodel(iphase))=1
									endif

									station_name=station(numstation)
									pathfilename='TT_OUTPUT/'//phase(1:lphase)//branch(iarr)(1:lnblnk(branch(iarr)))//
     #                       	        '/path_'//station_name(1:lnblnk(station_name))
       				    			open(8,file=pathfilename)
									write(8,"('#PATH: ',5f9.3)") eplageoc,eplo,depth,stlageoc,stlo
									lfilename=lnblnk(filename)
									write(8,"('#MODEL: ',a)") filename(1:lfilename)
									write(8,"('#IFANI: ',i2)") ianisoflag
									write(8,"('#IFSH: ',i2)") ishflag
									write(8,"('#PHASE: ',a)") phase(1:lphase)
									write(6,"('..... written ',a)") pathfilename(1:lnblnk(pathfilename))
       				    		else
       				    			includebranch=0
       				    		endif

       				    		!now include the branches as required
       				    		if(includebranch.gt.0) then
									write(8,"('#BRANCH: ',a)") branch(iarr)(1:lnblnk(branch(iarr)))
									write(8,"('#TTPARA: ',5f12.5)") delta,slowness(iarr),ttime(iarr),
     #                    			ddeltadp(iarr),tstar(iarr)
									write(8,"('#START#')")
									do i=1,nraypts(iarr)
										write(8,"(i5,i3,i4,5f10.3)")
     #      							i,iraytp(i,iarr),iraylay(i,iarr),rayrad(i,iarr),
     #      							raydelta(i,iarr),raylat(i,iarr),raylong(i,iarr),raytime(i,iarr)
									enddo
									write(8,"('#END#')")
									if(includebranch.eq.2) then
										close(8)
									endif

									if(includebranch.eq.1) iout=11
									if(includebranch.eq.2) iout=11+imodel(iphase)
									if (iellip.gt.0) then
         								write(iout,"(f10.4,f9.5,f9.3,f13.3,f7.3,f7.3,a3)")
C         								write(6,"(f10.4,f9.5,f9.3,f13.3,f7.3,f7.3,a3)")
     #          						delta,slowness(iarr),
     #          						ttime(iarr),ddeltadp(iarr),tstar(iarr),ellcorr(iarr),branch(iarr)(1:lnblnk(branch(iarr)))
									else
						   		 		write(iout,"(f10.4,f9.5,f9.3,f13.3,f7.3,a3)")
C						   		 		write(6,"(f10.4,f9.5,f9.3,f13.3,f7.3,a3)")
     #          						delta,slowness(iarr),
     #          						ttime(iarr),ddeltadp(iarr),tstar(iarr),branch(iarr)(1:lnblnk(branch(iarr)))
     								endif
								endif
							enddo
							if(includebranch.eq.1) then
								close(8)
							endif
					endif
				endif
      		endif
	  	enddo
	  	if(includebranch.eq.1) close(11)
		if(includebranch.eq.2) then
			if(iphasebranchopen(imodel(iphase)).ne.0) then
				close(11+imodel(iphase))
			endif
		endif

	  enddo
      end subroutine
