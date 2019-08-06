      subroutine get_branch(imodel,phase,narr,delta,slowness,ddeltadp,
     #branch,ierror)
      character*16 phase
      parameter (maxarr=20)
      character*3 branch(maxarr)
      dimension diffslow_temp(maxarr)
      integer branchtemp(maxarr)
      character*1 ch1
      dimension ddeltadp(maxarr),temp(maxarr),temp2(maxarr)
      dimension slowness(maxarr),slownesstemp(maxarr)
      parameter(maxmodel=350)
      dimension deltalast(maxmodel)
      character*16 phaselast(maxmodel)
      character*3 branchlast(maxmodel,maxarr)
      dimension slownesslast(maxmodel,maxarr)
      integer inumbranchlast(maxmodel),ibranchoffset(maxmodel)

c     deltalast is used to check that we are going from large to small delta. P, PP and SS phases require this for indexing      
      data deltalast /350*99999999.0/
      data phaselast /350*' '/
      data slownesslast /7000*0.0/
      data branchlast /7000*' '/
      data inumbranchlast /350*0/
      data ibranchoffset /350*0/
      data imodellast /0/
      save

	  if(imodel.gt.maxmodel) stop'imodel.gt.maxmodel'
	  if(imodel.le.0) stop'imodel.le.0'
c     This subroutne returns a character array specifying the branch
c     ab, bc, or numbers 1,2,3, etc. The branches are sorted accroeding to
c     to ddeltadp for core phases and by looking at slowness vairations with delta
c     for main phases with lots of upper mantle triplications like
c     S, SS, PP and SS. This causes the following results in the phases:
c     (1) P, PP and S: reurn to prograde after transition zone triplications
c     (2) PKP: ddeltadp positive for ab; ddeltadp negative for bc. bc indicates the prograde branch of the PKP caustic
c     (3) PKKP: decreasing order ddeltadp - ab followed by bc.
c     (4) P'P': decreasing order. use branch ab typically as in AK135
c     (5) SKP: use bc like in AK135  

	   ierror=0
	   do ii=1,narr
	   	branch(ii)=' '
	   enddo

	  	
c     fine phase and replace the name
	  if(phase(1:lnblnk(phase)).eq.'P'.or.phase(1:lnblnk(phase)).eq.'S'.or.
     #   phase(1:lnblnk(phase)).eq.'PP'.or.phase(1:lnblnk(phase)).eq.'SS'.or.
     #   phase(1:lnblnk(phase)).eq.'SP'.or.phase(1:lnblnk(phase)).eq.'PPP') then
		  if(delta.ge.25.0.and.delta.le.335.0) then
		  else
			write(6,"('Error: finding the ',a,' branch has not been tested at delta<25')") phase(1:lnblnk(phase))
			stop
		  endif

	   	  
c	   	  THE FOOLOWING IS NEEDED AS WE TRACK BRACH BASED ON A DECREASING DELTA
		  if(phase.eq.phaselast(imodel)) then  !if phase the smae as earlier call, then delta should be in decreasing order
			if(delta.gt.deltalast(imodel)) then
				write(6,"('Error: The current delta ',f7.2,' is greater than previous delta ',f7.2)") delta,deltalast(imodel)
				write(6,"('The delta should move from large to small values for indexing of phases like P,PP,S and SS')") 
				stop
			endif
		  else
			phaselast(imodel)=phase
			inumbranchlast(imodel)=0
			do ii=1,maxarr
				slownesslast(imodel,ii)=0. !intialize the last slowness of all branches for a new phase
				branchlast(imodel,ii)=' '
			enddo
			deltalast(imodel)=99999999.0
		  endif  

		  
		  if(inumbranchlast(imodel).eq.0) then !if this is the first time we are doing this for a phase
		  	if(narr.gt.1) then
		  		ierror=1
		  		write(6,"('Error: narr.gt.1 in the first delta in imodel = ',i3,'. maybe increase deltamax')")imodel
		  		return
		  	endif
		  	inumbranchlast(imodel)=1
		  	branch(1)='1'
		  	branchlast(imodel,1)=branch(1)
		  	slownesslast(imodel,1)=slowness(1)
		  else ! we already have one or more branches 
			if(narr.ge.inumbranchlast(imodel)) then  !if there are more branches than earlier 
				ilastmaxbranch=0
				do j=1,inumbranchlast(imodel)  !loop over previous branches and find the best current rays
			  	  	diff_slowness=999999.
					do i=1,narr   ! choose preferred branch with smallest difference with slownesslast				  
						if(abs(slownesslast(imodel,j)-slowness(i)).le.diff_slowness) then
							diff_slowness=abs(slownesslast(imodel,j)-slowness(i))
						    ipreflastbranch=i
						endif
					enddo

					read(branchlast(imodel,j),'(i1.1)') itemp 
					if(itemp.gt.ilastmaxbranch) ilastmaxbranch=itemp
					branch(ipreflastbranch)=branchlast(imodel,j)
C					print*,j,ipreflastbranch,' ',branch(ipreflastbranch),slownesslast(imodel,j)
			  	enddo
			  	!find the remaining branches
			  	inumbranchnew=ilastmaxbranch
			  	do i=1,narr
					if(branch(i).eq.' ') then
						inumbranchnew=inumbranchnew+1
						write(ch1, '(i1.1)') inumbranchnew
						branch(i)=ch1
C						print*,'NEW --- ',inumbranchnew,' ',branch(i),slowness(i)
					endif
					
				enddo				
				
				if(inumbranchnew.gt.9) stop'Error:get_branch cannot deal with more than 9 branches'

				! store the branches
			  	do i=1,narr
					branchlast(imodel,i)=branch(i)
			  		slownesslast(imodel,i)=slowness(i)
			  	enddo
			  	inumbranchlast(imodel)=narr
		  else !if there are less branches than earlier
			do i=1,narr   ! choose preferred branch with smallest difference with slownesslast	
				diff_slowness=999999.
				ilastmaxbranch=0
				do j=1,inumbranchlast(imodel)  !loop over previous branches and find the best current rays
					if(abs(slownesslast(imodel,j)-slowness(i)).le.diff_slowness) then
						diff_slowness=abs(slownesslast(imodel,j)-slowness(i))
						ipreflastbranch=j
						diffslow_temp(i)=diff_slowness
					endif
					read(branchlast(imodel,j),'(i1.1)') itemp 
					if(itemp.gt.ilastmaxbranch) ilastmaxbranch=itemp
				enddo
				branch(i)=branchlast(imodel,ipreflastbranch)
			enddo
			
			do i=2,narr
					do j=1,i-1
						if(branch(i).eq.branch(j)) then
							write(ch1, '(i1.1)') ilastmaxbranch+1
C							print*,branch(i),branch(j),diffslow_temp(i),diffslow_temp(j)
C							pause
							if(diffslow_temp(i).le.diffslow_temp(j)) then
								branch(j)=ch1
							else
								branch(i)=ch1
							endif
						endif	
					enddo
			enddo	
C			do i=1,narr 	
C				print*,i,ipreflastbranch,' ',branch(i),slownesslast(imodel,ipreflastbranch)
C		  	enddo
C		  	
			! store the branches
		  	do i=1,narr
				branchlast(imodel,i)=branch(i)
			  	slownesslast(imodel,i)=slowness(i)
		  	enddo
		  	inumbranchlast(imodel)=narr
		  	ibranchoffset(imodel)=inumbranchlast(imodel)-1
		  endif	
		 endif 
	  else if(phase(1:lnblnk(phase)).eq.'PKKP'.or.phase(1:lnblnk(phase)).eq.'SKP') then
		if(narr.eq.1) then
			if(ddeltadp(1).ge.0.) stop'Error: A single arrival PKKP or SKP at positive ddeltadp. ab does not exist '
		else if(narr.gt.2) then
			stop'Error: Cannot have more than 2 branches for PKKP or SKP. Inner core phases are named differently. '
		endif	
		do i=1,narr
			if(ddeltadp(i).ge.0.) then
				branch(i)='ab'
			else 
				branch(i)='bc'
			endif	
		enddo
	  else if(phase(1:lnblnk(phase)).eq.'PKP'.or.phase(1:lnblnk(phase)).eq.'PKPPKP') then
		if(narr.eq.1) then
	  	  	if(ddeltadp(1).le.0.) stop'Error: A single arrival PKP or PKPPKP at negative ddeltadp'
		else if(narr.gt.2) then
			stop'Error: Cannot have more than 2 branches for PKP or PKPPKP. Inner core phases are named differently. '
		endif	
		do i=1,narr
			if(ddeltadp(i).ge.0.) then
				branch(i)='ab'
			else 
				branch(i)='bc'
			endif	
		enddo
	  else
	  	if(narr.eq.1) then ! only one arrival for undefined phased so simply give a branch name 1
	  		branch(1)='1'
	  	else
		  	ierror=2
	  		write(6,"('Error: The phase ',a,' with multiple arrivals not specified in get_branch')") phase(1:lnblnk(phase))
	  	endif
	  endif
		
	  deltalast(imodel)=delta
	  imodellast=imodel
      return
      end
