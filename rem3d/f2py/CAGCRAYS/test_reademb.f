      implicit double precision (a-h,o-z)
      character*80 emodel
      parameter (maxlay=1000)
      double precision xb(maxlay)
      double precision xt(maxlay)
      dimension nlev(maxlay)
      dimension iflu(maxlay)
      dimension iani(maxlay)
      parameter (twopi=6.2831853072d0)
      logical isotropic
      common /plevel/iprtlv
c
      iprtlv=2
      call getarg(1,emodel)
      write(6,"('type name of model to read')")
      call getarg(1,emodel)
      if(lnblnk(emodel).le.0) then
        stop 'give the name of the earthmodel as an argument'
      endif
c      read(5,"(a)") emodel 
      write(6,"(a)") emodel
      open(1,file=emodel)
      call reademb(1,maxlay,nlay,xb,xt,nlev,iflu,iani,ierr)
      close(1)
      ilev=0
      do i=1,nlay
        do j=1,nlev(i)
	  ilev=ilev+1
	  thickness=xt(i)-xb(i)
	  step=thickness/float(nlev(i)-1)
	  drad=xb(i)+float(j-1)*step
          call evemb(i,drad,rho,grav,bigk,pressure,vaisala,bullen,ierr)
          write(6,"(2i4,3f10.3,f9.1,3e13.4)")i,ilev,drad,rho,grav,
     #          bigk*1.e-8,vaisala,bullen,pressure
       enddo
      enddo
      end
