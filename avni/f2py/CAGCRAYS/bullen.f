      subroutine getbullen(emodel,maxlay,capom,capg,grav,
     #                     vaisala,bullen,pressure)
      implicit double precision (a-h,o-z)
Cf2py intent(inout) emodel,maxlay,capom,capg
Cf2py intent(out) grav,vaisala,bullen,pressure
Cf2py depend(maxlay) grav,vaisala,bullen,pressure
      character*250 emodel
      include 'emcommonb.h'
      double precision xb(maxlay)
      double precision xt(maxlay)
      integer nlev(maxlay)
      integer iflu(maxlay)
      integer iani(maxlay)
      real*8 grav(maxlay)
      real*8 vaisala(maxlay)
      real*8 bullen(maxlay)
      real*8 pressure(maxlay)
      parameter (twopi=6.2831853072d0)
c     To make it generic for any planet, input argument instead
c      parameter (capom=7.292115d-5)
c      parameter (capg=6.6723d-11)

      if (maxlay.gt.maxlev) then
         write(6,"('Warning: Maximum layer maxlay cannot be greater than ',i5,' in emcommonb.h')") maxlev
        return
      endif
      open(1,file=emodel)
      call reademb(1,maxlay,nlay,xb,xt,nlev,iflu,iani,capom,capg,ierr)
      close(1)
      ilev=0
      do i=1,nlay
        do j=1,nlev(i)
	  ilev=ilev+1
	  thickness=xt(i)-xb(i)
	  step=thickness/float(nlev(i)-1)
	  drad=xb(i)+float(j-1)*step
          call evemb(i,drad,rho,gra,bigk,pres,vais,bull,ierr)
      grav(ilev) = gra
      vaisala(ilev) = vais
      bullen(ilev) = bull
      pressure(ilev) = pres
c          write(6,"(2i4,3f10.3,f9.1,3e13.4)")i,ilev,drad,rho,grav(i),
c     #          bigk*1.e-8,vaisala(i),bullen(i),pressure(i)
       enddo
      enddo
      end
