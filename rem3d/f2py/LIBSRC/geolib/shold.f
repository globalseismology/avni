c --- evaluate old basis functions: Miaki's spherical harmonics.
c --- at this point, we are using real sph har, so no conjugate is necesary

      subroutine shold(rlat,rlon,lold,numold,yreal)
Cf2py intent(inout) rlat,rlon,lold,numold
Cf2py intent(out) yreal
Cf2py depend(numold) yreal

      implicit none

      integer lold,numold,l,m,i
      real yreal(numold)
      real theta,phi,rlat,rlon

      integer lold2
      parameter(lold2=200)
      real Xlm(0:lold2, 0:lold2)
      real XP
      real Pi
      parameter (Pi = 3.14159265358979)
c
c ---
c
      if(lold2.lt.lold) stop 'lold2 must not be smaller than lold'
c
c --- convert degrees into radians (why rlat-90 ?????)
c
      if (rlat.gt.0.0) then
      theta = (rlat-90)*Pi/180.0
      elseif (rlat.lt.0.0) then
      theta = (90-rlat)*Pi/180.0
      elseif (rlat.eq.0) then
      theta = Pi/2.0
      end if
      phi = rlon*Pi/180.0
c      write(101,*) theta

c
c --- calculate associated Legendre functions for this theta
c
      do l = 0, lold
        do m = 0, l
          Xlm(l,m) = XP(l,m,theta)
        end do  ! end m loop
      end do  ! end l loop
c
c --- calculate real spherical harmonics
c
      i = 1
      do l = 0, lold
        do m = 0, l
        if (i.gt.numold) stop 'i.gt.numold'
        if (m.eq.0) then  ! zonal term
          yreal(i) = Xlm(l,m)
            i = i + 1
          elseif (m.ne.0) then ! non-zonal terms
            yreal(i) = ((-1.0)**m)*Sqrt(2.0)*
     #                  Xlm(l,m)*Cos(m*phi) ! cosine part
            yreal(i+1) = ((-1.0)**m)*Sqrt(2.0)*
     #                    Xlm(l,m)*Sin(m*phi) ! sine part
          i = i + 2
          end if  ! if (m.ne.0)...
        end do  ! end m loop
      end do  ! end l loop
      if (i-1.ne.numold) stop 'i-1.ne.numold'
      end
c
c=============================
c
      subroutine sub_indices(lmax,numcoe,indexl,indexm,indextyp)
      implicit none
      integer lmax,l,m,i,numcoe
      integer indexl(*),indexm(*),indextyp(*)
c
      i = 1
      do l = 0, lmax
        do m = 0, l
        if (i.gt.numcoe) stop 'i.gt.numcoe'
        if (m.eq.0) then  ! zonal term
          indexl(i)=l
          indexm(i)=m
          indextyp(i)=0
            i = i + 1
          elseif (m.ne.0) then ! non-zonal terms
          indexl(i)=l
          indexm(i)=m
          indextyp(i)=1
          indexl(i+1)=l
          indexm(i+1)=m
          indextyp(i+1)=2
          i = i + 2
          end if  ! if (m.ne.0)...
        end do  ! end m loop
      end do  ! end l loop
      if (i-1.ne.numcoe) then
      print*,i-1,numcoe
      stop 'i-1.ne.numcoe'
      endif
      end

cc*****************************************************
      real function XP(l,m,theta)

c  This function calculates the colatitude dependent part of the
c  spherical harmonics.

      integer l, m
c  l = spherical harmonics angular degree
c  m = spherical harmonics angular order

      real theta
c  colatitude

      real Pi
      parameter (Pi = 3.14159265358979)

      real factrl, plgndr

      XP = ((-1.0)**m)*Sqrt((2.0*float(l)+1.0)/(4.0*Pi))*
     +         Sqrt(FACTRL(l-m)/FACTRL(l+m))*
     +         PLGNDR(l,m,Cos(theta))

      return
      end


c*******************************************
      real FUNCTION FACTRL(N)
c  from the Numerical Recipes...

      integer N

      real A(33)
      integer NTOP
      DATA NTOP,A(1)/0,1./

      integer j
      real gammln

      IF (N.LT.0) THEN
        STOP 'negative factorial'
      ELSE IF (N.LE.NTOP) THEN
        FACTRL=A(N+1)
      ELSE IF (N.LE.32) THEN
        DO 11 J=NTOP+1,N
          A(J+1)=J*A(J)
11      CONTINUE
        NTOP=N
        FACTRL=A(N+1)
      ELSE
        FACTRL=EXP(GAMMLN(N+1.))
      ENDIF
      RETURN
      END



c*******************************************
      real FUNCTION GAMMLN(XX)
c  from the Numerical Recipes...

      REAL*8 COF(6),STP,HALF,ONE,FPF,X,TMP,SER
      DATA COF,STP/76.18009173D0,-86.50532033D0,24.01409822D0,
     *    -1.231739516D0,.120858003D-2,-.536382D-5,2.50662827465D0/
      DATA HALF,ONE,FPF/0.5D0,1.0D0,5.5D0/

      real xx

      integer j

      X=Dble(XX)-ONE
      TMP=X+FPF
      TMP=(X+HALF)*LOG(TMP)-TMP
      SER=ONE

      DO 11 J=1,6
        X=X+ONE
        SER=SER+COF(J)/X
11    CONTINUE

      GAMMLN=Sngl(TMP+LOG(STP*SER))

      RETURN
      END



c*******************************************
      real FUNCTION PLGNDR(L,M,X)
c  This is from Numerical Recipes...

      integer l, m
      real x

      real pmm, somx2, fact, pmmp1, pll

      integer i, ll

      IF(M.LT.0.OR.M.GT.L.OR.ABS(X).GT.1.) stop 'bad arguments'
      PMM=1.
      IF(M.GT.0) THEN
        SOMX2=SQRT((1.-X)*(1.+X))
        FACT=1.
        DO 11 I=1,M
          PMM=-PMM*FACT*SOMX2
          FACT=FACT+2.
11      CONTINUE
      ENDIF
      IF(L.EQ.M) THEN
        PLGNDR=PMM
      ELSE
        PMMP1=X*(2*M+1)*PMM
        IF(L.EQ.M+1) THEN
          PLGNDR=PMMP1
        ELSE
          DO 12 LL=M+2,L
            PLL=(X*(2*LL-1)*PMMP1-(LL+M-1)*PMM)/(LL-M)
            PMM=PMMP1
            PMMP1=PLL
12        CONTINUE
          PLGNDR=PLL
        ENDIF
      ENDIF
      RETURN
      END

