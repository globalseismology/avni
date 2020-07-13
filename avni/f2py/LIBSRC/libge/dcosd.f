      double precision function dcosd(darg)
      real*8 darg
      real*8 dtwopi
      real*8 ddegtorad
      parameter (dtwopi=6.28318530718d0)
      parameter (ddegtorad=dtwopi/360.d0)
      dcosd=dcos(darg*ddegtorad)
      return
      end
