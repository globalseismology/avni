      double precision function dsind(darg)
      real*8 darg
      real*8 dtwopi
      real*8 ddegtorad
      parameter (dtwopi=6.28318530718d0)
      parameter (ddegtorad=dtwopi/360.d0)
      dsind=dsin(darg*ddegtorad)
      return
      end
