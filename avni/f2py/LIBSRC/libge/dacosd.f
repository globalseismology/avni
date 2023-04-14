      double precision function dacosd(darg)
      real*8 darg
      real*8 dtwopi
      real*8 dradtodeg
      real*8 ddegtorad
      parameter (dtwopi=6.28318530718d0)
      parameter (ddegtorad=dtwopi/360.d0)
      parameter (dradtodeg=360.d0/dtwopi)
      dacosd=dradtodeg*dacos(darg)
      return
      end
