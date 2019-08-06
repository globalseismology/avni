      function asind(arg)
      parameter (twopi=6.28318530718)
      parameter (degtorad=twopi/360.)
      parameter (radtodeg=360./twopi)
      asind=radtodeg*asin(arg)
      return
      end
