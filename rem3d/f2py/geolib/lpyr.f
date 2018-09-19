      integer function lpyr(iy)
c
c---- returns 1 if iy is a leap year
c
      lpyr=0
      if(mod(iy,400).eq.0) then
        lpyr=1
      else if(mod(iy,4).eq.0) then
        lpyr=1
        if(mod(iy,100).eq.0) then
          lpyr=0
        endif
      endif
      return     
      end
