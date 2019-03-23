      subroutine tadder(jy,jd,ih,im,fs,jy1,jd1,ih1,im1,fs1,secs)
c---- subroutine to add secs seconds to a given time
c
      double precision secs,ss
      if(secs.gt.9.0d11) then
        jy1=0
        jd1=0
        ih1=0
        im1=0
        fs1=0.
        return
      endif
      ss=dble(fs)+secs
      mins=idint(ss/60.d0)
      if(ss.lt.0.0d0) mins=mins-1
      fs1=sngl(ss-60.d0*dfloat(mins))
      if(fs1.lt.0.) fs1=0.
      mm=im+mins
      ihrs=mm/60
      if(mm.lt.0.and.mod(mm,60).ne.0) ihrs=ihrs-1
      im1=mm-60*ihrs
      ihh=ih+ihrs
      idays=ihh/24
      if(ihh.lt.0.and.mod(ihh,24).ne.0) idays=idays-1
      ih1=ihh-24*idays
      jd1=jd+idays
      jy1=jy
   10 continue 
      leap=lpyr(jy1)
      ld=leap+365
      if(jd1.gt.ld) then
        jd1=jd1-ld
        jy1=jy1+1
        goto 10
      else if(jd1.le.0) then
        jy1=jy1-1
        leap=lpyr(jy1)
        ld=leap+365
        jd1=jd1+ld
        go to 10
      endif
      return
      end

