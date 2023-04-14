      subroutine monday(iy,jd,mo,id)
      integer*4 mon(12)
      data mon/0,31,59,90,120,151,181,212,243,273,304,334/
      ld=lpyr(iy)
      mo=1
      do 10 i=2,12
        imon=mon(i)
        if(i.gt.2) imon=imon+ld
        if(jd.le.imon) go to 11
        mo=mo+1
   10 continue
   11 continue
      id=jd-mon(mo)
      if(mo.gt.2) id=id-ld
      return
      end
