      subroutine delazgc(eplat,eplong,stlat,stlong,delta,azep,azst)
c
Cf2py intent(inout) eplat,eplong,stlat,stlong
Cf2py intent(out) delta,azep,azst
      real*8 deplat,deplong,dstlat,dstlong,ddelta,dazep,dazst
      deplat=dble(eplat)
      deplong=dble(eplong)
      dstlat=dble(stlat)
      dstlong=dble(stlong)
c      print*,deplat,deplong,dstlat,dstlong
c
      call ddelazgc(deplat,deplong,dstlat,dstlong,ddelta,dazep,dazst)
c
      delta=sngl(ddelta)
      azep=sngl(dazep)
      azst=sngl(dazst)
      return
      end

