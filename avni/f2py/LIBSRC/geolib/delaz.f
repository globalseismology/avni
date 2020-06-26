      subroutine delaz(eplat,eplong,stlat,stlong,delta,azep,azst)

      real*8 deplat,deplong,dstlat,dstlong,ddelta,dazep,dazst
c
      deplat=dble(eplat)
      deplong=dble(eplong)
      dstlat=dble(stlat)
      dstlong=dble(stlong)
c
      call ddelaz(deplat,deplong,dstlat,dstlong,ddelta,dazep,dazst)
c
      delta=sngl(ddelta)
      azep=sngl(dazep)
      azst=sngl(dazst)
      return
      end

