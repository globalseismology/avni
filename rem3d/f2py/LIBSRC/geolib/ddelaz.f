      subroutine ddelaz(eplat,eplong,stlat,stlong,delta,azep,azst)
c
      implicit double precision (a-h,o-z)
      data geoco,hpi,twopi,rad,reprad/0.993277d0,1.57079632675d0,
     16.283185307d0,.0174532925d0,57.2957795d0/
c
      dtan(x)=dsin(x)/dcos(x)
      darcos(x)=datan2(dsqrt(1.d0-x*x),x)
c
      if(dabs(eplat).ne.90.0d0) then
	el=datan(geoco*dtan(eplat*rad))
      else
        el=eplat*rad
      endif
      el=hpi-el
c
      if(dabs(stlat).ne.90.0d0) then
	stl=datan(geoco*dtan(stlat*rad))
      else
        stl=stlat*rad
      endif
      stl=hpi-stl
c
      elon=eplong*rad
      slon=stlong*rad
      as=dcos(stl)
      bs=dsin(stl)
      cs=dcos(slon)
      ds=dsin(slon)
      a=dcos(el)
      b=dsin(el)
      c=dcos(elon)
      d=dsin(elon)
      cdel=a*as+b*bs*(c*cs+d*ds)
      if(dabs(cdel).gt.1.d0) cdel=dsign(1.d0,cdel)
      delt=darcos(cdel)
      delta=delt*reprad
      sdel=dsin(delt)
      caze=(as-a*cdel)/(sdel*b)
      if(dabs(caze).gt.1.d0) caze=dsign(1.d0,caze)
      aze=darcos(caze)
      if(bs.gt.0.d0) cazs=(a-as*cdel)/(bs*sdel)
      if(bs.eq.0.d0) cazs=dsign(1.d0,cazs)
      if(dabs(cazs).gt.1.d0) cazs=dsign(1.d0,cazs)
      azs=darcos(cazs)
      dif=ds*c-cs*d
      if(dif.lt.0.d0) aze=twopi-aze
      azep=reprad*aze
      if(dif.gt.0.d0) azs=twopi-azs
      azst=reprad*azs
      return
      end
