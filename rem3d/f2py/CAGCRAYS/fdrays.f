      subroutine fdrays(delta,npps,p,t,d,dddp,q,turn,
     #           pr,tr,dddpr,qr,turad,maxray,nray)
c
c---- routine to interpolate p,delta,travel time,and dd/dp tables
c     dimensioned for 5 possible arrivals
c     input:  range  [rad]
c             npps    =  number of points in p table to search
c     output: nrays   =  number of geometrical arrivals
c             pr,tr,dddpr  =  values for found rays
c
c      dimension p(1),t(1),d(1),dddp(1),q(1),turn(1)
      dimension p(maxray)
      dimension t(maxray)
      dimension d(maxray)
      dimension dddp(maxray)
      dimension q(maxray)
      dimension turn(maxray)
      dimension pr(maxray)
      dimension tr(maxray)
      dimension dddpr(maxray)
      dimension qr(maxray)
      dimension turad(maxray)
      data caustic /1.e-6/
c
      rtod=180./3.141592
      range=delta/rtod
      i1=0
      nray=0
    1 continue
      i1=i1+1
      if(i1+1.gt.npps) go to 9
      if(d(i1).gt.0..and.d(i1+1).gt.0.) then
        if((d(i1).le.range.and.d(i1+1).gt.range).or.
     #    (d(i1).gt.range.and.d(i1+1).le.range)) then
          if(nray+1.le.maxray) then
            nray=nray+1
            if(d(i1).ne.d(i1+1)) then
              dd=abs((range-d(i1))/(d(i1+1)-d(i1)))
            else
             dd=0.5
            endif
            pr(nray)=p(i1)+dd*(p(i1+1)-p(i1))
            tr(nray)=t(i1)+dd*(t(i1+1)-t(i1))
            dddpr(nray)=dddp(i1)+dd*(dddp(i1+1)-dddp(i1))
            if(abs(dddpr(nray)).lt.caustic) dddpr(nray)=caustic
            dddpr(nray)=dddpr(nray)*rtod*rtod
            qr(nray)=q(i1)+dd*(q(i1+1)-q(i1))
            turad(nray)=turn(i1)+dd*(turn(i1+1)-turn(i1))
          else
            write(6,"('too many rays exist for this distance - fdrays')")
            stop
          endif
        endif
      endif
      go to 1
    9 continue
      return
      end
