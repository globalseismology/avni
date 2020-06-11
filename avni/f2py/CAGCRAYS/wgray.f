      subroutine wgray(pray,ilay,xbot,xtop,numlev,isotropic,iraytype,
     #                 delta,ttime,dddp,tstar,xturn,
     #                 maxpts,npts,runrad,rundel,runtim)
      implicit double precision (a-h,o-z)
      logical isotropic
      real*8 runrad(maxpts)
      real*8 rundel(maxpts)
      real*8 runtim(maxpts)
      common /plevel/iprtlv
c
      delta=0.d0
      ttime=0.d0
      dddp=0.d0
      tstar=0.d0
      npts=0
c
c---- first check if ray turns in this layer
c
      xturn=-dfloat(ilay)
      call qtau(ilay,xtop,pray,isotropic,iraytype,qray)
c
      if(qray.lt.0.d0) then
        return
      endif
c
      call qtau(ilay,xbot,pray,isotropic,iraytype,qray)
      if(qray.lt.0.d0) then
        call qtauzero(ilay,xbot,xtop,pray,isotropic,iraytype,xturn)
        xnorm=xturn
        iturn=0
      else
        xturn=0.d0
        xnorm=xbot
        iturn=-1
      endif
      if(iturn.ge.0) then
        call dxt(ilay,xturn,pray,isotropic,iraytype,dxtdp,fxtpp)
      else
        dxtdp=0.d0
        fxtpp=0.d0
      endif
c
      nl=2*numlev-1
      xdec=(xtop-xbot)/dfloat(nl)
      x1=xtop-xnorm
c
c---- loop on integration 
c
      npts=1
      rundel(npts)=0.d0
      runrad(npts)=xtop
c
      numint=4
      x2=x1
      do while (x2.gt.0.0001d0)
        x2=x1
        if(x2.gt.0.d0) then
          if(x1.lt.250.d0.and.xturn.gt.0.0d0) then
            x1=x1-0.25d0*xdec
          else
            x1=x1-xdec
          endif
          if(x1.lt.1.d-3) then
            x1=0.d0
          endif
            xx1=dsqrt(x1)
            xx2=dsqrt(x2)
            call wgint(ilay,xx1,xx2,xnorm,iturn,fxtpp,dxtdp,
     #               pray,isotropic,iraytype,numint,
     #               deltainc,ttimeinc,dddpinc,tstarinc)
            delta=delta+deltainc
            ttime=ttime+ttimeinc
            dddp=dddp+dddpinc
            tstar=tstar+tstarinc
c
            npts=npts+1
            runrad(npts)=x1+xnorm
            rundel(npts)=delta
            runtim(npts)=ttime
            if(ttimeinc.lt.1.d-3)then
            write(6,"('rad,del,tim',6g15.5)")
     #        x1+xnorm,delta,ttime,ttimeinc,xx1,xx2
            endif
c
        endif
      enddo
c
c---- here we modify dddp if the ray turns in this layer
c
      if(iturn.ge.0) then
        dddp=dddp-dxtdp*fxtpp/(dsqrt(xtop-xnorm))
      endif
      return
      end
