      subroutine fqs(ilay,x,xnorm,iturn,fxtpp,dxtdp,pray,isotropic,
     #iraytype,nvals,vals)
      implicit double precision (a-h,o-z)
      logical isotropic
      real*8 vals(4)
      common /plevel/iprtlv
c
      y=x*x+xnorm
      call qtauall(ilay,y,pray,isotropic,iraytype,qray,
     #          vpv,vph,vsv,vsh,eta,qmu,qkappa,s1,s2,s3,s4,s5,r)
      if(iprtlv.gt.4) then
        write(6,"('radius:',3g15.5)") x,xnorm,y
        write(6,"(14e12.3)") qray,vpv,vph,vsv,vsh,eta,qmu,qkappa,
     #        s1,s2,s3,s4,s5,r
      endif
      qrayrt=dsqrt(qray)
c
      p2=pray/(y*y)
      py2=p2*pray
c
c---- kernel for delta
c
      q1=p2/qrayrt
      f=1.0d0
      if(isotropic) then
      else
        if(iraytype.eq.1) then
          f=(vsh*vsh)/(vsv*vsv)
        else
          if(vsv.ne.0.d0) then
            a=(s4*py2+s5)/r
            if(iraytype.eq.3) then
              a=-a
            endif
            f=s3+a
          endif
        endif
        q1=q1*f
      endif
      vals(1)=2.d0*x*q1
c
c---- kernel for ttime
c
      if(isotropic) then
        if(iraytype.eq.2) then
          q2=1.d0/(vpv*vpv*qrayrt)
        else
          q2=1.d0/(vsv*vsv*qrayrt)
        endif
      else
        if(iraytype.eq.1) then
          q2=1.d0/(vsv*vsv*qrayrt)
        else
          b=(s5*py2+s2*s2)/r
          if(iraytype.eq.2) then
            q2=(s1-b)/qrayrt
          else if(iraytype.eq.3) then
            q2=(s1+b)/qrayrt
          endif
        endif
      endif
      vals(2)=2.d0*x*q2
c      
c---- kernel for dDelta/dp
c
      q3=f*(1.d0+pray*q1/qrayrt)/(y*y*qrayrt)
      if(isotropic) then
      else
        if(iraytype.eq.1) then
        else
          a=2.d0*py2*(s4*s2*s2-s5*s5)/(qrayrt*(r**3)*(y**2))
          if(iraytype.eq.3) a=-a
          q3=q3+a
        endif
      endif
      if(iturn.ge.0) then
        q3=q3-fxtpp*dxtdp/(2.d0*x**3)
      endif
      vals(3)=2.d0*x*q3
c
c---- kernel for tstar (vsv,vpv should properly be anisotropic wavespeed)
c
      if(iraytype.eq.1.or.iraytype.eq.3) then
        qbeta=1.d0/qmu
        vals(4)=vals(2)*qbeta
      else
        if(vsv.lt.0.001d0) then
          qalpha=1.d0/qkappa
        else
          if(isotropic) then
            qalpha=1.d0/qkappa+(4.d0/3.d0)*(1.d0/qmu-1.d0/qkappa)*((vsv/vpv)**2)
          else
            call evem(ilay,y,.true.,rh1,vp1,vp2,vs1,vs2,et1,q1,q2,ierr2)
            qalpha=1.d0/qkappa+(4.d0/3.d0)*(1.d0/qmu-1.d0/qkappa)*((vs1/vp1)**2)
          endif
        endif
        vals(4)=vals(2)*qalpha
      endif
c
      return
      end
