      subroutine dxt(ilay,x,pray,isotropic,iraytype,d,f)
      implicit double precision (a-h,o-z)
      logical isotropic
      dimension v(3),dv(3),ds(5),bigv(3,3),bigd(3,3)
c
      call evem(ilay,x,isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
      call evemdr(ilay,x,isotropic,rhodr,vpvdr,vphdr,vsvdr,vshdr,
     #            etadr,qmudr,qkappadr)
c
      d=0.d0
      f=0.d0
      if(x.lt.1.d-3) then
         return
      endif
      px2=(pray/x)**2
c
      if(iraytype.eq.1) then
        if(vsv.eq.0.d0) then
        else
          if(.not.isotropic) then
            dvsh=vshdr/vsh
            dvsv=vsvdr/vsv
            top=px2*vsh*vsh/(vsv*vsv*pray)
            bot=((px2/x)-dvsv/(vsh*vsh)+px2*(dvsv-dvsh))*(vsh*vsh)/(vsv*vsv)
          else
            top=px2/pray
            bot=px2/x-vsvdr/(vsv**3)
          endif
          d=top/bot
          f=0.5d0*d*dsqrt(2.d0*bot)
        endif
      else if(vsv.eq.0.d0) then
        top=px2/pray
        bot=px2/x-vpvdr/(vpv**3)
        d=top/bot
        f=0.5d0*d*dsqrt(2.d0*bot)
      else
        v(2)=vpv
        v(3)=vsv
        dv(2)=vpvdr
        dv(3)=vsvdr
        a=dv(2)/(v(2)**3)
        b=dv(3)/(v(3)**3)
        ds(1)=-a-b
        ds(2)=a-b
        if(isotropic) then
          top=2.d0*px2/pray
          bot=ds(1)+2.d0*px2/x
          if(iraytype.eq.2) bot=bot-ds(2)
          if(iraytype.eq.3) bot=bot+ds(2)
        else
          xa=rho*vph*vph
          xc=rho*vpv*vpv
          xl=rho*vsv*vsv
          xf=eta*(xa-2.0d0*xl)
          v(1)=vph
          dv(1)=vpvdr
          deta=etadr
          do 1 i=1,3
          do 1 j=1,3
    1       bigv(i,j)=v(i)*v(i)/(v(j)*v(j))
          do 2 i=1,3
          do 2 j=1,3
    2       bigd(i,j)=2.d0*bigv(i,j)*(dv(i)/v(i)-dv(j)/v(j))
          aa=0.5d0/(vsv*vsv)
          bb=0.5d0/(vpv*vpv)
          s1=aa+bb
          s2=aa-bb
          s3=(xa*xc-xf*xf-2.d0*xf*xl)/(2.d0*xc*xl)
          s4=s3*s3-xa/xc
          s5=0.5d0*(1.d0+xa/xl)/(vpv*vpv)-s1*s3
          r=sqrt((s4*px2+2.d0*s5)*px2+s2*s2)
c
c*** ds(1--5) are derivs of s() wrt x (similarly with rdx and r)
c
          b1=-0.5d0*eta*bigd(1,3)+deta*(2.d0-0.5d0*bigv(1,3))
          a2=-1.d0+eta*(2.d0-0.5d0*bigv(1,3))
          b2=eta*bigd(1,2)+deta*bigv(1,2)
          a3=-2.d0*eta*deta*bigv(3,2)
          a4=2.d0*(1.d0-eta)*(eta*bigd(3,2)+deta*bigv(3,2))
          ds(3)=0.5d0*bigd(1,3)+eta*bigv(1,2)*b1+a2*b2+a3+a4
          ds(4)=2.d0*s3*ds(3)-bigd(1,2)
          ds(5)=(0.5d0*bigd(1,3)-dv(2)*(1.d0+bigv(1,3))/vpv)/(vpv*vpv)
     +          -s1*ds(3)-s3*ds(1)
          a=0.5d0*px2*px2*(ds(4)-4.d0*s4/x)
          b=px2*(ds(5)-2.d0*s5/x)+s2*ds(2)
          rdx=(a+b)/r
          rdp=2.d0*px2*(s4*px2+s5)/(r*pray)
          top=2.d0*s3*px2/pray
          bot=ds(1)+px2*(2.d0*s3/x-ds(3))
          if(iraytype.eq.2)then
            top=top+rdp
            bot=bot-rdx
          else
            bot=bot+rdx
            top=top-rdp
          end if
        endif
        d=top/bot
        f=0.5d0*d*dsqrt(bot)
      end if
      return
      end
