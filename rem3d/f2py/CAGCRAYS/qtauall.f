      subroutine qtauall(ilay,x,pray,isotropic,iraytype,qray,
     #          vpv,vph,vsv,vsh,eta,qmu,qkappa,s1,s2,s3,s4,s5,r)
      implicit double precision (a-h,o-z)
      logical isotropic
c
      call evem(ilay,x,isotropic,rho,vpv,vph,vsv,vsh,eta,qmu,qkappa,ierr)
c
c---- avoid the center of the earth
c
      if(x.lt.0.001d0) then
        x=0.001d0
      endif
c
      px2=(pray/x)**2
c
      if(iraytype.eq.1) then
        if(vsv.eq.0.d0) then
          qray=-1.0d0
          return
        endif
        if(isotropic) then
          qray=(1.d0/(vsv*vsv)) - px2
        else
          qray=(1.d0-vsh*vsh*px2)/(vsv*vsv)
        endif
      else if(iraytype.eq.2.or.iraytype.eq.3) then
        if(vsv.eq.0.d0) then
          if(iraytype.eq.3) then
            qray=-1.0d0
            return
          else
            qray=(1.d0/(vpv*vpv)) - px2
          endif
        else
          if(isotropic) then
            if(iraytype.eq.2) then
              qray=(1.d0/(vpv*vpv)) - px2
            else if(iraytype.eq.3) then
              qray=(1.d0/(vsv*vsv)) - px2
            endif
          else
            aa=0.5d0/(vsv*vsv)
            bb=0.5d0/(vpv*vpv)
            f=eta*(vph*vph-2.d0*vsv*vsv)
            cl2=2.d0*((vsv*vpv)**2)
            ac=((vph*vpv)**2)
            xl2=2.d0*vsv*vsv
c
            s1=aa+bb
            s2=aa-bb
            s3=(ac-f*f-f*xl2)/cl2
            s4=s3*s3-(vph*vph)/(vpv*vpv)
            s5=0.5d0*((1.d0+(vph*vph)/(vsv*vsv))/(vpv*vpv))-s1*s3
            r=dsqrt(s4*((pray/x)**4)+2.d0*s5*((pray/x)**2)+s2*s2)
            if(iraytype.eq.2) then
              qray=s1-s3*px2-r
            else if(iraytype.eq.3) then
              qray=s1-s3*px2+r
            endif
          endif
        endif
      endif
      return
      end
