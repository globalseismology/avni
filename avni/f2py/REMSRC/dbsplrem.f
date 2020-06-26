      subroutine dbsplrem(x,x0,xright,np,yarr,yarrp)
      implicit double precision (a-h,o-z)
      dimension yarr(np)
      dimension yarrp(np)
Cf2py intent(inout) x,x0,xright,np
Cf2py intent(out) yarr,yarrp

c
      dx=(xright-x0)/dfloat(np-1)
      nx=np-1
      xdiff=x-x0
      do i=1,np
        yarr(i)=0.d0
        yarrp(i)=0.d0
      enddo
c---- extrapolation to smaller values
      if(xdiff.le.0.) then
        value=1.0-xdiff/dx
        yarr(1)=value
        yarr(2)=xdiff/dx
        yarrp(1)=-1./dx
        yarrp(2)=1./dx
c---- extrapolation to bigger values
      else if(xdiff.ge.dfloat(nx)*dx) then
        value=1.0+(xdiff-dfloat(nx)*dx)/dx
        yarr(np)=value
        yarr(np-1)=-(xdiff-dfloat(nx)*dx)/dx
        yarrp(np)=1./dx
        yarrp(np-1)=-1./dx
c---- value falls between limits
      else
        interval=1+int((x-x0)/dx)
        if(dabs(xdiff-dfloat(nx)*dx).lt.0.0000001) interval=nx
        if(interval.ge.np) then
          stop 'something wrong in dbspl 1'
        endif
        if(interval.le.0) then
	  stop 'something wrong in dbspl 2'
        endif
        xd=x-x0-dfloat(interval-1)*dx
        h=1./dx
        hsq=1./dx**2
        hcu=1./dx**3
        xdsq=xd**2
        xdcu=xd**3
        xdsqp=2.*xd
        xdcup=3.*xdsq
c---- loop on basis cubic B-spline functions
        do i=0,nx
          value=0.
          valuep=0.
          if(i.eq.0) then
            if(interval.eq.1) then
              value=0.25*hcu*xdcu-1.5*h*xd+1.5
              valuep=0.25*hcu*xdcup-1.5*h
            else if(interval.eq.2) then
              value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
              valuep=-0.25*hcu*xdcup+0.75*hsq*xdsqp-0.75*h
            else
              value=0.
              valuep=0.
            endif
          else if(i.eq.1) then
            if(interval.eq.1) then
              value=-0.5*hcu*xdcu+1.5*h*xd
              valuep=-0.5*hcu*xdcup+1.5*h
            else if(interval.eq.2) then
              value=0.75*hcu*xdcu-1.5*hsq*xdsq+1.
              valuep=0.75*hcu*xdcup-1.5*hsq*xdsqp
              if(nx.eq.2) value=0.5*hcu*xdcu-1.5*hsq*xdsq+1.
              if(nx.eq.2) valuep=0.5*hcu*xdcup-1.5*hsq*xdsqp
            else if(interval.eq.3) then
              value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
              valuep=-0.25*hcu*xdcup+0.75*hsq*xdsqp-0.75*h
            else
              value=0.
              valuep=0.
            endif
          else if(i.gt.1.and.i.lt.nx-1) then
            if(interval.eq.i-1) then
              value=0.25*hcu*xdcu
              valuep=0.25*hcu*xdcup
            else if(interval.eq.i) then
              value=-0.75*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
              valuep=-0.75*hcu*xdcup+0.75*hsq*xdsqp+0.75*h
            else if(interval.eq.i+1) then
              value=0.75*hcu*xdcu-1.5*hsq*xdsq+1.
              valuep=0.75*hcu*xdcup-1.5*hsq*xdsqp
            else if(interval.eq.i+2) then
              value=-0.25*hcu*xdcu+0.75*hsq*xdsq-0.75*h*xd+0.25
              valuep=-0.25*hcu*xdcup+0.75*hsq*xdsqp-0.75*h
            else
              value=0.
              valuep=0.
            endif
          else if(i.eq.nx-1) then
            if(interval.eq.nx-2) then
              value=0.25*hcu*xdcu
              valuep=0.25*hcu*xdcup
            else if(interval.eq.nx-1) then
              value=-0.75*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
              valuep=-0.75*hcu*xdcup+0.75*hsq*xdsqp+0.75*h
            else if(interval.eq.nx) then
              value=0.5*hcu*xdcu-1.5*hsq*xdsq+1.
              valuep=0.5*hcu*xdcup-1.5*hsq*xdsqp
            else
              value=0.
              valuep=0.
            endif
          else if(i.eq.nx) then
            if(interval.eq.nx-1) then
              value=0.25*hcu*xdcu
              valuep=0.25*hcu*xdcup
            else if(interval.eq.nx) then
              value=-0.25*hcu*xdcu+0.75*hsq*xdsq+0.75*h*xd+0.25
              valuep=-0.25*hcu*xdcup+0.75*hsq*xdsqp+0.75*h
            else
              value=0.
              valuep=0.
            endif
          endif
          yarr(i+1)=value/1.5
          yarrp(i+1)=valuep/1.5
        enddo
      endif
      return
      end
