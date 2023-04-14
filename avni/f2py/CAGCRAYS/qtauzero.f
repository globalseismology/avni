      subroutine qtauzero(ilay,a1,b1,pray,isotropic,iraytype,z)
c  'z' is the zero of qtau (ie turning radius)  placed between 'a1'
c  and 'b1'. it is found by repeated linear interpolation between the
c  bounds of smaller and smaller intervals around the zero.
      implicit double precision(a-h,o-z)
      logical isotropic
      data re/1.d-14/
      a=a1
cc      call qtau(a,iq,fa)
c     write(6,"('radius:',f10.2)") a
      call qtau(ilay,a,pray,isotropic,iraytype,fa)
      if(fa.eq.0.d0) then
         z=a
         return
      end if
      b=b1
c      call qtau(b,iq,fb)
c     write(6,"('radius:',f10.2)") b
      call qtau(ilay,b,pray,isotropic,iraytype,fb)
      if(fb.eq.0.d0) then
         z=b
         return
      end if
c******* no zero crossing or not monotonic
      if(fa*fb.ge.0.d0) then
         z=0.d0
         ierr=1
         return
      end if
      ierr=0
      c=a
      fc=fa
      s=c
      fs=fc
c =======================================================
c--this block to be repeated until a satisfactory estimate is reached
   10 h=0.5d0*(b+c)
      t=dabs(h*re)
      if(dabs(h-b).le.t) then
        z=h
        return
      end if
      if(dabs(fb).gt.dabs(fc)) then
        y=b
        fy=fb
        g=b
        fg=fb
        s=c
        fs=fc
      else
        y=s
        fy=fs
        g=c
        fg=fc
        s=b
        fs=fb
      end if
      if(fy.eq.fs) then
        b=h
      else 
        e=(s*fy-y*fs)/(fy-fs)
        if(dabs(e-s).le.t) e=s+dsign(t,g-s)
        if((e-h)*(s-e).lt.0.d0) then
          b=h
        else
          b=e
        end if
      end if
c     write(6,"('radius:',f10.2)") b
      call qtau(ilay,b,pray,isotropic,iraytype,fb)
c      call qtau(b,iq,fb)
      if(fg*fb.ge.0.d0) then
        c=s
        fc=fs
      else
        c=g
        fc=fg
      end if
      goto 10
c ================================================================
      end
