      subroutine reademb(lu,ml,nl,xb,xt,nlev,iflu,iani,capom,capg,ierr)
      implicit double precision (a-h,o-z)
      real*8 xb(ml)
      real*8 xt(ml)
      integer iflu(ml)
      integer iani(ml)
      integer nlev(ml)
c
      include 'emcommonb.h'
c
      dimension y(3),y1(3),y2(3),y3(3),ydot(3)
c      parameter (capom=7.292115d-5)
c      parameter (capg=6.6723d-11)
      parameter (twopi=6.2831853072)
c
      read(lu,"(a)") emtitle
      read(lu,*) ifanis,tref,ifdeck
      read(lu,*) n,nic,noc,moho
c
      ios=0
      numemlev=0
      do while(ios.eq.0) 
        read(lu,*,iostat=ios) v1,v2,v3,v4,v5,v6,v7,v8,v9
        if(ios.eq.0) then
          numemlev=numemlev+1
        endif
        iemlev=numemlev
        radlev(iemlev)=v1*0.001d0
        rholev(iemlev)=v2*0.001d0
        vpvlev(iemlev)=v3*0.001d0
        vsvlev(iemlev)=v4*0.001d0
        bigqklev(iemlev)=v5
        bigqmlev(iemlev)=v6
        vphlev(iemlev)=v7*0.001d0
        vshlev(iemlev)=v8*0.001d0
        etalev(iemlev)=v9
        xc=v2*v3*v3
        xl=v2*v4*v4
        xa=v2*v7*v7
        xn=v2*v8*v8
        xf=v9*(xa-2.d0*xl)
        xkapa=(4.d0*xa+xc+4.d0*xf-4.d0*xn)/9.d0
        xmu=(xa+xc-2.d0*xf+5.0d0*xn+6.0d0*xl)/15.d0
        xa=xkapa+4.0d0*xmu/3.0d0
        xn=xmu
        vpisolev(iemlev)=dsqrt(xa/v2)*0.001d0
        vsisolev(iemlev)=dsqrt(xn/v2)*0.001d0
	bigklev(iemlev)=xkapa
      enddo
c      write(6,"('number of levels defined:',2i5)") n,numemlev
c
c---- spline the variables
c
      call drspln(1,numemlev,radlev,rholev,rhospl,work)
      call drspln(1,numemlev,radlev,vpvlev,vpvspl,work)
      call drspln(1,numemlev,radlev,vsvlev,vsvspl,work)
      call drspln(1,numemlev,radlev,bigqklev,bigqkspl,work)
      call drspln(1,numemlev,radlev,bigqmlev,bigqmspl,work)
      call drspln(1,numemlev,radlev,vphlev,vphspl,work)
      call drspln(1,numemlev,radlev,vshlev,vshspl,work)
      call drspln(1,numemlev,radlev,etalev,etaspl,work)
      call drspln(1,numemlev,radlev,vpisolev,vpisospl,work)
      call drspln(1,numemlev,radlev,vsisolev,vsisospl,work)
      call drspln(1,numemlev,radlev,bigklev,bigkspl,work)
c
c---- find the number of regions
c
      numemreg=0
      radprev=radlev(1)
      ibotprev=1
      do iemlev=2,numemlev
        if(radlev(iemlev).eq.radprev) then
          numemreg=numemreg+1
          ibotlev(numemreg)=ibotprev
          itoplev(numemreg)=iemlev-1
          ibotprev=iemlev
        endif
        radprev=radlev(iemlev)
      enddo
      numemreg=numemreg+1
      ibotlev(numemreg)=ibotprev
      itoplev(numemreg)=numemlev
c
c      write(6,"('found ',i3,' regions')") numemreg
c      do iemreg=1,numemreg
c        write(6,"(i3,':',i4,i4,f12.4,f12.4)") iemreg,ibotlev(iemreg),
c     #     itoplev(iemreg),radlev(ibotlev(iemreg)),radlev(itoplev(iemreg))
c      enddo
c
c---- find the inner-core boundary (top layer in inner core),
c---- the outer-core boundary (top layer in outer core),
c---- and the top of the mantle
c
      itopic=0
      itopoc=0
      itopmantle=0
      i670=0
      i400=0
      do iemreg=1,numemreg-1
        if(dabs(radlev(itoplev(iemreg))-1221.5d0).lt.25.0d0) then
          itopic=itoplev(iemreg)
        endif
        if(dabs(radlev(itoplev(iemreg))-3480.d0).lt.25.d0) then
          itopoc=itoplev(iemreg)
        endif
        if(dabs(radlev(itoplev(iemreg))-(6371.d0-670.d0)).lt.25.d0) then
          i670=itoplev(iemreg)
        endif
        if(dabs(radlev(itoplev(iemreg))-(6371.d0-400.d0)).lt.25.d0) then
          i400=itoplev(iemreg)
        endif
        if(vpvlev(itoplev(iemreg)).gt.7.5d0.and.
     #     vpvlev(ibotlev(iemreg+1)).lt.7.5d0) then
c        write(6,"(i4,2f10.4)") iemreg,vpvlev(itoplev(iemreg)),vpvlev(ibotlev(iemreg+1))
            itopmantle=itoplev(iemreg)
        endif
      enddo
C     write(6,"('top level in inner core:',2i4)") nic,itopic
C     write(6,"('top level in outer core:',2i4)") noc,itopoc
C     write(6,"('top level in lower mantle (670):',i4,f10.4)") i670,6371.d0-radlev(i670)
C     write(6,"('top level in transition zone (400):',i4,f10.1)") i400,6371.d0-radlev(i400)
C     write(6,"('top level in mantle    :',2i4)") moho,itopmantle
c
      do iemreg=1,numemreg
        xb(iemreg)=radlev(ibotlev(iemreg))
        xt(iemreg)=radlev(itoplev(iemreg))
        iflu(iemreg)=0
        if(vsvlev(ibotlev(iemreg)).eq.0.d0) then
          iflu(iemreg)=1
        endif
        iani(iemreg)=0
        if(vpvlev(itoplev(iemreg)).ne.vphlev(itoplev(iemreg))) then
          iani(iemreg)=1
        endif
        nlev(iemreg)=itoplev(iemreg)-ibotlev(iemreg)+1
      enddo
      nl=numemreg
c
c---- calculate the ellipticity function
c
      y(1)=1.d0
      y(2)=0.d0
      y(3)=rholev(1)
      rr=1.d-5
      icount=0
c
      do 100 ilev=2,numemlev
        rstep=radlev(ilev)-radlev(ilev-1)
        if(dabs(rstep).gt.1.d-5) then
          go to 110
        endif
        elllev(ilev)=elllev(ilev-1)
        eta2lev(ilev)=eta2lev(ilev-1)
        gravlev(ilev)=gravlev(ilev-1)
        go to 100
c
  110   continue
        do i=1,3
          y1(i)=y(i)
        enddo
        rr1=rr
        if(rr+rstep.gt.radlev(ilev)) then
          rstep=radlev(ilev)-rr
        endif
        iback=1
        go to 1100
c
 1201   continue
        do i=1,3
          y2(i)=y(i)
        enddo
        rr2=rr
  901   continue
        rstep=rstep*0.5d0
        do i=1,3
          y(i)=y1(i)
        enddo
        rr=rr1
        iback=2
        go to 1100
c
 1202   continue
        do i=1,3
          y3(i)=y(i)
        enddo
        rr3=rr
        iback=3
        go to 1100
c
 1203   continue
        do i=1,3
          if(dabs(y(i)-y2(i)).gt.0.5d-5) then
            go to 601
          endif
        enddo
        go to 701
c
  601   continue
        do i=1,3
          y2(i)=y3(i)
        enddo
        rr2=rr3
        go to 901
c
  701   continue
        if(dabs(rr-radlev(ilev)).lt.1.d-5) then 
          go to 1300
        endif
        rstep=4.d0*rstep
        go to 110
c
 1300   continue
        elllev(ilev)=y(1)
        eta2lev(ilev)=rr*y(2)/y(1)
        gravlev(ilev)=4.*rr*y(3)/3.
        go to 100
c
 1100   continue
        icount=icount+1
        k=krunge(3,y,ydot,rr,rstep)
        if(k.ne.1) then
          go to (1201,1202,1203), iback
        endif
        ydot(1)=y(2)
        rhot=drsple(1,numemlev,radlev,rholev,rhospl,rr)
        ydot(3)=3.d0*(rhot-y(3))/rr
        ydot(2)=-2.d0*(ydot(3)*y(1)+3.d0*rhot*y(2))/(y(3)*rr)
        go to 1100
c
  100 continue
c      write(6,"('ell,eta,rhob',3g15.5)") elllev(numemlev),eta2lev(numemlev),y(3)
c
c---- scaling for ellipticity
c
      const=2.5d0*(capom**2)*3.d0/
     #         (2.d0*twopi*capg*y(3)*1000.d0*elllev(numemlev)*(eta2lev(numemlev)+2.d0))
c      write(6,"('const',g15.5)") const 
      do ilev=1,numemlev
        elllev(ilev)=elllev(ilev)*const
      enddo
c
c---- spline the gravity also
c
      do ilev=1,numemlev
        gravlev(ilev)=gravlev(ilev)*capg*twopi*0.5*1.e6
      enddo
c
c---- spline the ellipticity, eta parameter, and gravity
c
      call drspln(1,numemlev,radlev,elllev,ellspl,work)
      call drspln(1,numemlev,radlev,eta2lev,eta2spl,work)
      call drspln(1,numemlev,radlev,gravlev,gravspl,work)
c      write(6,"('ell:',g15.5)") elllev(numemlev),1.d0/elllev(numemlev)
c
c---- integrate to get the pressure, starting from CoE.
c
      prelev(1)=0.
      do ilev=2,numemlev 
        rstep=radlev(ilev)-radlev(ilev-1)
	if(dabs(rstep).gt.1.d-5) then
	  prod1=rholev(ilev-1)*gravlev(ilev-1)
	  prod2=rholev(ilev-1)*gravspl(1,ilev-1)+
     #          rhospl(1,ilev-1)*gravlev(ilev-1)
          prod3=rholev(ilev-1)*gravspl(2,ilev-1)+
     #          rhospl(1,ilev-1)*gravspl(1,ilev-1)+
     #          rhospl(2,ilev-1)*gravlev(ilev-1)
          prod4=rholev(ilev-1)*gravspl(3,ilev-1)+
     #          rhospl(1,ilev-1)*gravspl(2,ilev-1)+
     #          rhospl(2,ilev-1)*gravspl(1,ilev-1)+
     #          rhospl(3,ilev-1)*gravlev(ilev-1)
          prod5=rhospl(1,ilev-1)*gravspl(3,ilev-1)+
     #          rhospl(2,ilev-1)*gravspl(2,ilev-1)+
     #          rhospl(3,ilev-1)*gravspl(1,ilev-1)
          prod6=rhospl(2,ilev-1)*gravspl(3,ilev-1)+
     #          rhospl(3,ilev-1)*gravspl(2,ilev-1)
          prod7=rhospl(3,ilev-1)*gravspl(3,ilev-1)
c
          sum7=prod1*rstep+prod2*(rstep**2)/2.+prod3*(rstep**3)/3.+
     #         prod4*(rstep**4)/4.+prod5*(rstep**5)/5. +
     #         prod6*(rstep**6)/6.+prod7*(rstep**7)/7.
          sum7=sum7*1.d6
          prelev(ilev)=prelev(ilev-1)+sum7
	else
	  prelev(ilev)=prelev(ilev-1)
        endif
      enddo
      prelev(1)=prelev(numemlev)
      do ilev=2,numemlev 
        prelev(ilev)=prelev(1)-prelev(ilev)
      enddo
      call drspln(1,numemlev,radlev,prelev,prespl,work)
      return
      end
c
      function krunge(n,y,f,x,h)
      implicit double precision (a-h,o-z)
      dimension phi(6),savey(6),y(6),f(6)
      data m/0/
      save
c
      m = m + 1
      goto(1,2,3,4,5), m
c***
1     krunge=1
      return
c***
2     do 22 j=1,n
      savey(j) = y(j)
      phi(j)   = f(j)
22    y(j) = savey(j)+0.5d0*h*f(j)
      x = x + 0.5d0*h
      krunge = 1
      return
c***
3     do 33 j=1,n
      phi(j) = phi(j) + 2.d0*f(j)
33    y(j)   = savey(j)+0.5d0*h*f(j)
      krunge = 1
      return
c***
4     do 44 j=1,n
      phi(j) = phi(j)+2.0d0*f(j)
44    y(j)   = savey(j)+h*f(j)
      x      = x + 0.5d0*h
      krunge = 1
      return
c***
5     do 55 j=1,n
55    y(j) = savey(j)+(phi(j)+f(j))*h/6.d0
      m    = 0
      krunge = 0
      return
      end

