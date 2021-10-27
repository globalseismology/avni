      subroutine reademfl(lu,ml,nl,xb,xt,icm,ici,nlev,iflu,iani,ierr)
      implicit double precision (a-h,o-z)
      real*8 xb(ml)
      real*8 xt(ml)
      integer iflu(ml)
      integer iani(ml)
      integer nlev(ml)
      common /plevel/iprtlv
c
      include 'emcommon.h'
c
      call bffi(lu,1,iversion,4,nr,irr,0)
c
      call bffi(lu,1,emtitle,80,ierr,nr,0)
c
c      write(6,"('file version:',i4)") iversion
c      write(6,"('earth model:',a)") emtitle(1:lnblnk(emtitle))
c
      call bffi(lu,1,numemlev,4,ierr,nr,0)
c
      nbytes=numemlev*8
      call bffi(lu,1,radlev,nbytes,ierr,nr,0)
      call bffi(lu,1,rholev,nbytes,ierr,nr,0)
      call bffi(lu,1,vpvlev,nbytes,ierr,nr,0)
      call bffi(lu,1,vsvlev,nbytes,ierr,nr,0)
      call bffi(lu,1,bigqklev,nbytes,ierr,nr,0)
      call bffi(lu,1,bigqmlev,nbytes,ierr,nr,0)
      call bffi(lu,1,vphlev,nbytes,ierr,nr,0)
      call bffi(lu,1,vshlev,nbytes,ierr,nr,0)
      call bffi(lu,1,etalev,nbytes,ierr,nr,0)
      call bffi(lu,1,vsisolev,nbytes,ierr,nr,0)
      call bffi(lu,1,vpisolev,nbytes,ierr,nr,0)
      call bffi(lu,1,elllev,nbytes,ierr,nr,0)
      call bffi(lu,1,eta2lev,nbytes,ierr,nr,0)
      call bffi(lu,1,gravlev,nbytes,ierr,nr,0)
c
      nbytes=numemlev*8*3     
      call bffi(lu,1,rhospl,nbytes,ierr,nr,0)
      call bffi(lu,1,vpvspl,nbytes,ierr,nr,0)
      call bffi(lu,1,vsvspl,nbytes,ierr,nr,0)
      call bffi(lu,1,bigqkspl,nbytes,ierr,nr,0)
      call bffi(lu,1,bigqmspl,nbytes,ierr,nr,0)
      call bffi(lu,1,vphspl,nbytes,ierr,nr,0)
      call bffi(lu,1,vshspl,nbytes,ierr,nr,0)
      call bffi(lu,1,etaspl,nbytes,ierr,nr,0)
      call bffi(lu,1,vsisospl,nbytes,ierr,nr,0)
      call bffi(lu,1,vpisospl,nbytes,ierr,nr,0)
      call bffi(lu,1,ellspl,nbytes,ierr,nr,0)
      call bffi(lu,1,eta2spl,nbytes,ierr,nr,0)
      call bffi(lu,1,gravspl,nbytes,ierr,nr,0)
c
      call bffi(lu,1,numemreg,4,ierr,nr,0)
c
      nbytes=numemreg*8
      call bffi(lu,1,topemreg,nbytes,ierr,nr,0)
      call bffi(lu,1,botemreg,nbytes,ierr,nr,0)
c
      nbytes=numemreg*4
      call bffi(lu,1,itoplev,nbytes,ierr,nr,0)
      call bffi(lu,1,ibotlev,nbytes,ierr,nr,0)
c
      call getflpos(lu,iposition,ierr)
c
      if(iprtlv.gt.0) then
      write(6,"('number of levels defined:',i5)") numemlev
      write(6,"('found ',i3,' regions')") numemreg
      do iemreg=1,numemreg
        write(6,"(i3,':',i4,i4,f12.4,f12.4)") iemreg,ibotlev(iemreg),
     #     itoplev(iemreg),radlev(ibotlev(iemreg)),radlev(itoplev(iemreg))
      enddo
      endif
c
c---- find the inner-core boundary (top layer in inner core),
c---- the outer-core boundary (top layer in outer core),
c---- and the top of the mantle
c
      nic=0
      noc=0
      moho=0
      itopic=0
      itopoc=0
      itopmantle=0
      i670=0
      i400=0
      ici=0
      do iemreg=1,numemreg-1
        if(dabs(radlev(itoplev(iemreg))-1221.5d0).lt.25.0d0) then
          itopic=itoplev(iemreg)
          ici=iemreg
        endif
        if(dabs(radlev(itoplev(iemreg))-3480.d0).lt.25.d0) then
          itopoc=itoplev(iemreg)
          icm=iemreg
        endif
        if(dabs(radlev(itoplev(iemreg))-(6371.d0-670.d0)).lt.25.d0) then
          i670=itoplev(iemreg)
        endif
        if(dabs(radlev(itoplev(iemreg))-(6371.d0-400.d0)).lt.25.d0) then
          i400=itoplev(iemreg)
        endif
        if(vpvlev(itoplev(iemreg)).gt.7.5d0.and.
     #     vpvlev(ibotlev(iemreg+1)).lt.7.5d0) then
            itopmantle=itoplev(iemreg)
        endif
      enddo
      if(iprtlv.gt.0) then
      write(6,"('top level in inner core:',2i4)") nic,itopic
      write(6,"('top level in outer core:',2i4)") noc,itopoc
      write(6,"('top level in lower mantle (670):',i4,f10.4)") i670,6371.d0-radlev(i670)
      write(6,"('top level in transition zone (400):',i4,f10.1)") i400,6371.d0-radlev(i400)
      write(6,"('top level in mantle    :',2i4)") moho,itopmantle
      endif
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
      call getflpos(lu,iposition,ierr)
      if(iprtlv.gt.0) then
      write(6,"('file position at end of reademfl:',i6)") iposition
      endif
      return
      end
