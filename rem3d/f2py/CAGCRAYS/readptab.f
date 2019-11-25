      subroutine readptab(lu,initpos,rad4,isl,ifanis,ierr)
c----------------------------------------------------------------------------
c---- This subroutine reads the ray-parameter tables and interpolates
c---- a new table for the appropriate source depth
c
      common /plevel/ iprtlv
c
      parameter (maxrec=500)
      dimension idum(10)
      parameter (twopi=6.2831853072)
c     
      include 'rays.h'
c
      dimension ilayer(maxrec)
      dimension ianilay(maxrec)
      dimension ipartial(maxrec)
      dimension itypwave(maxrec)
      dimension toprad(maxrec)
      dimension botrad(maxrec)
      dimension nplay(maxrec)
c
      dimension t1(maxp)
      dimension t2(maxp)
      dimension d1(maxp)
      dimension d2(maxp)
      dimension dddp1(maxp)
      dimension dddp2(maxp)
      dimension tstar1(maxp)
      dimension tstar2(maxp)
c-----------------------------------------------------------------------------
c
c---- read the number of p values
c
      call bffi(lu,1,np,4,ierr,nr,initpos+1)
      call bffi(lu,1,pfull,4*np,ierr,nr,0)
      call bffi(lu,1,nlaytabs,4,ierr,nr,0)
c
      if(iprtlv.gt.0) then
        write(6,"('first and last p', i5,2f10.4)") np,pfull(1),pfull(np)
        write(6,"('number of layer tables',i5)") nlaytabs
        write(6,"('isl, ifanis, rad4:',2i5,f10.4)") isl,ifanis,rad4
      endif
c      
c----bug fixed by Raj, 2014. Was reading lu=1 everytime      
C      call getflpos(1,itabpos,ierr)
      call getflpos(lu,itabpos,ierr)
      nbytes=np*4
c
c---- read in the full layers
c
      nlayers=0
      do ilaytabs=1,nlaytabs
        ipos=itabpos+1+(ilaytabs-1)*(4*7+nbytes*6)
        call bffi(lu,1,ilay,4,ierr,nr,ipos)
        call bffi(lu,1,iani,4,ierr,nr,0)
        call bffi(lu,1,iraytype,4,ierr,nr,0)
        call bffi(lu,1,isub,4,ierr,nr,0)
        call bffi(lu,1,top4,4,ierr,nr,0)
        call bffi(lu,1,bot4,4,ierr,nr,0)
        call bffi(lu,1,npout,4,ierr,nr,0)
            if(iprtlv.gt.0) then
              write(6,"(5i5,2f12.5)") 
     #              iani,ilay,isub,iraytype,npout,top4,bot4
            endif
        if(ifanis.eq.0.and.iani.eq.1) then
        else
          if(isub.eq.1) then
            call bffi(lu,1,tfull(1,ilay,iraytype),nbytes,ierr,nr,0)
            call bffi(lu,1,dfull(1,ilay,iraytype),nbytes,ierr,nr,0)
            call bffi(lu,1,dddpfull(1,ilay,iraytype),nbytes,ierr,nr,0)
            call bffi(lu,1,tstarfull(1,ilay,iraytype),nbytes,ierr,nr,0)
            call bffi(lu,1,xturnfull(1,ilay,iraytype),nbytes,ierr,nr,0)
            if(iprtlv.gt.0) then
              write(6,"(5i5,2f12.5)") 
     #              iani,ilay,isub,iraytype,npout,top4,bot4
            endif
            npfull(ilay,iraytype)=npout
            if(iraytype.eq.1.and.iani.eq.0) then
              npfull(ilay,iraytype+2)=npfull(ilay,iraytype)
              do i=1,np
                tfull(i,ilay,3)=tfull(i,ilay,1)
                dfull(i,ilay,3)=dfull(i,ilay,1)
                dddpfull(i,ilay,3)=dddpfull(i,ilay,1)
                tstarfull(i,ilay,3)=tstarfull(i,ilay,1)
                xturnfull(i,ilay,3)=xturnfull(i,ilay,1)
              enddo
            endif
            topfull(ilay)=top4
            botfull(ilay)=bot4
          endif
        endif
        ilayer(ilaytabs)=ilay
        ianilay(ilaytabs)=iani
        ipartial(ilaytabs)=isub
        itypwave(ilaytabs)=iraytype
        toprad(ilaytabs)=top4
        botrad(ilaytabs)=bot4
        nplay(ilaytabs)=npout
        nlayers=max0(nlayers,ilay)
      enddo
c
c---- now find the layers that have a top above and below the source depth
c
      nrt=2
      if(ifanis.eq.1) then
        nrt=3
      endif
      do irt=1,nrt
        iabove=0
        ibelow=0
        if(ifanis.eq.1) then
          do il=1,nlaytabs
            if(ilayer(il).eq.isl.and.
     #         itypwave(il).eq.irt.and.ianilay(il).eq.1) then
              if(toprad(il).ge.rad4) then
                iabove=il
              endif
              if(toprad(il).lt.rad4.and.ibelow.eq.0) then
                ibelow=il
              endif
            endif
          enddo
        endif
        if(iabove.eq.0) then
          do il=1,nlaytabs
            if(ilayer(il).eq.isl.and.
     #       ((itypwave(il).eq.irt.and.ianilay(il).eq.0).or.
     #       (itypwave(il)+2.eq.irt.and.ianilay(il).eq.0))) then
              if(toprad(il).ge.rad4) then
                iabove=il
              endif
              if(toprad(il).lt.rad4.and.ibelow.eq.0) then
                ibelow=il
              endif
            endif
          enddo
        endif
        if(iprtlv.gt.0) then
        write(6,"('irt, iabove, ibelow:',3i4)") irt,iabove,ibelow
        write(6,"('toptop,topbot,source',3f10.2,2i4)") 
     #      toprad(iabove),toprad(ibelow),rad4,ilayer(iabove),ilayer(ibelow)
        write(6,"('bottop,botbot,source',3f10.2,2i4)") 
     #      botrad(iabove),botrad(ibelow),rad4,ilayer(iabove),ilayer(ibelow)
        endif
c
        if(iabove.gt.0) then
          ipos=itabpos+1+(iabove-1)*(4*7+nbytes*6)
          call bffi(lu,1,ilay,4,ierr,nr,ipos)
          call bffi(lu,1,iani,4,ierr,nr,0)
          call bffi(lu,1,iraytype,4,ierr,nr,0)
          call bffi(lu,1,isub,4,ierr,nr,0)
          call bffi(lu,1,top4,4,ierr,nr,0)
          call bffi(lu,1,bot4,4,ierr,nr,0)
          call bffi(lu,1,npout,4,ierr,nr,0)
          call bffi(lu,1,t1(1),nbytes,ierr,nr,0)
          call bffi(lu,1,d1(1),nbytes,ierr,nr,0)
          call bffi(lu,1,dddp1(1),nbytes,ierr,nr,0)
          call bffi(lu,1,tstar1(1),nbytes,ierr,nr,0)
        else
          write(6,"('something went wrong')") 
          pause
        endif
c
        if(ibelow.gt.0) then
          ipos=itabpos+1+(ibelow-1)*(4*7+nbytes*6)
          call bffi(lu,1,ilay,4,ierr,nr,ipos)
          call bffi(lu,1,iani,4,ierr,nr,0)
          call bffi(lu,1,iraytype,4,ierr,nr,0)
          call bffi(lu,1,isub,4,ierr,nr,0)
          call bffi(lu,1,top4,4,ierr,nr,0)
          call bffi(lu,1,bot4,4,ierr,nr,0)
          call bffi(lu,1,npout,4,ierr,nr,0)
          call bffi(lu,1,t2(1),nbytes,ierr,nr,0)
          call bffi(lu,1,d2(1),nbytes,ierr,nr,0)
          call bffi(lu,1,dddp2(1),nbytes,ierr,nr,0)
          call bffi(lu,1,tstar2(1),nbytes,ierr,nr,0)
          deltarad=toprad(iabove)-toprad(ibelow)
          delsrad=rad4-toprad(ibelow)
          frac=delsrad/deltarad
          nps=min0(nplay(ibelow),nplay(iabove))
          do i=1,nps
            tbelow(i,irt)=t2(i)+frac*(t1(i)-t2(i))
            dbelow(i,irt)=d2(i)+frac*(d1(i)-d2(i))
            dddpbelow(i,irt)=dddp2(i)+frac*(dddp1(i)-dddp2(i))
            tstarbelow(i,irt)=tstar2(i)+frac*(tstar1(i)-tstar2(i))
          enddo
          npbelow(irt)=nps
          if(iprtlv.gt.0) then
          write(6,"('frac:',f10.3)") frac
          endif
        else
          deltarad=toprad(iabove)-botfull(isl)
          delsrad=rad4-botfull(isl)
          frac=0.
          if(delsrad.gt.0.) then
            frac=delsrad/deltarad
          endif
          nps=min0(nplay(iabove),npfull(isl,irt))
          do i=1,nps
            tbelow(i,irt)=frac*t1(i)
            dbelow(i,irt)=frac*d1(i)
            dddpbelow(i,irt)=frac*dddp1(i)
            tstarbelow(i,irt)=frac*tstar1(i)
          enddo
          npbelow(irt)=nps
        endif
        nps=npfull(isl,irt)
        do i=1,nps
          if(i.le.npbelow(irt)) then
            tabove(i,irt)=tfull(i,isl,irt)-tbelow(i,irt)
            dabove(i,irt)=dfull(i,isl,irt)-dbelow(i,irt)
            dddpabove(i,irt)=dddpfull(i,isl,irt)-dddpbelow(i,irt)
            tstarabove(i,irt)=tstarfull(i,isl,irt)-tstarbelow(i,irt)
          else
            tabove(i,irt)=tfull(i,isl,irt)
            dabove(i,irt)=dfull(i,isl,irt)
            dddpabove(i,irt)=dddpfull(i,isl,irt)
            tstarabove(i,irt)=tstarfull(i,isl,irt)
          endif
        enddo
        npabove(irt)=nps
        if(iprtlv.gt.0) then
          write(6,"(i3,2f10.2, 2i6)") irt,tabove(1,irt),tbelow(1,irt),
     #       npabove(irt),npbelow(irt)
        endif
      enddo
      if(iprtlv.gt.0) then
        write(6,"('bef',i3,2f10.2, 2i6)") 1,tabove(1,1),tbelow(1,1),
     #       npabove(1),npbelow(1)
      endif
      if(nrt.eq.2) then
        npabove(3)=npabove(1)
        npbelow(3)=npbelow(1)
        do i=1,npabove(3)
          tabove(i,3)=tabove(i,1)
          dabove(i,3)=dabove(i,1)
          dddpabove(i,3)=dddpabove(i,1)
          tstarabove(i,3)=tstarabove(i,1)
        enddo
        do i=1,npbelow(3)
          tbelow(i,3)=tbelow(i,1)
          dbelow(i,3)=dbelow(i,1)
          dddpbelow(i,3)=dddpbelow(i,1)
          tstarbelow(i,3)=tstarbelow(i,1)
        enddo
        if(iprtlv.gt.0) then
          write(6,"(i3,2f10.2, 2i6)") 3,tabove(1,3),tbelow(1,3),
     #       npabove(3),npbelow(3)
        endif
      endif
      return
      end
