      subroutine getptab(ianisoflag,ishflag,
     #           iss,isl,np,ns,id,pmin,pmax,nps,p,t,d,dddp,tstar,turn)
      common /plevel/ iprtlv
      dimension p(1)
      dimension t(1)
      dimension d(1)
      dimension dddp(1)
      dimension tstar(1)
      dimension turn(1)
c
      dimension np(1),ns(1),id(10,1)
c
      include 'rays.h'
c
      if(iprtlv.gt.1) then
        write(6,"('nlayers:',i4)") nlayers
        write(6,"('source index',i2)")iss
        write(6,"('multipliers:')")
        write(6,"(20i2)") (np(i),i=1,nlayers)
        write(6,"(20i2)") (ns(i),i=1,nlayers)
      endif
c
      ifirst=0
ccc      ilast=npfull(1,2)
      ilast=0
      do i=1,maxp
        t(i)=0.
        d(i)=0.
        dddp(i)=0.
        tstar(i)=0.
        turn(i)=0.
        if(pfull(i).gt.pmin.and.ifirst.eq.0) ifirst=i
        if(pfull(i).lt.pmax) ilast=i
      enddo
      nps=ilast-ifirst+1
      if(iprtlv.gt.0) then
      write(6,"('getptab: ifirst,ilast,nps',3i6)") ifirst,ilast,nps
      endif
c
      do i=ifirst,ilast
        p(i-ifirst+1)=pfull(i)
      enddo
c
      do it=1,2
        iarr=it
        if(it.eq.1.and.ishflag.ne.1) then
          iarr=iarr+2
        endif
        do il=1,nlayers
          xleg=0.
          if(it.eq.1) xleg=float(ns(nlayers+1-il))
          if(it.eq.2) xleg=float(np(nlayers+1-il))
          do i=ifirst,ilast
            if(xleg.gt.0.) then
              if(i.le.npfull(il,iarr)) then
                t(i-ifirst+1)=t(i-ifirst+1)+xleg*tfull(i,il,iarr)
                d(i-ifirst+1)=d(i-ifirst+1)+xleg*dfull(i,il,iarr)
                dddp(i-ifirst+1)=dddp(i-ifirst+1)+xleg*dddpfull(i,il,iarr)
                tstar(i-ifirst+1)=tstar(i-ifirst+1)+xleg*tstarfull(i,il,iarr)
                xxx=xturnfull(i,il,iarr)
              else
                xxx=-float(il)
              endif
              if(turn(i-ifirst+1).gt.0.5) then
              else
                if(xturnfull(i,il,iarr).gt.turn(i-ifirst+1)) then
                  if(xturnfull(i,il,iarr).gt.0.) then
                    turn(i-ifirst+1)=xturnfull(i,il,iarr) 
                  endif
                else if(xturnfull(i,il,iarr).lt.turn(i-ifirst+1)) then
                  turn(i-ifirst+1)=xturnfull(i,il,iarr)
                endif
              endif               
            endif
          enddo
        enddo
      enddo
c
      if(iprtlv.gt.1) then
        write(6,"('iss is:',i3)") iss
      endif
c
c-----p wave going up - subtract lower half
      if(iss.eq.5) then
          do i=ifirst,ilast
            if(i.le.npbelow(2)) then
              t(i-ifirst+1)=t(i-ifirst+1)-tbelow(i,2)
              tstar(i-ifirst+1)=tstar(i-ifirst+1)-tstarbelow(i,2)
              d(i-ifirst+1)=d(i-ifirst+1)-dbelow(i,2)
              dddp(i-ifirst+1)=dddp(i-ifirst+1)-dddpbelow(i,2)
            endif
          enddo
c-----p wave going down - subtract upper half
      else if(iss.eq.1) then
          do i=ifirst,ilast
            if(i.le.npabove(2)) then
              t(i-ifirst+1)=t(i-ifirst+1)-tabove(i,2)
              tstar(i-ifirst+1)=tstar(i-ifirst+1)-tstarabove(i,2)
              d(i-ifirst+1)=d(i-ifirst+1)-dabove(i,2)
              dddp(i-ifirst+1)=dddp(i-ifirst+1)-dddpabove(i,2)
            endif
          enddo
c-----s wave going up - subtract lower half
      else if(iss.eq.6) then
          iarr=1
          if(ishflag.ne.1) then
            iarr=3
          endif
          do i=ifirst,ilast
            if(i.le.npbelow(iarr)) then
              t(i-ifirst+1)=t(i-ifirst+1)-tbelow(i,iarr)
              tstar(i-ifirst+1)=tstar(i-ifirst+1)-tstarbelow(i,iarr)
              d(i-ifirst+1)=d(i-ifirst+1)-dbelow(i,iarr)
              dddp(i-ifirst+1)=dddp(i-ifirst+1)-dddpbelow(i,iarr)
            endif
          enddo
c-----s wave going down - subtract upper half
      else if(iss.eq.2) then
          iarr=1
          if(ishflag.ne.1) then
            iarr=3
          endif      
          do i=ifirst,ilast
            if(i.le.npabove(iarr)) then
              t(i-ifirst+1)=t(i-ifirst+1)-tabove(i,iarr)
              tstar(i-ifirst+1)=tstar(i-ifirst+1)-tstarabove(i,iarr)
              d(i-ifirst+1)=d(i-ifirst+1)-dabove(i,iarr)
              dddp(i-ifirst+1)=dddp(i-ifirst+1)-dddpabove(i,iarr)
            endif
          enddo
      else
        stop 'something is rotten in getray'
      endif
c
c---- now check the status of bottom reflections or turning points
c---- if xturn(p) is greater than 0, this wave has a turning point
c---- if turn(p) is less than 0, the value gives the reflecting layer
c
      do i=ifirst,ilast,100
c        write(6,"(i5,2f10.5)") i,p(i-ifirst+1),turn(i-ifirst+1)
      enddo
      do i=ifirst,ilast
        if(turn(i-ifirst+1).ge.0.) then
        else
          irefl=nlayers+1-nint(-turn(i-ifirst+1))
          if(id(1,irefl).gt.0.or.id(4,irefl).gt.0.or.id(7,irefl).gt.0) then
          else
c            write(6,"('not a turning ray:',i4,f10.5)") irefl,p(i-ifirst+1)
            d(i-ifirst+1)=0.
          endif
        endif
      enddo
      if(iprtlv.gt.0) then
        write(6,"(3f10.3)") t(1),d(1),p(1)
        write(6,"(3f10.3)") t(nps),d(nps),p(nps)
      endif
      return
      end
