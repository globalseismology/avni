      subroutine wgint(ilay,r1,r2,xnorm,iturn,fxtpp,dxtdp,pray,isotropic
     #               ,iraytype,nint,deltainc,ttimeinc,dddpinc,tstarinc)
c action: performs gauss-legendre integrations
c input parameters (supplied by routine integr):
c     r1,r2  - integration bounds
c     nin    - (idleg in calling routine) order of the legendre polynomial
c     f(external)- fqanis (integral kernel)
c     iq     - number of the layer
c     fint   - array, contains the values (delta,t,d(delta)/d(p) and
c              partials derivatives) requested by the user (qints(3)=
c              fint(3) will undergo a later scaling in routine integr)
c     nint   - number of elements of fint requested by the user
c  ***** nint must be less than nmax *****
c
c the choice of the order of the legendre integration is limited to
c a few values (recorded via a data statement in array nord). through
c nnord and nord the nearest higher selected order is choosen. maximum
c order is 24. array x contains the zeros and w contains the weights.
c array iad contains the addresses of the first zero and weight for
c each polynomial order. the symmetry of zeros and weights is welcomed.
c
c    the integration bounds are r1,r2. they are transformed to -1,1
c with the new variable t: x=((r1+r2)/2)+t*((r2-r1)/2). the integral
c is computed in t and then multiplied by (r2-r1)/2.
c
      implicit double precision(a-h,o-z)
      parameter (nmax=4)
      dimension fint(nmax),vals(nmax),vals1(nmax),sum(nmax)
      dimension iad(13),nord(24),x(65),w(65),w1(29),w2(36)
     1     ,x1(29),x2(36),nnord(13)
      equivalence (x(1),x1(1)),(x(30),x2(1)),(w(1),w1(1))
     1     ,(w(30),w2(1))
      data  iad/ 1, 2, 4, 6, 9,12,16,20,25,30,36,44,54/
      data nord/1,1,2,3,4,5,6,7,8,9,10,10,11,11,11,11
     1     ,12,12,12,12,13,13,13,13/
      data nnord/ 2, 3, 4, 5, 6, 7, 8, 9,10,12,16,20,24/
      data x1/ .57735 02691 89626d0,.00000 00000 00000d0
     1        ,.77459 66692 41483d0,.33998 10435 84856d0
     2        ,.86113 63115 94053d0,.00000 00000 00000d0
     3        ,.53846 93101 05683d0,.90617 98459 38664d0
     4        ,.23861 91860 83197d0,.66120 93864 66265d0
     5        ,.93246 95142 03152d0,.00000 00000 00000d0
     6        ,.40584 51513 77397d0,.74153 11855 99394d0
     7        ,.94910 79123 42759d0,.18343 46424 95650d0
     8        ,.52553 24099 16329d0,.79666 64774 13627d0
     9        ,.96028 98564 97536d0,.00000 00000 00000d0
     x        ,.32425 34234 03809d0
     1        ,.61337 14327 00590d0,.83603 11073 26636d0
     2        ,.96816 02395 07626d0,.14887 43389 81631d0
     3        ,.43339 53941 29247d0,.67940 95682 99024d0
     4        ,.86506 33666 88985d0,.97390 65285 17172d0/
      data x2/ .12523 34085 11469d0,.36783 14989 98180d0
     1        ,.58731 79542 86617d0,.76990 26741 94305d0
     2        ,.90411 72563 70475d0,.98156 06342 46719d0
     3        ,.09501 25098 37637d0,.28160 35507 79259d0
     4        ,.45801 67776 57227d0,.61787 62444 02644d0
     5        ,.75540 44083 55003d0,.86563 12023 87831d0
     6        ,.94457 50230 73233d0,.98940 09349 91650d0
     7        ,.07652 65211 33497d0,.22778 58511 41645d0
     8        ,.37370 60887 15420d0,.51086 70019 50827d0
     9        ,.63605 36807 26515d0,.74633 19064 60151d0
     x        ,.83911 69718 22218d0,.91223 44282 51325d0
     1        ,.96397 19272 77913d0,.99312 85991 85094d0
     2        ,.06405 68928 62605d0,.19111 88674 73616d0
     3        ,.31504 26796 96163d0,.43379 35076 26045d0
     4        ,.54542 14713 88840d0,.64809 36519 36976d0
     5        ,.74012 41915 78554d0,.82000 19859 73902d0
     6        ,.88641 55270 04401d0,.93827 45520 02733d0
     7        ,.97472 85559 71309d0,.99518 72199 97021d0/
      data w1/1.00000 00000 00000d0,.88888 88888 88889d0
     1        ,.55555 55555 55556d0,.65214 51548 62546d0
     2        ,.34785 48451 37454d0,.56888 88888 88889d0
     3        ,.47862 86704 99366d0,.23692 68850 56189d0
     4        ,.46791 39345 72691d0,.36076 15730 48139d0
     5        ,.17132 44923 79170d0,.41795 91836 73469d0
     6        ,.38183 00505 05119d0,.27970 53914 89277d0
     7        ,.12948 49661 68870d0,.36268 37833 78362d0
     8        ,.31370 66458 77887d0,.22238 10344 53374d0
     9        ,.10122 85362 90376d0,.33023 93550 01260d0
     x        ,.31234 70770 40003d0
     1        ,.26061 06964 02935d0,.18064 81606 94857d0
     2        ,.08127 43883 61574d0,.29552 42247 14753d0
     3        ,.26926 67193 09996d0,.21908 63625 15982d0
     4        ,.14945 13491 50581d0,.06667 13443 08688d0/
      data w2/ .24914 70458 13403d0,.23349 25365 38355d0
     1        ,.20316 74267 23066d0,.16007 83285 43346d0
     2        ,.10693 93259 95318d0,.04717 53363 86512d0
     3        ,.18945 06104 55068d0,.18260 34150 44923d0
     4        ,.16915 65193 95003d0,.14959 59888 16576d0
     5        ,.12462 89712 55534d0,.09515 85116 82493d0
     6        ,.06225 35239 38648d0,.02715 24594 11754d0
     7        ,.15275 33871 30726d0,.14917 29864 72604d0
     8        ,.14209 61093 18382d0,.13168 86384 49176d0
     9        ,.11819 45319 61518d0,.10193 01198 17240d0
     x        ,.08327 67415 76705d0,.06267 20483 34109d0
     1        ,.04060 14298 00387d0,.01761 40071 39152d0
     2        ,.12793 81953 46752d0,.12583 74563 46828d0
     3        ,.12167 04729 27803d0,.11550 56680 53726d0
     4        ,.10744 42701 15966d0,.09761 86521 04114d0
     5        ,.08619 01615 31953d0,.07334 64814 11080d0
     6        ,.05929 85849 15437d0,.04427 74388 17420d0
     7        ,.02853 13886 28934d0,.01234 12297 99987d0/
c
      if(nint.gt.nmax) then
         write(6,"('nmax in routine gauslv needs increasing',2i6)") nint,nmax
         stop
      endif
c
c---- requested degree of integration
c
c---- nin is the requested degree of the integration
c---- the actual order is nnord(nord(n)); the zeros are from
c---- x(iad(nord(n)) on, the weights start at w(iad(nord(n)))=w(ia)
c
      nin=6
      n=nin
      if(nin.lt.1) then
        n=1
      endif
      if(nin.gt.24) then
        n=24
      endif
      ind=nord(n)
      n=nnord(ind)
      ia=iad(ind)
      nc=n/2
      nc2=2*nc
c
      y1=.5d0*(r2+r1)
      y2=.5d0*(r2-r1)
      do j=1,nint
        sum(j)=0.d0
      enddo
c
      if(n.ne.nc2) then
        call fqs(ilay,y1,xnorm,iturn,fxtpp,dxtdp,pray,isotropic,iraytype,nint,vals)
        do j=1,nint
          sum(j)=w(ia)*vals(j)
        enddo
        ia=1+ia
      end if
c
      if(nc.ne.0) then
        do i=1,nc
          t1=x(ia)*y2
          call fqs(ilay,y1+t1,xnorm,iturn,fxtpp,dxtdp,pray,isotropic,iraytype,nint,vals)
          call fqs(ilay,y1-t1,xnorm,iturn,fxtpp,dxtdp,pray,isotropic,iraytype,nint,vals1)
c
          do j=1,nint
            sum(j)=sum(j)+w(ia)*(vals(j)+vals1(j))
          enddo
          ia=1+ia
        enddo
      endif
c
      do j=1,nint
        fint(j)=y2*sum(j)
      enddo
c
      if(nint.ge.1) then
        deltainc=fint(1)
      endif
      if(nint.ge.2) then
        ttimeinc=fint(2)
      endif
      if(nint.ge.3) then
        dddpinc=fint(3)
      endif
      if(nint.ge.4) then
        tstarinc=fint(4)
      endif
c      
      return
      end
