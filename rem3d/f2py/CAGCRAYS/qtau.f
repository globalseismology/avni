      subroutine qtau(ilay,x,pray,isotropic,iraytype,qray)
      implicit double precision (a-h,o-z)
      logical isotropic
c
      call qtauall(ilay,x,pray,isotropic,iraytype,qray,
     #          vpv,vph,vsv,vsh,eta,qmu,qkappa,s1,s2,s3,s4,s5,r)
      return
      end
