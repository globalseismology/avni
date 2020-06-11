      subroutine evemell(rad,ell,eta2,ierr)
      implicit double precision (a-h,o-z)
      include 'emcommon.h'
      ell=drsple(1,numemlev,radlev,elllev,ellspl,rad)
      eta2=drsple(1,numemlev,radlev,eta2lev,eta2spl,rad)
      ierr=0
      return
      end
