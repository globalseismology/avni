      subroutine pdaz(xla,xlo,az,del,yla,ylo)
Cf2py intent(inout) xla,xlo,del,az
Cf2py intent(out) yla,ylo
      double precision xla8,xlo8,del8,az8,yla8,ylo8
c
      xla8=dble(xla)
      xlo8=dble(xlo)
      del8=dble(del)
      az8=dble(az)
c
      call dpdaz(xla8,xlo8,az8,del8,yla8,ylo8)
c
      yla=sngl(yla8)
      ylo=sngl(ylo8)
c
      return
      end
