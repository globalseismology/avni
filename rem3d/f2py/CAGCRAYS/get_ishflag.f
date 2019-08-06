      subroutine get_ishflag(phase,ishflag)
      character*16 phase
c     returns whether athe phase is SH (ishflag=1) or PSV  
c     SKS and SKS are on SV as SH cannot interact with core.
c     Can make a significant difference (>~2s) for anisotropic models
      if(phase(1:lnblnk(phase)).eq.'S'.or.phase(1:lnblnk(phase)).eq.'SS'.or.
     #   phase(1:lnblnk(phase)).eq.'ScS') then
			ishflag=1
			write(6,"('Note: The phase ',a,' is a SH wave defined in get_ishflag. Using ishflag=1 in cagcrays_pm')") 
     #          phase(1:lnblnk(phase))
	  else if(phase(1:lnblnk(phase)).eq.'P'.or.phase(1:lnblnk(phase)).eq.'PP'.or.
     #   phase(1:lnblnk(phase)).eq.'PPP'.or.
     #   phase(1:lnblnk(phase)).eq.'PcP'.or.phase(1:lnblnk(phase)).eq.'ScP'.or.
     #   phase(1:lnblnk(phase)).eq.'SP'.or.phase(1:lnblnk(phase)).eq.'PKP'.or.
     #   phase(1:lnblnk(phase)).eq.'PKIKP'.or.phase(1:lnblnk(phase)).eq.'PKKP'.or.
     #   phase(1:lnblnk(phase)).eq.'SKIKP'.or.phase(1:lnblnk(phase)).eq.'SKP'.or.
     #   phase(1:lnblnk(phase)).eq.'SKS'.or.phase(1:lnblnk(phase)).eq.'SKKS'.or.
     #   phase(1:lnblnk(phase)).eq.'PKKP'.or.phase(1:lnblnk(phase)).eq.'PKPPKP') then
			
			ishflag=0
	  else
	  	write(6,"('Error: The phase ',a,' not specified in get_branch')") phase(1:lnblnk(phase))
	  endif
					
	  return 
	  end			

