      subroutine swap4of4(iword)
c
c---- swaps the order of 4 bytes of the 4-byte word
c
      byte iword(4)
      byte idummy
c
      idummy=iword(1)
      iword(1)=iword(4)
      iword(4)=idummy
c
      idummy=iword(2)
      iword(2)=iword(3)
      iword(3)=idummy
c
      return
      end
