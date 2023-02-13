      integer*4 function julday(iy,imo,idm)
      integer*4 mon(12)
      data mon/0,31,59,90,120,151,181,212,243,273,304,334/
      julday=idm+mon(imo)
      if(imo.gt.2) julday=julday+lpyr(iy)
      return
      end

