c      parameter (maxp=14000)
c      parameter (maxp=20500)
      parameter (maxp=16000)
      parameter (maxl=25)
      common /ptablex/ nlayers,pfull(maxp),
     #        tfull(maxp,maxl,3),dfull(maxp,maxl,3),dddpfull(maxp,maxl,3),
     #        tstarfull(maxp,maxl,3),xturnfull(maxp,maxl,3),
     #        npfull(maxl,3),topfull(maxl),botfull(maxl),
     #        tabove(maxp,3),tbelow(maxp,3),
     #        dabove(maxp,3),dbelow(maxp,3),
     #        dddpabove(maxp,3),dddpbelow(maxp,3),
     #        tstarabove(maxp,3),tstarbelow(maxp,3),
     #        npabove(3),npbelow(3)
