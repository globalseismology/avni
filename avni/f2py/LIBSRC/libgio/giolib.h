      parameter (mxfiles=20)
      integer isize(mxfiles)
      integer ipos(mxfiles)
      integer irecl(mxfiles)
      integer irw(mxfiles)
      integer iopen(mxfiles)
      common /gioheader/ iopen,
     #       ipos,
     #       irw,
     #       isize,
     #       irecl
      save /gioheader/
