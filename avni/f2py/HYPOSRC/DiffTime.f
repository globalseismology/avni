       subroutine  DiffTime(time1,time2,timediff)
c
c !OUTPUT PARAMETERS:
c
c                 Integer function returns number of seconds between the
c                 the times given as input.  -1 is returned in the event 
c                 of an error.  
c
c !DESCRIPTION:   This function returns the number of seconds between two
c                 times.  This function 
c                 determines the Julian day of each date using the "julday"
c                 function from the book "Numerical Recipes in FORTRAN, the 
c                 art of scientific computing (2nd Ed.), by William H. Press, 
c                 Saul A. Teukolsky, William T. Vetterling, and Brian P. 
c                 Flannery (Cambridge University Press, 1992).  This julian
c                 day is reduced by a constant and converted to seconds.  The
c                 reduction is required to allow the conversion to seconds to
c                 fit in a 32 bit integer.  The difference between the two 
c                 times is then calculated and returned.  The times need not
c                 be in chronological order as the function returns the abs
c                 value.  -1 is returned in the event of an error.
c
c !REVISION HISTORY:
c
c                 Modified after http://map.nasa.gov/GEOSGCMbrowser/html_code/GMAO_gfio/diffdate.f.html
c-------------------------------------------------------------------------

       integer StartDate, julday
       parameter (StartDate = 2433282)   
c      Use 01/01/1950 as base date
        
       dimension time1(6)
       dimension time2(6)	
       integer*4 year1,mon1,day1,year2,mon2,day2
       real*4 hour1,min1,sec1,hour2,min2,sec2
       integer*4 julian1,julian2
       real*4 julsec1, julsec2

       year1=int(time1(1))
       mon1=int(time1(2))
       day1=int(time1(3))
       year2=int(time2(1))
       mon2=int(time2(2))
       day2=int(time2(3))
       hour1=time1(4)
       min1=time1(5)
       sec1=time1(6)
       hour2=time2(4)
       min2=time2(5)
       sec2=time2(6)
       
c Get Julian Days and subtract off a constant (Julian days since 01/01/1950)
 
       julian1 = julday (mon1, day1, year1)
       julian1 = julian1 - StartDate
       julian2 = julday (mon2, day2, year2)
       julian2 = julian2 - StartDate
      
c Calculcate Julian seconds

       julsec1 = (real(julian1)-1)*86400 + hour1*3600 + min1*60 + sec1
       julsec2 = (real(julian2)-1)*86400 + hour2*3600 + min2*60 + sec2
       
       timediff = julsec2 - julsec1

       return
       end
