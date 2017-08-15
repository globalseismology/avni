c     This subroutine converts polar co-ordinates to rectangular 
c     coordinates. 

c     Inputs   : Polar coordinates R and Theta
c     Outputs  : Rectangular coordinates X and Y

      subroutine pc2rc(r,theta,x,y)

         implicit none

         real r, theta, x, y

        intent(in) r, theta
        intent(out) x, y

         x = r*cos(theta)
         y = r*sin(theta)

      end subroutine pc2rc
      
      subroutine hwd()

          print*, 'Hello World !'

      end subroutine hwd
