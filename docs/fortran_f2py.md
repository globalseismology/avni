F2PY troubleshooting
--------------------
A way to troubleshoot f2py issues trouble shooting was to i) copy all of the source code to a single directory ii) make sure that the code will actually compile without f2py (either with a makefile or with a simple script) and iii) Use f2py to create a python signature file (.pyf file). See the makefile attached below which will compile all of the code, and generate the file *rem3d_forward.pyf*. The .pyf file can be used to figure out what problems are going on. In our case, there were a number of subroutines that were not being properly translated to the signature file. These subroutines are called "unknown_subroutine" in the signature file. For example, the signature file complained about a subroutine in rayseq.f:

~~~fortran
subroutine unknown_subroutine ! in rayseq.f
    integer :: iprtlv
    common /plevel/ iprtlv
end subroutine unknown_subroutine
~~~
Since it wasn't really giving me any reasons why the subroutine was unknown, it was a little hard to figure out how to fix. It turns out the problem (in this case) was that the subroutine definition in rayseq.f exceeded the fortran character limit of 72. So to fix it, I just added a continuation line to the source code:
~~~
diff rayseq.f ../ALLCODES/rayseq.f
1c1,2
<       SUBROUTINE RAYSEQ (NP,NS,NN,ISS,IRR,ICNSx,ICNR,ICN,ISEQ,idirseg,NSEG)
---
>       SUBROUTINE RAYSEQ (NP,NS,NN,ISS,IRR,ICNSx,ICNR,ICN,ISEQ,
>      #idirseg,NSEG)
~~~
There were a number of other source files that I had to modify, but as far as I remember, the only changes I made were to fix character limit issues. Once there were no more problems in the python signature file, I copied the source codes back to their respective directories, and everything worked fine with pip. I can't say that this method will work for troubleshooting other f2py issues, but it worked for me.