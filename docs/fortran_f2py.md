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

~~~fortran
# makefile
#LOBJ = ./objects
LOBJ = ./
FFLAGS = -ffixed-line-length-none -fPIC -O3 -fbounds-check
SRCS = tt_predict.f acosd.f asind.f atand.f bffi.f bullen.f cagcrays_pm.f closfl.f conv2geocen.f cosd.f dacosd.f dbsplrem.f dcosd.f ddelaz.f ddelazgc.f \
    delaz.f delazgc.f DiffTime.f drspleder.f drspledr.f drsple.f drspln.f dsind.f dxt.f evemb.f evemdr.f evemell.f evem.f fdrays.f fqs.f \
    get_branch.f getcmt5.f getcmtbyname.f getcmtdatetime.f getflpos.f get_ishflag.f getptab.f julday.f legndr.f lpyr.f monday.f openfl.f \
    pdaz.f prange.f psvrayin.f qtauall.f qtau.f qtauzero.f rayseq.f reademb.f reademfl.f readnbn.f readptab.f sind.f swap4of4.f tadder.f \
    tand.f vbspl.f wgint.f wgray.f ylm.f
OBJS = tt_predict.o acosd.o asind.o atand.o bffi.o bullen.o cagcrays_pm.o closfl.o conv2geocen.o cosd.o dacosd.o dbsplrem.o dcosd.o ddelaz.o ddelazgc.o \
    delaz.o delazgc.o DiffTime.o drspleder.o drspledr.o drsple.o drspln.o dsind.o dxt.o evemb.o evemdr.o evemell.o evem.o fdrays.o fqs.o \
    get_branch.o getcmt5.o getcmtbyname.o getcmtdatetime.o getflpos.o get_ishflag.o getptab.o julday.o legndr.o lpyr.o monday.o openfl.o \
    pdaz.o prange.o psvrayin.o qtauall.o qtau.o qtauzero.o rayseq.o reademb.o reademfl.o readnbn.o readptab.o sind.o swap4of4.o tadder.o \
    tand.o vbspl.o wgint.o wgray.o ylm.o
F2PY = f2py
F77 = gfortran
F90 = fortran
CC = gcc
rem3d_forward.so: $(OBJS) rem3d_forward.pyf
    $(F2PY) -c rem3d_forward.pyf $(OBJS)
rem3d_forward.pyf: $(SRCS)
    $(F2PY) --overwrite-signature -m rem3d_forward -h rem3d_forward.pyf $(SRCS)
clean:
    $(RM) $(LOBJ)/*.o $(LOBJ)/*.pyf ./rem3d_forward.so
$(LOBJ)/%.o: %.f
    $(F77) $(FFLAGS) -c $*.f -o $(LOBJ)/$*.o
~~~
