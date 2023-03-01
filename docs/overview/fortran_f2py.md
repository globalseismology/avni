(fortran)=

# Fortran Guidelines


## Coding style

When modifying an existing file, try to maintain consistency with its original style.  If the code you add looks drastically different from the original code, it may be difficult for readers to follow. Try to avoid this. As a general guideline, we recommend the following code formatting style:

**give space for breathing:**

good
~~~fortran
  dx = 0.5 * fac * (a - b)
~~~

bad
~~~fortran
  dx=1/2*fac*(a-b)
~~~

Note that in performance critical sections, please use multiplication by 0.5 rather than divide by 2 for floating-points.

**use consistent 2-space indents:**

good
~~~fortran
  if (i == 1) then
    print *,'great'
  endif
~~~

bad
~~~fortran
  if(i == 1)then
        print *,'not so great'
  endif
~~~

**start your code with an indent:**

good
~~~fortran
  subroutine vbspl()
  implicit none
  ..
~~~

bad
~~~fortran
subroutine vbspl
  implicit none
  ..
~~~

The line beginning should only be used for *very important* sections, as it makes the line *very prominent* to read.
For example, only use it for function descriptions, important comments, or file headers. For comments, see also next point...

*exception, module definitions start at beginning:*

good
~~~fortran
module models
  integer :: count
end module
~~~

bad
~~~fortran
  module models
  integer :: count
  end module
~~~

**comment, comment, comment your code:**

good
~~~fortran
  ! gets associated values
  fg = vbspl(4,2)

  ! find values
  gt = fg
~~~

bad
~~~fortran
  fg = vbspl(4,2)
~~~

Note we prefer indenting the comments as well to make it easier for reading the code, e.g., when inside multiple if-then statements. Putting the comment at the beginning breaks the flow.

**comment, comment, comment your functions:**

good
~~~fortran
  subroutine blah_blah_function()

! calculates TI gradient based on a conjugate gradient method
!
! based on: Tarantola, inverse problem theory, 2005.
!                  section 6.22.7 conjugate directions, page 217.
!                  formula for alpha_n based on Polak & Ribiere (1969)
!
! note: we use a preconditioner F_0 = 1, thus lambda_n = gamma_n in (6.322)
!          and use gamma_n as the smoothed kernel (for bulk_c, bulk_betav,..).

  ..
~~~

bad
~~~fortran
  subroutine blah_blah_function()

! computation step

  ..
~~~

Note that we haven't been very strict in adopting a doxygen-readable function declaration.

**use double-colons for parameter declarations:**

good
~~~fortran
  integer :: i,j,k
~~~

bad
~~~fortran
  integer i,j,k
~~~

**use separators between subroutines:**

good
~~~fortran
  ..
  end subroutine

!
!----------------------------------------------------------
!

  subroutine get_color(icolor)
  ..
~~~

bad
~~~fortran
  ..
end subroutine

subroutine get_color(icolor)
  ..
~~~

## F2PY Troubleshooting

:::{warning}
AVNI used [F2PY](https://numpy.org/doc/stable/f2py/) to wrap legacy Python code by building from source code using the `numpy.distutils`. This requires restricting the versions to `numpy<=1.22` and `setuptools<60.0` in `setup.py` following [this link](https://numpy.org/doc/stable/reference/distutils_status_migration.html#moving-to-setuptools). Ultimately, AVNI will be migrated to use the latest version of the `setuptools` package.
:::

The purpose of the `F2PY` –Fortran to Python interface generator– utility is to provide a connection between Python and Fortran. `F2PY` is a part of NumPy (`numpy.f2py`) and also available as a standalone command line tool. We store legacy Fortran code in the `avni/f2py` folder. A way to troubleshoot `F2PY` issues involves the following steps:

1. Copy all of the source code to a single directory
2. Make sure that the code will actually compile without f2py (either with a makefile or with a simple script)
3. Use f2py to create a python signature file (`.pyf` file).

See the makefile attached below which will compile all of the code, and generate the file `avni_forward.pyf`. The `.pyf` file can be used to figure out what problems are going on. In our case, there were a number of subroutines that were not being properly translated to the signature file. These subroutines are called "unknown_subroutine" in the signature file. For example, the signature file complained about a subroutine in `rayseq.f`:

```{code-block} Fortran
---
caption: Fortran code block giving an unknown tag issue.
---
subroutine unknown_subroutine ! in rayseq.f
    integer :: iprtlv
    common /plevel/ iprtlv
end subroutine unknown_subroutine
```

Since it wasn't really giving us any reasons why the subroutine was unknown, it was a little hard to figure out how to fix. It turns out the problem (in this case) was that the subroutine definition in rayseq.f exceeded the fortran character limit of 72. So to fix it, we just added a continuation line to the source code:

```{code-block} Fortran
diff rayseq.f ../ALLCODES/rayseq.f
1c1,2
<       SUBROUTINE RAYSEQ (NP,NS,NN,ISS,IRR,ICNSx,ICNR,ICN,ISEQ,idirseg,NSEG)
---
>       SUBROUTINE RAYSEQ (NP,NS,NN,ISS,IRR,ICNSx,ICNR,ICN,ISEQ,
>      #idirseg,NSEG)
```

There were a number of other source files that we had to modify, but as far as we remember, the only changes we made were to fix character limit issues. Once there were no more problems in the python signature file, we copied the source codes back to their respective directories, and everything worked fine with `pip`. We can't say that this method will work for troubleshooting other `F2PY` issues.

```{code-block} Fortran
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
avni_forward.so: $(OBJS) avni_forward.pyf
    $(F2PY) -c avni_forward.pyf $(OBJS)
avni_forward.pyf: $(SRCS)
    $(F2PY) --overwrite-signature -m avni_forward -h avni_forward.pyf $(SRCS)
clean:
    $(RM) $(LOBJ)/*.o $(LOBJ)/*.pyf ./avni_forward.so
$(LOBJ)/%.o: %.f
    $(F77) $(FFLAGS) -c $*.f -o $(LOBJ)/$*.o
```
