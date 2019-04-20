Table of Contents
-----------------

  1. [Use *issues* to discuss intended modifications](#use-issues-to-discuss-any-intended-modifications)
  2. [Document your code](#document-your-code)
  3. [Coding style](#coding-style)
  	* [Python](#python-formatting)
  	* [Fortran](#fortran-formatting)

Use *issues* to discuss intended modifications
--------------------------------------------------

GitHub provides a [system](https://github.com/globalseismology/rem3d/issues) to track issues. It should be a central place to monitor REM3D evolution. In particular:

-   report bug as they occur

-   plan modifications.

REM3D's issue tracker interface lets us track bugs being fixed and enhancements being made by assigning labels. We will reserve the labels

-   *bug* for issue fixing

-   *enhancement* for feature development.


Document your code
------------------

Any new code should be fully Doxygen commented in [Python](#python-formatting) or [fortran](#fortran-formatting). If you have some free time, feel free to comment any code you modify.


Coding style
------------

When modifying an existing file, try to maintain consistency with its original style.  If the code you add looks drastically different from the original code, it may be difficult for readers to follow. Try to avoid this. As a general guideline, we recommend the following code formatting style:

  - [Python](#python-formatting)
  - [Fortran](#fortran-formatting)

Python formatting
------------------

**give space for breathing, use 4 spaces instead of tabs:**

good
~~~python
    dx = 0.5 * fac * (a - b)
~~~

bad
~~~python
	dx=1/2*fac*(a-b)
~~~

**comment, comment, comment your functions using Doxygen convention:**

good
~~~python
"""@package docstring
Documentation for this module.
More details.
"""
def func():
    """Documentation for a function.
    More details.
    """
    pass
class PyClass:
    """Documentation for a class.
    More details.
    """
   
    def __init__(self):
        """The constructor."""
        self._memVar = 0;
        pass
~~~


Fortran formatting
------------------

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
  ! gets associated normal on GLL point
  ! (note convention: pointing outwards of acoustic element)
  nx = coupling_ac_el_normal(1,igll,iface)

  ! continuity of displacement and pressure on global point
  accel(1,iglob) = accel(1,iglob) + jacobianw*nx*pressure  
~~~

bad
~~~fortran
  nx = coupling_ac_el_normal(1,igll,iface)
  accel(1,iglob) = accel(1,iglob) + jacobianw*nx*pressure
~~~

Note we prefer indenting the comments as well to make it easier for reading the code, e.g., when inside multiple if-then statements. Putting the comment at the beginning breaks the flow.

**comment, comment, comment your functions:**

good
~~~fortran
  subroutine get_cg_direction_tiso()

! calculates TI gradient based on a conjugate gradient method
!
! based on: Tarantola, inverse problem theory, 2005.
!                  section 6.22.7 conjugate directions, page 217.
!                  formula for alpha_n based on Polak & Ribiere (1969)
!
! note: we use a preconditioner F_0 = 1, thus lambda_n = gamma_n in (6.322)
!          and use gamma_n as the smoothed kernel (for bulk_c, bulk_betav,..).
!
!          however, one could see smoothing as preconditioner F_0, thus
!          gamma_n would be un-smoothed kernel and lambda_n would be smoothed one...
!          i'm not sure if this makes a difference.

  ..
~~~

bad
~~~fortran
  subroutine get_cg_direction_tiso()

! CG step

  ..
~~~

Note that we haven't been very strict in adopting a doxygen-readable function declaration. It is nice though to have it.
For example, you would do:
~~~fortran
!> Define an ADIOS scalar integer variable and autoincrement the adios
!! group size by (4).
!! \param adios_group The adios group where the variables belongs
!! \param group_size_inc The inout adios group size to increment
!!                       with the size of the variable
!! \param path The logical path structuring the data and containing
!!             the variable
!! \param name The variable name in the ADIOS file.
!! \param var The variable to be defined. Used for type inference. Can be
!             ignored.
!!
!! \note See define_adios_double_scalar()
  subroutine define_adios_integer_scalar(adios_group, group_size_inc, path, name, var)
  ..
~~~


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

