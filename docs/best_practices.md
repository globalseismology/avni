Table of Contents
-----------------

  1. [Use *issues* to discuss intended modifications](#use-issues-to-discuss-any-intended-modifications)
  2. [Document your code](#document-your-code)
  3. [Test your code](#test-your-code)
  4. [F2PY troubleshooting](#f2py-troubleshooting)
  6. [Coding style](#coding-style)
  	* [Python](#python-formatting)
  	* [Fortran](#fortran-formatting)

Use *issues* to discuss intended modifications
----------------------------------------------

GitHub provides a [system](https://github.com/geodynamics/avni/issues) to track issues. It should be a central place to monitor AVNI evolution. In particular:

-   report bug as they occur

-   plan modifications.

AVNI's issue tracker interface lets us track bugs being fixed and enhancements being made by assigning labels. We will reserve the labels

-   *bug* for issue fixing

-   *enhancement* for feature development.


Document your code
------------------

Any new code should be fully Doxygen commented in [Python](#python-formatting), [MATLAB](#matlab-formatting) or [fortran](#fortran-formatting). If you have some free time, feel free to comment any code you modify.

F2PY troubleshooting
------------------

F2PY is used by setup.py to compile fortran routines during installation. As a developer, you might need to troubleshoot the fortran codes to make them compile correctly. The error handling and exceptions are not straightforward; a few suggestions are provided [here](./fortran_f2py.md)


Test your code
------------------

[![Build Status](https://travis-ci.com/globalseismology/avni.svg?token=Z1JjFn7SrxG1nGGE9y1u&branch=devel)](https://travis-ci.com/globalseismology/avni) [![codecov](https://codecov.io/gh/globalseismology/avni/branch/devel/graph/badge.svg?token=NTCVjCUfJm)](https://codecov.io/gh/globalseismology/avni)

Any new code should be tested with unit tests kept in the [tests](../tests) folder. Note that we use [coverage.py](https://coverage.readthedocs.io) for testing so all files and routines should be named test_*. We attempt to keep coverage above 90% on our development builds.

Coding style
------------

When modifying an existing file, try to maintain consistency with its original style.  If the code you add looks drastically different from the original code, it may be difficult for readers to follow. Try to avoid this. As a general guideline, we recommend the following code formatting style:

  - [Python](#python-formatting)
  - [Fortran](#fortran-formatting)
  - [MATLAB](#matlab-formatting)

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

MATLAB formatting
------------------

**comment, comment, comment your functions using Doxygen convention:**

good
~~~matlab
function out=function(input)
% cm_data=rem3dcolorscal(m, reverse)
% Loads the color scale
%
% INPUT:
%
% m            number of contours (default: same as in cm below : 512)
% reverse      0 for standard, 1 for reversing color palatte (default: 0)
%
% OUTPUT:
%
% d            interpolated values for use with colormap
%
% EXAMPLE:
%
% colormap(rem3dcolorscal(20));
%
% TESTED ON:
%
% R2019a (9.6.0.1072779)
~~~

bad
~~~matlab
function out=function(input)
~~~

**give space for breathing:**

good
~~~matlab
  dx = 0.5 * fac * (a - b);
~~~

bad
~~~matlab
  dx=1/2*fac*(a-b);
~~~

Note that in performance critical sections, please use multiplication by 0.5 rather than divide by 2 for floating-points.

**provide the default value of arguments:**

good
~~~matlab
function out=function(input)

% Comment on what input does
defval('input',0)
~~~

**use optional output arguments:**

good
~~~matlab
function varargout=function(input)

% optional output
varns = {out};
varargout = varns(1:nargout);

if nargout == 0:
   %plot something
end
~~~

bad
~~~matlab
function out=function(input)
~~~




