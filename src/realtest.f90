!! A code to test precision in gfortran
!!
!! run this code with:
!! gfortran -fdefault-real-8 realtest.f90 -o realtest
!! 
!! Simply: functions like sin, sqrt, abs takes the precision of the argument.
!!          if the argumrnt is real (like 27.0) then the answer is real
!!          if the argumrnt is double (like 27.0d0) then the answer is double
!!          The output can be controlled via compiler flags like "-fdefault-real-8".
!!          I think, the best practice is to use d0 by default 
!!          which gives more control over the result.
program real_test
implicit none

real      :: r, d
real*4    :: r4
real*8    :: r8
double precision :: d8


  r  = sin(27.0)
  d = sin(27.0d0)
  r4 = sin(27.0)
  r8 = sin(27.0d0)
  d8 = dsin(27.0d0) 

print*, r, d, r4, r8, d8



end program real_test