!> \file mod_parameters.f90
!! \brief File with iteration parameters.
module parameters

implicit none

!integer           :: i,j,k,iPoint,jpoint

!> @param nExtIter  number of external iterations
integer           :: nExtIter
!> @param nPIter  number of iterations for Poisson solver
integer           :: nPIter
!> @param nMIter  number of iterations for momentum equations
integer           :: nMIter

logical           :: implicit_time
logical           :: upwind
logical           :: muscl
logical           :: restart


end module parameters
