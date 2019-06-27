!> \file sub_preprocessing.f90
!! \brief Subroutine to perform some preprocessing operations.

!--- Set old variables, compute gradient ---!
subroutine preprocessing_flow
use global_vars

implicit none

integer           :: i,j
integer           :: iPoint, nPoint

  nPoint = Nx*Ny
  
   R = 0.0
   Tot_Jac = 0.0
   Fc = 0.0
   Fv = 0.0
   
   !--- Set old solution ---!
   
   do iPoint=1,nPoint
     U_old(1:2,iPoint) = U(1:2,iPoint)
   enddo   
   
   !--- Set Primitive ---!
   !--- 1 -> u, 2 -> v, 3 ->P ---!
   do i=1,Nx
    do j=1,Ny
    iPoint = i + (j-1)*Nx
      V(1:2,i,j) = U_old(1:2,iPoint)
      V(3,i,j) = P(i,j)
    enddo
   enddo
   
   !--- Compute gradients at nodes ---!
   call compute_gradient


end subroutine preprocessing_flow
