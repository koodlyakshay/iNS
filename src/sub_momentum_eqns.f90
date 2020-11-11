subroutine momentum_eqns(nmiter, extiter, implicit_time, upwind, muscl)

use global_vars

implicit none

integer     :: i,j
integer     :: ExtIter, miter, nmiter
logical     :: implicit_time, upwind, muscl

!----------------------------------------------------------------------!
   !--- SU2 equivalent: CPBFludIteration()::Iterate()->SinglegridIteration(momentum)
   call preprocessing_flow
      
   !--- No set time step as we use fixed time step given as input ---!
   
   !--- Compute spatial discretization (Space Integration) ---!
   
   !--- Convective terms (1st order upwind/Central) ---!
   call convective_residual(upwind, muscl)
       
   !--- Viscous terms ---!
   call viscous_residual
   
   !--- Source term (pressure only) ---!
   call pressure_residual
   
   !--- Apply BC ---!
   !--- Right (i=Nx) ---!
   call right_boundary_flow
   
   !--- Left (i=1) ---!
   call left_boundary_flow
   
   !--- Lower (j=1) ---!
   call lower_boundary_flow
      
   !--- Upper (j=Ny) ---!
   call upper_boundary_flow
   
   if (implicit_time) then
     !--- Time Integration (Implicit) ---!
     call implicit_euler(miter, nmiter, extiter)
   else   
     !--- Time Integration (Explicit) ---!
     call explicit_euler(miter, nmiter, extiter)
   endif
   

end subroutine momentum_eqns


subroutine left_boundary_flow
!--- Left (i=1) ---!

use global_vars

implicit none

  integer   :: i, j, ipoint

 i=1
 do j=1,Ny
   iPoint = i + (j-1)*Nx
   !--- Specified velocity ---!
   U_old(1,iPoint) = 1.0
   U_old(2,iPoint) = 0.0
  
  !--- Update residual ---!
  R(1:2,iPoint) = 0.0
     
  Tot_Jac((iPoint-1)*nVar+1,:) = 0.0
  Tot_Jac((iPoint-1)*nVar+2,:) = 0.0
  Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 1.0
  Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 1.0
enddo   

end subroutine left_boundary_flow


subroutine right_boundary_flow
!--- Right (i=Nx) ---!
use global_vars

implicit none

  integer   :: i, j, ipoint

 !Zero gradient at the wall. Only a convective flux is added. 
 !The visc contribution is zero because of zero gradient.
   
 i=Nx
 do j=1,Ny
  iPoint = i + (j-1)*Nx
  !--- Velocity gradient is zero ---!
  !--- Have to add convective flux as is ---!     
  !--- Update residual ---!
  if ((j.eq.1) .or. (j.eq.Ny)) then
     R(1,iPoint) = R(1,iPoint) + max(rho*U_old(1,iPoint)*dy/2.0,0.0)*U_old(1,iPoint) ! max(rho*U_e*dy,0.0)*U_old(1,iPoint)
     R(2,iPoint) = R(2,iPoint) + max(rho*U_old(1,iPoint)*dy/2.0,0.0)*U_old(2,iPoint) ! max(rho*U_e*dy,0.0)*U_old(2,iPoint)        
  else
     R(1,iPoint) = R(1,iPoint) + max(rho*U_old(1,iPoint)*dy,0.0)*U_old(1,iPoint) ! max(rho*U_e*dy,0.0)*U_old(1,iPoint)
     R(2,iPoint) = R(2,iPoint) + max(rho*U_old(1,iPoint)*dy,0.0)*U_old(2,iPoint) ! max(rho*U_e*dy,0.0)*U_old(2,iPoint)      
  endif
     
  Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 1.0*U(1,iPoint)
  Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
  Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.5*U(2,iPoint)
  Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 0.5*U(1,iPoint)
 enddo

end subroutine right_boundary_flow


subroutine lower_boundary_flow
!--- Lower (j=1) ---!
use global_vars
implicit none

  integer   :: i, j, ipoint
 
 
 j=1
 do i=1,Nx
   iPoint = i + (j-1)*Nx
   if (x(i,j) .ge. 2.0e-2) then
     !--- Zero velocity ---!
     U_old(1:2,iPoint) = 0.0
     
     !--- Update residual ---!
     R(1:2,iPoint) = 0.0
     
     Tot_Jac((iPoint-1)*nVar+1,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+2,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 1.0
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 1.0
   else
   !--- Update residual ---!
     !Symmetry plane --> Add zero flux 
   endif
 enddo

end subroutine lower_boundary_flow


subroutine upper_boundary_flow
!--- Upper (j=Ny) ---!
use global_vars
implicit none

  integer   :: i, j, ipoint
 
 j=Ny
 do i=1,Nx
   iPoint = i + (j-1)*Nx
   if (x(i,j) .ge. 2.0e-2) then
    !--- Fixed wall ---!
     U_old(1:2,iPoint) = 0.0
     
     !--- Update residual ---!
     R(1:2,iPoint) = 0.0
     
     Tot_Jac((iPoint-1)*nVar+1,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+2,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 1.0
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 1.0
   else
     !--- Update residual ---!
     !Symmetry plane --> Add zero flux 
   endif
 enddo

end subroutine upper_boundary_flow
