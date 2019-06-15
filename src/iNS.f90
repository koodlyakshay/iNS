
!> \file iNS.f90
!! \brief Program to solve channel flow using SIMPLE
!!		Implementation based on SU2 
!!		Governing Eqn:
!       \f$ \partial_t U + \nabla F_c - \nabla F_v = Q \f$
!!		U -> Solution vector, F_c -> Convective flux, F_v -> Viscous flux, Q -> Source term
!!		Discretized Eqn:
!!		(\frac{\Omega_i}{\Delat t^n} + \frac{\partial R_i(U^n)}{\partial U_j}).\Delata U_j^n = -R_i(U^n)
!!		Residual, R_i(U) = \Sigma_{j\in N(i)} (F_{c_{ij}} - F_{v_{ij}})\Delta S_{ij} - Q_i|\Omega_i|
!!		Jacobian, J_{ij} = \frac{\partial R_i(U^n)}{\partial U_j}
!! \author Akshay KR
!! \version 0.1 

!---------------------------------------------------------------------------------------
!To do list
!2. Make subroutine for boundary conditions.
!3. Make seperate subroutines for momentum and pressure correction eqn, can remain in this file only.
!4. Add some logicals to clarify the physics, ex - viscous(to check  if euler/NS).
!5. Change sizes in allocate vars to accommodate 3d if necessary.
!6. Create functions like SetTotJacZero, AddJac, SubtractJac etc and pass the necessary 
!   value and indicies to the function instead of existing code.
!---------------------------------------------------------------------------------------

program iNS

use global_vars

implicit none

logical           :: old

inquire(file = '../out/',exist=old)

if (.not. old) then
  call system('mkdir ../out/')
endif  

call read_input

call allocate_vars

call initialize_vars

!--- Begin Solver ---!

do ExtIter = 1,nExtIter

!---------------------- Solve momentum equation -----------------------!

call momentum_eqns

!--- Compute residual mass flux from vel solution ---!
call compute_massflux

!---------------- Solve pressure corrrection equation -----------------!

call pressure_correc_eqn

!--- Correct pressure and velocities ---!
  
call sub_flow_correction
  

!--- Convergence monitoring ---!
   

!--- Output solution ---!

if (modulo(ExtIter,p_out) .eq. 0) then
     call sub_output
endif


enddo !ExtIter

!--- Call output at the end of loop regardless of nExtIter ---!
call sub_output

!--- Deallocate variables ---!
call delete_vars

end program iNS







