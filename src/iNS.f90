
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
!1. Add velocity correction array and write it to output file.
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

!--- SU2 equivalent: Run()->CPBFludIteration()::Iterate() ---!
!----------------------------------------------------------------------!
!---------------------- Solve momentum equation -----------------------!
!----------------------------------------------------------------------!
  do MIter = 1,nMIter
   !--- SU2 equivalent: CPBFludIteration()::Iterate()->SinglegridIteration(momentum)
   !--- Set old variables, compute gradient ---!
   call preprocessing
      
   !--- No set time step as we use fixed time step given as input ---!
   
   !--- Compute spatial discretization (Space Integration) ---!
   
   !--- Convective terms (1st order upwind/Central) ---!
   call convective_residual
       
   !--- Viscous terms ---!
    
   call viscous_residual
   
   !--- Source term (pressure only) ---!
   call pressure_residual
   
   !--- Apply BC ---!
   !--- Right outlet (i=Nx) ---!
   !Zero gradient at the wall. Only a convective flux is added. 
   !The visc contribution is zero because of zero gradient.
   
   i=Nx
   do j=1,Ny
   iPoint = i + (j-1)*Nx
    !--- Velocity gradient is zero ---!
    !--- Have to add convective flux as is ---!     
     !--- Update residual ---!
     if ((j.eq.1) .or. (j.eq.Ny)) then
        R(1,iPoint) = R(1,iPoint) + max(rho*U_old(1,iPoint)*dy/2.d0,0.d0)*U_old(1,iPoint) ! max(rho*U_e*dy,0.d0)*U_old(1,iPoint)
        R(2,iPoint) = R(2,iPoint) + max(rho*U_old(1,iPoint)*dy/2.d0,0.d0)*U_old(2,iPoint) ! max(rho*U_e*dy,0.d0)*U_old(2,iPoint)  
       
     else
        R(1,iPoint) = R(1,iPoint) + max(rho*U_old(1,iPoint)*dy,0.d0)*U_old(1,iPoint) ! max(rho*U_e*dy,0.d0)*U_old(1,iPoint)
        R(2,iPoint) = R(2,iPoint) + max(rho*U_old(1,iPoint)*dy,0.d0)*U_old(2,iPoint) ! max(rho*U_e*dy,0.d0)*U_old(2,iPoint)
        
     endif
     
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 1.0*U(1,iPoint)
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.5*U(2,iPoint)
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 0.5*U(1,iPoint)
   enddo
   !--- Left inlet (i=1) ---!
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


   !--- Lower wall (j=1) ---!
   j=1
   do i=1,Nx
   iPoint = i + (j-1)*Nx
   if (x(i,j).ge.10.d0) then
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
      
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=1,Nx
   iPoint = i + (j-1)*Nx
   if (x(i,j).ge.10.d0) then
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
   
   Res_l2 = 0.d0

   if (implicit_time) then
   !--- Time Integration (Implicit) ---!
   call implicit_euler
   
   else   
   !--- Time Integration (Explicit) ---!
   call explicit_euler
     
   endif
   
enddo   !MIter

!--- Compute residual mass flux from vel solution (SetPoissonSource()) ---!
call compute_massflux

!--- SU2 equivalent: CPBFludIteration()::Iterate()->SinglegridIteration(Poisson)
!----------------------------------------------------------------------!
!---------------- Solve Pressure correction equation now --------------!
!----------------------------------------------------------------------!
   Sol(3,:) = 0.d0
   P_Correc = 0.0
   GradPc = 0.0
do PIter = 1,nPIter

   R(3,:) = 0.d0  
   
     call compute_pcgradientgg
   
   
!--- Viscous Residual ---!
   !--- Assemble pressure equation ---!
   do i=1,Nx
    do j=1,Ny
      iPoint = i + (j-1)*Nx
      Jac(iPoint,:) = 0.0
      !--- Residual of the form ---!
      !--- R = -rho*\underbar{D}_f\nabla p'*area ---!
      !--- where \underbar{D}_f is the inverse of momentum coefficient of the two points calculated at the face ---!
      !East
      if (i .ne. Nx) then
      jPoint = i+1 + (j-1)*Nx
      if ((j.eq.1).or.(j.eq.Ny)) then
         R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*(GradPc(1,i,j) + GradPc(1,i+1,j)))*dy/2.0
      else 
         R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*(GradPc(1,i,j) + GradPc(1,i+1,j)))*dy
      endif
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dy/dx)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
       endif
      !West
      if (i .ne. 1) then
      jPoint = i-1 + (j-1)*Nx
      if ((j.eq.1).or.(j.eq.Ny)) then
          R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*(GradPc(1,i,j) + GradPc(1,i-1,j)))*dy/2.0
      else
          R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*(GradPc(1,i,j) + GradPc(1,i-1,j)))*dy
      endif
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dy/dx)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      endif
      !North
      if (j.ne.Ny) then
      jPoint = i + (j+1-1)*Nx
      if ((i.eq.1).or.(i.eq.Nx)) then
         R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*(GradPc(2,i,j) + GradPc(2,i,j+1)))*dx/2.0
      else
         R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*(GradPc(2,i,j) + GradPc(2,i,j+1)))*dx
      endif
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dx/dy)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      endif
      !South
      if (j.ne.1) then
      jPoint = i + (j-1-1)*Nx
      if ((i.eq.1).or.(i.eq.Nx)) then
         R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*(GradPc(2,i,j) + GradPc(2,i,j-1)))*dx/2.0
      else
         R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*(GradPc(2,i,j) + GradPc(2,i,j-1)))*dx
      endif
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dx/dy)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      endif
    enddo
   enddo
   
   !--- Boundary elements ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
      !South
      R(3,iPoint) = R(3,iPoint) + D(1,iPoint)*GradPc(2,i,j)*dx
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)

   enddo
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
     !North
     R(3,iPoint) = R(3,iPoint) - D(1,iPoint)*GradPc(2,i,j)*dx
     
     Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)

   enddo
   
   !--- Left inlet (i=1) ---!
   i=1
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
    !West
    R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*GradPc(1,i,j)*dy
    Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)

   enddo
   
   
   !--- Right outlet (i=Nx) ---!
   i=Nx
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    !East
    R(3,iPoint) = R(3,iPoint) - (D(1,iPoint))*GradPc(1,i,j)*dy
    Jac(iPoint,iPoint) = Jac(iPoint,iPoint) - Jac(iPoint,jPoint)

   enddo
   
   
   i=1
   j=1
   iPoint = i + (j-1)*Nx
 
   !Inlet face(West)
   R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*GradPc(1,i,j)*dy/2.0

   Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)  

   
   !South
   R(3,iPoint) = R(3,iPoint) + D(1,iPoint)*GradPc(2,i,j)*dx/2.0
   Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   

   i=1
   j=Ny
   iPoint = i + (j-1)*Nx

   !Inlet face(West)
   R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*GradPc(1,i,j)*dy/2.0

   Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
          
   !North
   R(3,iPoint) = R(3,iPoint) - D(1,iPoint)*GradPc(2,i,j)*dx/2.0
      
   Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      
   i=Nx
   j=Ny
   iPoint = i + (j-1)*Nx

   !East (outlet)
   R(3,iPoint) = R(3,iPoint) - (D(1,iPoint))*GradPc(1,i,j)*dy/2.0
   
   Jac(iPoint,iPoint) = Jac(iPoint,iPoint) - Jac(iPoint,jPoint)
      
   !North
   R(3,iPoint) = R(3,iPoint) - D(1,iPoint)*GradPc(2,i,j)*dx/2.0
   
   Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      
    
   i=Nx
   j=1
   iPoint = i + (j-1)*Nx

   !East
   R(3,iPoint) = R(3,iPoint) - (D(1,iPoint))*GradPc(1,i,j)*dy/2.0
   
   Jac(iPoint,iPoint) = Jac(iPoint,iPoint) - Jac(iPoint,jPoint)
      
   !South
   R(3,iPoint) = R(3,iPoint) + D(1,iPoint)*GradPc(2,i,j)*dx/2.0
         
   Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      
   if (implicit_time) then
   
   mass_l2 = 0.d0
   Res_l2 = 0.0
   do iPoint=1,nPoint
     R(3,iPoint) =  R(3,iPoint) + Mass(iPoint)
     R(3,iPoint) = -R(3,iPoint)
     
     Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Vol(iPoint)/dt_p
          
     Res_l2 = Res_l2 + R(3,iPoint)**2.0
   enddo
   close(19)

   Res_l2 = sqrt(Res_l2/nPoint)
   if ((modulo(PIter,p_screen2)).eq.0) print*,'Res(p): ',log10(Res_l2),PIter,ExtIter
      
   !--- Solve pressure correction equation ---!
   call seidel(0,nPoint,Jac(:,:),R(3,:),1.d0,Sol(3,:),Res(:),liniter,convergence)
   if (convergence .ne. 0) print*, 'Error in p',convergence,ExtIter
   
   do iPoint = 1,nPoint
     P_Correc(iPoint) = P_Correc(iPoint) + Sol(3,iPoint)
   enddo
   
   !--- Time Integration (Explicit) ---!
   else  
     Res_l2 = 0.d0
     do iPoint=1,nPoint
        R(3,iPoint) =  R(3,iPoint) + Mass(iPoint)
        P_Correc(iPoint) = P_Correc(iPoint) - R(3,iPoint)*dt_p/Vol(iPoint)
        Res_l2 = Res_l2 + R(3,iPoint)**2.0
     enddo
     
     Res_l2 = sqrt(Res_l2/nPoint)
   if ((modulo(PIter,p_screen2)).eq.0) print*,'Res(p): ',log10(Res_l2),PIter,ExtIter
   endif 
   
   if (Piter .eq. nPIter)call compute_pcgradientgg
    
    
 enddo   !PIter
  
 !--- Correct_Velocity(), Correct_Pressure() ---!
 !--- Correct pressure and velocities ---!
  
 call sub_flow_correction
  

   !--- Convergence monitoring ---!
   
   !--- Output solution ---!
   if (modulo(ExtIter,p_out) .eq. 0) then
     call sub_output
   endif
enddo !ExtIter

call sub_output

close(11)
close(12)

call delete_vars

end program iNS







