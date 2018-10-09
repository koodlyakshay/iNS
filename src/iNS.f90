!Program to solve channel flow using SIMPLE 
!Implementation based on SU2 
!Governing Eqn:
!             \partial_t U + \nabla F_c - \nabla F_v = Q
!U -> Solution vector, F_c -> Convective flux, F_v -> Viscous flux, Q -> Source term
!Discretized Eqn:
!(\frac{\Omega_i}{\Delat t^n} + \frac{\partial R_i(U^n)}{\partial U_j}).\Delata U_j^n = -R_i(U^n)
!Residual, R_i(U) = \Sigma_{j\in N(i)} (F_{c_{ij}} - F_{v_{ij}})\Delta S_{ij} - Q_i|\Omega_i|
!Jacobian, J_{ij} = \frac{\partial R_i(U^n)}{\partial U_j}
program iNS

use global_vars

implicit none

call read_input

call allocate_vars

call initialize_vars

!--- Output Files ---!
if (wrt_data .eq. 1) then
  open(unit=16,file='../out/Convbudgetspb.txt',status='unknown')
  open(unit=17,file='../out/Viscbudgetspb.txt',status='unknown')
  open(unit=18,file='../out/presbudgetspb.txt',status='unknown')
endif

!--- Begin Solver ---!

do ExtIter = 1,nExtIter
!----------------------------------------------------------------------!
!---------------------- Solve momentum equation -----------------------!
!----------------------------------------------------------------------!
  do MIter = 1,nMIter
   
   !--- Set old variables, compute gradient ---!
   call preprocessing
      
   if (wrt_data .eq. 1) then
   write(16,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   write(17,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   write(18,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   write(20,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   endif
   !--- Compute spatial discretization ---!
   
   !--- Convective terms (1st order upwind/Central) ---!
   call convective_residual
       
   !--- Viscous terms ---!
    
   call viscous_residual
   
   !--- Source term (pressure only) ---!
   call pressure_residual
   
   !--- Apply BC ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=1,Nx
   iPoint = i + (j-1)*Nx
   if (x(i,j).ge.0.d0) then
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
   
   !--- Right outlet (i=Nx) ---!
   !Zero gradient at the wall. Only a convective flux is added. 
   !The visc contribution is zero because of zero gradient.
   i=Nx
   do j=1,Ny
   iPoint = i + (j-1)*Nx
    !--- Velocity gradient is zero ---!
    !--- Have to add convective flux as is ---!
    ! U_old(1:2,iPoint) = 0.0
     
     !--- Update residual ---!
     R(1,iPoint) = R(1,iPoint) + max(rho*U_old(1,iPoint)*dy,0.d0)*U_old(1,iPoint) ! max(rho*U_e*dy,0.d0)*U_old(1,iPoint)
     R(2,iPoint) = R(2,iPoint) + max(rho*U_old(1,iPoint)*dy,0.d0)*U_old(2,iPoint) ! max(rho*U_e*dy,0.d0)*U_old(2,iPoint)
     
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 1.0*U(1,iPoint)
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.5*U(2,iPoint)
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 0.5*U(1,iPoint)
   enddo
   
   
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=1,Nx
   iPoint = i + (j-1)*Nx
   if (x(i,j).ge.0.d0) then
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

!--- Compute residual mass flux from vel solution ---!
call compute_massflux


!----------------------------------------------------------------------!
!---------------- Solve Pressure correction equation now --------------!
!----------------------------------------------------------------------!


  

   Sol(3,:) = 0.d0
   P_Correc = 0.0
do PIter = 1,nPIter   

   R(3,:) = 0.d0  
   
   !--- Assemble pressure equation ---!
   do i=2,Nx-1
    do j=2,Ny-1
      iPoint = i + (j-1)*Nx
      Jac(iPoint,:) = 0.0
      !--- Residual of the form ---!
      !--- R = -rho*\underbar{D}_f\nabla p'*area ---!
      !--- where \underbar{D}_f is the inverse of momentum coefficient of the two points calculated at the face ---!
      !East
      jPoint = i+1 + (j-1)*Nx
      R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(dy/dx)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dy/dx)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !West
      jPoint = i-1 + (j-1)*Nx
      R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(dy/dx)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dy/dx)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !North
      R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(dx/dy)
      jPoint = i + (j+1-1)*Nx
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dx/dy)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !South
      R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(dx/dy)
      jPoint = i + (j-1-1)*Nx
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dx/dy)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
    enddo
   enddo
   
   !--- Boundary elements ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
      !East
      jPoint = i+1 + (j-1)*Nx
      R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dy/dx)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dy/dx)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !West
      jPoint = i-1 + (j-1)*Nx
      R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dy/dx)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dy/dx)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !North
      jPoint = i + (j+1-1)*Nx
      R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(dx/dy)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dx/dy)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !South
      !jPoint = i + (j-1-1)*Nx
      !Jac(iPoint,jPoint) = 
      
      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) +  
   enddo
   
   
   !--- Left inlet (i=1) ---!
   i=1
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
      !East
      jPoint = i+1 + (j-1)*Nx
      R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(dy/dx)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dy/dx)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !West
      jPoint = i+1 + (j-1)*Nx
      R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(dy/dx)
      Jac(iPoint,jPoint) = D(1,iPoint)*(dy/dx)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !North
      jPoint = i + (j+1-1)*Nx
      R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dx/dy)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dx/dy)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !South
      jPoint = i + (j-1-1)*Nx
      R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dx/dy)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dx/dy)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   enddo
   
   
   !--- Right outlet (i=Nx) ---!
   i=Nx
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
    Jac(iPoint,iPoint) = 1.0
    R(3,iPoint) = 0.d0
      !East
      !jPoint = i+1 + (j-1)*Nx
      !Jac(iPoint,jPoint) = 0.d0

      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) - Jac(iPoint,jPoint)
      !West
      !jPoint = i-1 + (j-1)*Nx
      !R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))
      !Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(Area/dx)
      
      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !North
      !jPoint = i + (j+1-1)*Nx
      !R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))
      !Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*Area/dx)
      
      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !South
      !jPoint = i + (j-1-1)*Nx
      !R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))
      !Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*Area/dx)
      
      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   enddo
   
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
      !East
      jPoint = i+1 + (j-1)*Nx
      R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dy/dx)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dy/dx)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !West
      jPoint = i-1 + (j-1)*Nx
      R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dy/dx)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dy/dx)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
      !North
      !jPoint = i + (j+1-1)*Nx
      !Jac(iPoint,jPoint) = 0.d0
      
      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) - Jac(iPoint,jPoint)
      !South
      jPoint = i + (j-1-1)*Nx
      R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(dx/dy)
      Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(dx/dy)
      
      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   enddo
   
   i=1
   j=1
   iPoint = i + (j-1)*Nx
   !East
   jPoint = i+1 + (j-1)*Nx
   R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dy/dx)
   Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dy/dx)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)   
   !North 
   jPoint = i + (j+1-1)*Nx
   R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dx/dy)
   Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dx/dy)
   
   !Inlet face(West)
   jPoint = i+1 + (j-1)*Nx !(only to compute derivative)
   R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dx/dy)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)  

   i=1
   j=Ny
   iPoint = i + (j-1)*Nx
   !East
   jPoint = i+1 + (j-1)*Nx
   R(3,iPoint) = R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dy/dx)
   Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dy/dx)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   
   !South 
   jPoint = i + (j-1-1)*Nx
   R(3,iPoint) = R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dx/dy)
   Jac(iPoint,jPoint) = 0.5*(D(1,iPoint)+D(1,jPoint))*(0.5*dx/dy)
   
   !Inlet face(West)
   jPoint = i+1 + (j-1)*Nx !(only to compute derivative)
   R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*(P_Correc(iPoint)-P_Correc(jPoint))*(0.5*dy/dx)

      Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   
   i=Nx
   j=Ny
   iPoint = i + (j-1)*Nx
   !West
   jPoint = i-1 + (j-1)*Nx
   R(3,iPoint) = 0.0!R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))
   Jac(iPoint,jPoint) = 0.0!0.5*(D(1,iPoint)+D(1,jPoint))*(Area/dx)

      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   !South 
   jPoint = i + (j-1-1)*Nx
   R(3,iPoint) = 0.0!R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))
   Jac(iPoint,jPoint) = 0.0!0.5*(D(1,iPoint)+D(1,jPoint))*(Area/dx)

      Jac(iPoint,iPoint) = 1.0!Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   
    
   i=Nx
   j=1
   iPoint = i + (j-1)*Nx
   !West
   jPoint = i-1 + (j-1)*Nx
   R(3,iPoint) = 0.0!R(3,iPoint) + 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))
   Jac(iPoint,jPoint) = 0.0!0.5*(D(1,iPoint)+D(1,jPoint))*(Area/dx)

      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   !North 
   jPoint = i + (j+1-1)*Nx
   R(3,iPoint) = 0.0!R(3,iPoint) - 0.5*(D(1,iPoint)+D(1,jPoint))*(P_Correc(iPoint)-P_Correc(jPoint))
   Jac(iPoint,jPoint) = 0.0!0.5*(D(1,iPoint)+D(1,jPoint))*(Area/dx)
   
      Jac(iPoint,iPoint) = 1.0!Jac(iPoint,iPoint) + Jac(iPoint,jPoint)
   !Jac(iPoint,:) = 0.0
   !R(3,iPoint) = 0.0
   !Jac(iPoint,iPoint) = 1.0
   !Mass(iPoint) = 0.0
   
   if (implicit_time) then
   
   mass_l2 = 0.d0
   Res_l2 = 0.0
   do iPoint=1,nPoint
     R(3,iPoint) =  R(3,iPoint) + Mass(iPoint)
     R(3,iPoint) = -R(3,iPoint)
     
     Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Vol/dt_p
          
     Res_l2 = Res_l2 + R(3,iPoint)**2.0
   enddo
   close(19)

   Res_l2 = sqrt(Res_l2/nPoint)
   if ((modulo(PIter,p_screen2)).eq.0) print*,'Res(p): ',Res_l2,PIter,ExtIter
   
   
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
        P_Correc(iPoint) = P_Correc(iPoint) - R(3,iPoint)*dt_p/Vol
        Res_l2 = Res_l2 + R(3,iPoint)**2.0
     enddo
     
     Res_l2 = sqrt(Res_l2/nPoint)
   if ((modulo(PIter,p_screen2)).eq.0) print*,'Res(p): ',Res_l2,PIter,ExtIter
   endif 
 enddo   !PIter
 
 
 !--- Correct pressure and velocities ---!
   do i=2,Nx-1
    do j=2,Ny-1
      iPoint = i + (j-1)*Nx
      P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
      !East
      jPoint = i+1 + (j-1)*Nx
      F_e(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy
      !West
      jPoint = i-1 + (j-1)*Nx
      F_w(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy
      !North
      jPoint = i + (j+1-1)*Nx
      F_n(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx
      !South
      jPoint = i + (j-1-1)*Nx
      F_s(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx
      
      U(1,iPoint) = U(1,iPoint) - (F_e(1) - F_w(1))/&
                               Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
      U(2,iPoint) = U(2,iPoint) - (F_n(1) - F_s(1))/&
                               Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)
    enddo
   enddo
   
   !--- Apply BC ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=1,Nx
    iPoint = i + (j-1)*Nx
    P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
    !P(i,j) = P(i,j+1)
    !print*,i,j,i,j+1,P(i,j+1)
   enddo
   
   !--- Left inlet (i=1) ---!
   i=1
   do j=2,Ny
    iPoint = i + (j-1)*Nx
    P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
    !P(i,j) = P(i+1,j)
    !print*,i,j,i+1,j,P(i+1,j)
   enddo
   !--- Move down perhaps ---!
   !--- Right outlet (i=Nx) ---!
   i=Nx
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    !P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
    P(i,j) = P_outlet
    
    !Should add vel here
    !East
    jPoint = i+1 + (j-1)*Nx
    F_e(1) = 0.0!0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy
    !West
    jPoint = i-1 + (j-1)*Nx
    F_w(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy
    !North
    jPoint = i + (j+1-1)*Nx
    F_n(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx/2.0
    !South
    jPoint = i + (j-1-1)*Nx
    F_s(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx/2.0
      
    U(1,iPoint) = U(1,iPoint) - (F_e(1) - F_w(1))/&
                               Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
    U(2,iPoint) = U(2,iPoint) - (F_n(1) - F_s(1))/&
                               Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)
    !print*,i,j,i-1,j,P(i-1,j)
   enddo
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
    !P(i,j) = P(i,j-1)
    !print*,i,j,i,j-1,P(i,j-1)
   enddo
   !print*,'----Next iteration----'
   !--- Convergence monitoring ---!
   
   !--- Output solution ---!
   if (modulo(ExtIter,1) .eq. 0) then
   
   open(unit=14,file='../out/Centerline_channel.txt',status='unknown')
     i = (Nx+1)/2
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(14,*) y(i,j),U(1,iPoint),U(2,iPoint),P(i,j),iPoint,x(i,j)
     enddo
    close(14)
    open(unit=24,file='../out/Start_channel.txt',status='unknown')
     i = k+1
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(24,*) y(i,j),U(1,iPoint),U(2,iPoint),P(i,j),iPoint,x(i,j)
     enddo
    close(24)
    open(unit=34,file='../out/Outlet_channel.txt',status='unknown')
     i = Nx
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(34,*) y(i,j),U(1,iPoint),U(2,iPoint),P(i,j),iPoint,x(i,j)
     enddo
    close(34)
    open(unit=34,file='../out/Interior_channel.txt',status='unknown')
     i = Nx-5
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(34,*) y(i,j),U(1,iPoint),U(2,iPoint),P(i,j),iPoint,x(i,j)
     enddo
    close(34)
   endif
enddo !ExtIter


open(unit=13,file='../out/Solution_channel.txt',status='unknown')
do j = 1,Ny
 do i=1,Nx
  iPoint = i+(j-1)*Nx
  write(13,*) x(i,j),y(i,j),P(i,j),U(1,iPoint),U(2,iPoint)
 enddo
enddo

close(11)
close(12)
if (wrt_data .eq. 1) then
close(16)
close(17)
close(18)
endif
end program iNS







