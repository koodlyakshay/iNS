subroutine pressure_correc_eqn(npiter, extiter, implicit_time)

use global_vars
use output_vars
use flow_parmaters


implicit none

integer           :: iPoint, nPoint, PIter
integer           :: ExtIter, npiter
integer           :: convergence, liniter
real              :: res_l2
logical           :: implicit_time

!---------------- Solve Pressure correction equation now --------------!
!----------------------------------------------------------------------!
  nPoint = Nx*Ny
  Sol(3,:) = 0.0
  P_Correc = 0.0
  GradPc = 0.0
  liniter = 0
  
  do PIter = 1,nPIter
  
     R(3,:) = 0.0  
     
     !--- Compute pressure correction gradient ---!
     call compute_pcgradientgg
     
     !--- Compute residual for the pressure correction equation ---!
     call compute_pcresidual
  
     call boundary_conditions
     
        
     if (implicit_time) then
     
     Res_l2 = 0.0

     do iPoint=1,nPoint
       R(3,iPoint) =  R(3,iPoint) + Mass(iPoint)
       R(3,iPoint) = -R(3,iPoint)
       
       Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + Vol(iPoint)/dt_p
            
       Res_l2 = Res_l2 + R(3,iPoint)**2.0
     enddo
  
     Res_l2 = sqrt(Res_l2/nPoint)
     if ((modulo(PIter,p_screen2)).eq.0) print*,'Res(p): ',log10(Res_l2),PIter,ExtIter
        
     !--- Solve pressure correction equation ---!
     call seidel(0,nPoint,Jac(:,:),R(3,:),1.0,Sol(3,:),Res(:),liniter,convergence)
     if ((convergence .ne. 0) .and. (convergence .ne. 4)) print*, 'Error in p',convergence,ExtIter
     
     do iPoint = 1,nPoint
       P_Correc(iPoint) = P_Correc(iPoint) + Sol(3,iPoint)
     enddo
     
     !--- Time Integration (Explicit) ---!
     else  
       Res_l2 = 0.0
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

end subroutine pressure_correc_eqn


subroutine compute_pcresidual

!--- Residual of the form ---!
!--- R = -rho*\underbar{D}_f\nabla p'*area ---!
!--- where \underbar{D}_f is the inverse of momentum coefficient of the two points calculated at the face ---!

use global_vars
implicit none

integer     :: i, j, iPoint, jPoint

!--- Assemble pressure equation ---!
 do i=1,Nx
   do j=1,Ny
     iPoint = i + (j-1)*Nx
     Jac(iPoint,:) = 0.0
      
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

end subroutine compute_pcresidual


subroutine boundary_conditions

use global_vars

implicit none

integer     :: i, j, iPoint, jPoint

!--- Boundary elements ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
      !South
      R(3,iPoint) = R(3,iPoint) + D(1,iPoint)*GradPc(2,i,j)*dx
      
      !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)

   enddo
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
     !North
     R(3,iPoint) = R(3,iPoint) - D(1,iPoint)*GradPc(2,i,j)*dx
     
     !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)

   enddo
   
   !--- Left inlet (i=1) ---!
   i=1
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    Jac(iPoint,:) = 0.0
    !West
    R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*GradPc(1,i,j)*dy
    !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)

   enddo
   
   
   !--- Right outlet (i=Nx) ---!
   i=Nx
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    !East
    R(3,iPoint) = R(3,iPoint) - (D(1,iPoint))*GradPc(1,i,j)*dy
    !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)

   enddo
   
   
   i=1
   j=1
   iPoint = i + (j-1)*Nx
 
   !Inlet face(West)
   R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*GradPc(1,i,j)*dy/2.0

   !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy) 

   
   !South
   R(3,iPoint) = R(3,iPoint) + D(1,iPoint)*GradPc(2,i,j)*dx/2.0
   !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)
   

   i=1
   j=Ny
   iPoint = i + (j-1)*Nx

   !Inlet face(West)
   R(3,iPoint) = R(3,iPoint) + (D(1,iPoint))*GradPc(1,i,j)*dy/2.0

   !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)
          
   !North
   R(3,iPoint) = R(3,iPoint) - D(1,iPoint)*GradPc(2,i,j)*dx/2.0
      
   !Jac(iPoint,iPoint) = Jac(iPoint,iPoint)+ D(2,iPoint)*(dx/dy)
      
   i=Nx
   j=Ny
   iPoint = i + (j-1)*Nx

   !East (outlet)
   R(3,iPoint) = R(3,iPoint) - (D(1,iPoint))*GradPc(1,i,j)*dy/2.0
   
   !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)
      
   !North
   R(3,iPoint) = R(3,iPoint) - D(1,iPoint)*GradPc(2,i,j)*dx/2.0
   
   !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)
      
    
   i=Nx
   j=1
   iPoint = i + (j-1)*Nx

   !East
   R(3,iPoint) = R(3,iPoint) - (D(1,iPoint))*GradPc(1,i,j)*dy/2.0
   
   !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)
      
   !South
   R(3,iPoint) = R(3,iPoint) + D(1,iPoint)*GradPc(2,i,j)*dx/2.0
         
   !Jac(iPoint,iPoint) = Jac(iPoint,iPoint) + D(2,iPoint)*(dx/dy)


end subroutine boundary_conditions
