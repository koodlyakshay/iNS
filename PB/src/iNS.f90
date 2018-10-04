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
  open(unit=16,file='out/Convbudgetspb.txt',status='unknown')
  open(unit=17,file='out/Viscbudgetspb.txt',status='unknown')
  open(unit=18,file='out/presbudgetspb.txt',status='unknown')
endif

!--- Begin Solver ---!

do ExtIter = 1,nExtIter
!----------------------------------------------------------------------!
!---------------------- Solve momentum equation -----------------------!
!----------------------------------------------------------------------!
  do MIter = 1,nMIter
   R = 0.0
   Tot_Jac = 0.0
   
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
   
   
   if (wrt_data .eq. 1) then
   write(16,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   write(17,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   write(18,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   write(20,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   endif
   !--- Compute spatial discretization ---!
   
   do i=1,Nx
    do j=1,Ny
    !--- Node definition ---!
    iPoint = i + (j-1)*Nx
    
    !--- Convective terms (1st order upwind/Central) ---!    
    F_e = 0.d0
    F_w = 0.d0
    F_n = 0.d0
    F_s = 0.d0
    if (i.ne.Nx) then
       !--- Compute flux and update Jacobians ---!
       !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       V_e = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       if (upwind) then
       F_e(1) = max(rho*U_e*dy,0.d0)*U_old(1,iPoint) + min(rho*U_e*dy,0.d0)*U_old(1,jPoint) ! (rho*U)*U|_e
       F_e(2) = max(rho*U_e*dy,0.d0)*U_old(2,iPoint) + min(rho*U_e*dy,0.d0)*U_old(2,jPoint) ! (rho*V)*U|_e
       else
       F_e(1) = rho*U_e*U_e*dy          ! (rho*U)*U|_e
       F_e(2) = rho*V_e*U_e*dy          ! (rho*V)*U|_e
       endif
       lambda_i = abs(2.d0*U_old(1,iPoint)*dy)
       lambda_j = abs(2.d0*U_old(1,jPoint)*dy)
       
       lambda_mean = 0.5*(lambda_i + lambda_j)
       if (lambda_mean .eq. 0.d0) then
        lambda_i = abs(U_old(1,iPoint)*dy) + sqrt(U_old(1,iPoint)*U_old(1,iPoint)*dy*dy + (artvisc/rho)*dy*dy) 
        lambda_j = abs(U_old(1,jPoint)*dy) + sqrt(U_old(1,jPoint)*U_old(1,jPoint)*dy*dy + (artvisc/rho)*dy*dy) 
        lambda_mean = 0.5*(lambda_i + lambda_j)
       endif
       Phi_i = (lambda_i/(4.d0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.d0*lambda_mean))**Param_p
       if ((Phi_i + Phi_j).eq.0.d0) then
           SF = 0.d0
       else 
           SF = 4.d0*Phi_i*Phi_j/(Phi_i + Phi_j)
       endif
       sc0 = 3.d0*(4.d0+4.d0)/(4.d0*4.d0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.d0/3.d0
       F_e(1) = F_e(1) + E_0*(U_old(1,iPoint) - U_old(1,jPoint))*SF*lambda_mean
       F_e(2) = F_e(2) + E_0*(U_old(2,iPoint) - U_old(2,jPoint))*SF*lambda_mean
       
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) + 1.0*U_e*dy - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) + 0.5*V_e*dy
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) + 0.5*U_e*dy - E_0*SF*lambda_mean
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 1.0*U_e*dy + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.5*V_e*dy
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 0.5*U_e*dy + E_0*SF*lambda_mean
    endif  
    
    if (i.ne.1) then
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       V_w = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       if (upwind) then
       F_w(1) = -max(rho*U_w*dy,0.d0)*U_old(1,jPoint) - min(rho*U_w*dy,0.d0)*U_old(1,iPoint) ! (rho*U)*U|_w
       F_w(2) = -max(rho*U_w*dy,0.d0)*U_old(2,jPoint) - min(rho*U_w*dy,0.d0)*U_old(2,iPoint) ! (rho*V)*U|_w
       else
       F_w(1) = -rho*U_w*U_w*dy          ! (rho*U)*U|_e
       F_w(2) = -rho*V_w*U_w*dy          ! (rho*V)*U|_e
       endif
       lambda_i = abs(2.d0*U_old(1,iPoint)*dy)
       lambda_j = abs(2.d0*U_old(1,jPoint)*dy)
       
       lambda_mean = 0.5*(lambda_i + lambda_j)
       if (lambda_mean .eq. 0.d0) then
        lambda_i = abs(U_old(1,iPoint)*dy) + sqrt(U_old(1,iPoint)*U_old(1,iPoint)*dy*dy + (artvisc/rho)*dy*dy) 
        lambda_j = abs(U_old(1,jPoint)*dy) + sqrt(U_old(1,jPoint)*U_old(1,jPoint)*dy*dy + (artvisc/rho)*dy*dy) 
        lambda_mean = 0.5*(lambda_i + lambda_j)
       endif
       Phi_i = (lambda_i/(4.d0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.d0*lambda_mean))**Param_p
       if ((Phi_i + Phi_j).eq.0.d0) then 
           SF = 0.d0
       else 
            SF = 4.d0*Phi_i*Phi_j/(Phi_i + Phi_j)
       endif
       sc0 = 3.d0*(4.d0+4.d0)/(4.d0*4.d0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.d0/3.d0
       F_w(1) = F_w(1) + E_0*(U_old(1,iPoint) - U_old(1,jPoint))*SF*lambda_mean
       F_w(2) = F_w(2) + E_0*(U_old(2,iPoint) - U_old(2,jPoint))*SF*lambda_mean
       
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 1.0*U_w*dy - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.5*V_w*dy
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - 0.5*U_w*dy - E_0*SF*lambda_mean
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) - 1.0*U_w*dy + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) - 0.5*V_w*dy
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) - 0.5*U_w*dy + E_0*SF*lambda_mean
    endif  
    
    if (j.ne.Ny) then
       !North 
       jPoint = i + (j+1-1)*Nx
       U_n = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       V_n = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       if (upwind) then
       F_n(1) = max(rho*V_n*dx,0.d0)*U_old(1,iPoint) + min(rho*V_n*dx,0.d0)*U_old(1,jPoint) ! (rho*U)*V|_n
       F_n(2) = max(rho*V_n*dx,0.d0)*U_old(2,iPoint) + min(rho*V_n*dx,0.d0)*U_old(2,jPoint) ! (rho*V)*V|_n
       else
       F_n(1) = rho*V_n*U_n*dx          ! (rho*U)*U|_e
       F_n(2) = rho*V_n*V_n*dx          ! (rho*V)*U|_e
       endif
       lambda_i = abs(2.d0*U_old(2,iPoint)*dx)
       lambda_j = abs(2.d0*U_old(2,jPoint)*dx)
       
       lambda_mean = 0.5*(lambda_i + lambda_j)
       if (lambda_mean .eq. 0.d0) then
        lambda_i = abs(U_old(2,iPoint)*dx) + sqrt(U_old(2,iPoint)*U_old(2,iPoint)*dx*dx + (artvisc/rho)*dx*dx) 
        lambda_j = abs(U_old(2,jPoint)*dx) + sqrt(U_old(2,jPoint)*U_old(2,jPoint)*dx*dx + (artvisc/rho)*dx*dx)
        lambda_mean = 0.5*(lambda_i + lambda_j)
       endif
       Phi_i = (lambda_i/(4.d0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.d0*lambda_mean))**Param_p
       if ((Phi_i + Phi_j).eq.0.d0) then 
           SF = 0.d0
       else 
            SF = 4.d0*Phi_i*Phi_j/(Phi_i + Phi_j)
       endif
       sc0 = 3.d0*(4.d0+4.d0)/(4.d0*4.d0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.d0/3.d0
       F_n(1) = F_n(1) + E_0*(U_old(1,iPoint) - U_old(1,jPoint))*SF*lambda_mean
       F_n(2) = F_n(2) + E_0*(U_old(2,iPoint) - U_old(2,jPoint))*SF*lambda_mean
       
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) + 0.5*V_n*dx - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) + 0.5*U_n*dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) + 1.0*V_n*dx - E_0*SF*lambda_mean
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.5*V_n*dx + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.5*U_n*dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 1.0*V_n*dx + E_0*SF*lambda_mean
    endif
    
    if (j.ne.1) then   
       !South
       jPoint = i + (j-1-1)*Nx
       U_s = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       V_s = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       if (upwind) then
       F_s(1) = -max(rho*V_s*dx,0.d0)*U_old(1,jPoint) - min(rho*V_s*dx,0.d0)*U_old(1,iPoint) ! (rho*U)*V|_s
       F_s(2) = -max(rho*V_s*dx,0.d0)*U_old(2,jPoint) - min(rho*V_s*dx,0.d0)*U_old(2,iPoint) ! (rho*V)*V|_s
       else
       F_s(1) = -rho*V_s*U_s*dx          ! (rho*U)*U|_e
       F_s(2) = -rho*V_s*V_s*dx          ! (rho*V)*U|_e
       endif
       lambda_i = abs(2.d0*U_old(2,iPoint)*dx)
       lambda_j = abs(2.d0*U_old(2,jPoint)*dx)
       
       lambda_mean = 0.5*(lambda_i + lambda_j)
       if (lambda_mean .eq. 0.d0) then
        lambda_i = abs(U_old(2,iPoint)*dx) + sqrt(U_old(2,iPoint)*U_old(2,iPoint)*dx*dx + (artvisc/rho)*dx*dx) 
        lambda_j = abs(U_old(2,jPoint)*dx) + sqrt(U_old(2,jPoint)*U_old(2,jPoint)*dx*dx + (artvisc/rho)*dx*dx)
        lambda_mean = 0.5*(lambda_i + lambda_j)
       endif
       Phi_i = (lambda_i/(4.d0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.d0*lambda_mean))**Param_p
       if ((Phi_i + Phi_j).eq.0.d0) then 
           SF = 0.d0
       else 
            SF = 4.d0*Phi_i*Phi_j/(Phi_i + Phi_j)
       endif
       sc0 = 3.d0*(4.d0+4.d0)/(4.d0*4.d0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.d0/3.d0
       F_s(1) = F_s(1) + E_0*(U_old(1,iPoint) - U_old(1,jPoint))*SF*lambda_mean
       F_s(2) = F_s(2) + E_0*(U_old(2,iPoint) - U_old(2,jPoint))*SF*lambda_mean
       
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.5*V_s*dx - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.5*U_s*dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - 1.0*V_s*dx - E_0*SF*lambda_mean
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) - 0.5*V_s*dx + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) - 0.5*U_s*dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) - 1.0*V_s*dx + E_0*SF*lambda_mean
    endif    
       !--- Update residuals ---!
       
       Fc(1,i,j) = F_e(1) + F_w(1) + F_n(1) + F_s(1) ! (rho*U)*U|_e - (rho*U)*U|_w + (rho*U)*V|_n - (rho*U)*V|_s
       Fc(2,i,j) = F_e(2) + F_w(2) + F_n(2) + F_s(2) ! (rho*V)*U|_e - (rho*V)*U|_w + (rho*V)*V|_n - (rho*V)*V|_s
       
       if (wrt_data .eq. 1) then
       write(16,*)F_e(1:2),'iPoint',iPoint,'jPoint',i+1 + (j-1)*Nx,U_e,V_e
       write(16,*)F_w(1:2),'iPoint',iPoint,'jPoint',i-1 + (j-1)*Nx,U_w,V_w
       write(16,*)F_n(1:2),'iPoint',iPoint,'jPoint',i + (j+1-1)*Nx,U_n,V_n
       write(16,*)F_s(1:2),'iPoint',iPoint,'jPoint',i + (j-1-1)*Nx,U_s,V_s
       write(16,*)
       endif
       
    !--- Viscous terms ---!
    
    F_e = 0.d0
    F_w = 0.d0
    F_n = 0.d0
    F_s = 0.d0
    
    if (i.ne.Nx) then
       !East
       Gr_e = 2.0*0.5*( GradU(1,1,i+1,j) + GradU(1,1,i,j))
       F_e(1) = mu*Gr_e*dy ! du/dx|_e
       Gr_e = 1.0*0.5*( GradU(2,1,i+1,j) + GradU(2,1,i,j)) + 1.0*0.5*( GradU(1,2,i+1,j) + GradU(1,2,i,j))
       F_e(2) = mu*Gr_e*dy ! dv/dx|_e
    endif
    if (i.ne.1) then
       !West
       Gr_w = 2.0*0.5*( GradU(1,1,i-1,j) + GradU(1,1,i,j)) 
       F_w(1) = -mu*Gr_w*dy ! du/dx|_w
       Gr_w = 1.0*0.5*( GradU(2,1,i-1,j) + GradU(2,1,i,j)) + 1.0*0.5*( GradU(1,2,i-1,j) + GradU(1,2,i,j)) 
       F_w(2) = -mu*Gr_w*dy ! dv/dx|_w
    endif
    if (j.ne.Ny) then
       !North
       Gr_n = 1.0*0.5*( GradU(1,2,i,j+1) + GradU(1,2,i,j)) + 1.0*0.5*( GradU(2,1,i,j+1) + GradU(2,1,i,j)) 
       F_n(1) = mu*Gr_n*dx ! du/dy|_n
       Gr_n = 2.0*0.5*( GradU(2,2,i,j+1) + GradU(2,2,i,j)) 
       F_n(2) = mu*Gr_n*dx ! dv/dy|_n
    endif
    if (j.ne.1) then
       !South
       Gr_s = 1.0*0.5*( GradU(1,2,i,j-1) + GradU(1,2,i,j)) + 1.0*0.5*( GradU(2,1,i,j-1) + GradU(2,1,i,j))
       F_s(1) = -mu*Gr_s*dx ! du/dy|_s
       Gr_s = 2.0*0.5*( GradU(2,2,i,j-1) + GradU(2,2,i,j)) 
       F_s(2) = -mu*Gr_s*dx ! dv/dy|_s
    endif       
       Fv(1,i,j) = F_e(1) + F_w(1) + F_n(1) + F_s(1) ! du/dx|_e - du/dx|_w + du/dy|_n - du/dy|_s
       Fv(2,i,j) = F_e(2) + F_w(2) + F_n(2) + F_s(2) ! dv/dx|_e - dv/dx|_w + dv/dy|_n - dv/dy|_s
       
       if(wrt_data .eq. 1) then
       write(17,*)F_e(1:2),'iPoint',iPoint,'jPoint',i+1 + (j-1)*Nx,Gr_e
       write(17,*)F_w(1:2),'iPoint',iPoint,'jPoint',i-1 + (j-1)*Nx,Gr_w
       write(17,*)F_n(1:2),'iPoint',iPoint,'jPoint',i + (j+1-1)*Nx,Gr_n
       write(17,*)F_s(1:2),'iPoint',iPoint,'jPoint',i + (j-1-1)*Nx,Gr_s
       write(17,*)
       endif
       
       !--- Update Jacobians ---!
    if (i.ne.Nx) then
       !East
       jPoint = i+1 + (j-1)*Nx
       
       !print*,iPoint,jPoint,F_e,1.0*mu*area/dx,'ve',x(i,j),y(i,j) 
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - mu*area/dx
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - 0.5*mu*area/dx
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + mu*area/dx
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 0.5*mu*area/dx
    endif
    if (i.ne.1) then
       !West
       jPoint = i-1 + (j-1)*Nx
       
       !print*,iPoint,jPoint,F_w,1.0*mu*area/dx,'vw',x(i,j),y(i,j)
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - mu*area/dx
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - 0.5*mu*area/dx
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + mu*area/dx
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 0.5*mu*area/dx
    endif   
    if (j.ne.Ny) then
       !North
       jPoint = i + (j+1-1)*Nx
       
       !print*,iPoint,jPoint,F_n,1.0*mu*area/dx,'vn',x(i,j),y(i,j)
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - mu*area/dx
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + mu*area/dx
    endif   
    if (j.ne.1) then
       !South
       jPoint = i + (j-1-1)*Nx
       
       !print*,iPoint,jPoint,F_s,1.0*mu*area/dx,'vs',x(i,j),y(i,j)
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - mu*area/dx
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.5*mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + mu*area/dx
    endif
       !--- Update residual ---!
       
       R(1:2,iPoint) = Fc(1:2,i,j) !- Fv(1:2,i,j)
       
       !--- Source term (pressure only) ---!
    if (i.ne.Nx) then
       !East
       F_e(1) = 0.5*(P(i+1,j) + P(i,j))*dy ! p*Area|_e
       F_e(2) = 0.0
       !print*,iPoint,(i+1 + (j-1)*Nx),F_e,'p'
    endif
    if (i.ne.1) then
       !West
       F_w(1) = -0.5*(P(i,j) + P(i-1,j))*dy ! p*Area|_w
       F_s(2) = 0.0
       !print*,iPoint,(i-1 + (j-1)*Nx),F_w,'p'
    endif
    if (j.ne.Ny) then
       !North
       F_n(1) = 0.0
       F_n(2) = 0.5*(P(i,j+1) + P(i,j))*dx ! p*Area|_n
       !print*,iPoint,(i + (j+1-1)*Nx),F_n,'p'
    endif
    if (j.ne.1) then
       !South
       F_s(1) = 0.0
       F_s(2) = -0.5*(P(i,j) + P(i,j-1))*dx ! p*Area|_s
       !print*,iPoint,(i + (j-1-1)*Nx),F_s,'p'
    endif
       R(1,iPoint) = R(1,iPoint) + GradU(3,1,i,j)*Vol!F_e(1) + F_w(1) + F_n(1) + F_s(1) ! p*Area|_e - p*Area|_w
       R(2,iPoint) = R(2,iPoint) + GradU(3,2,i,j)*Vol!F_e(2) + F_w(2) + F_n(2) + F_s(2) ! p*Area|_n - p*Area|_s
       
       if (wrt_data .eq. 1) then
       write(18,*)F_e(1:2),'iPoint',iPoint,'jPoint',i+1 + (j-1)*Nx
       write(18,*)F_w(1:2),'iPoint',iPoint,'jPoint',i-1 + (j-1)*Nx
       write(18,*)F_n(1:2),'iPoint',iPoint,'jPoint',i + (j+1-1)*Nx
       write(18,*)F_s(1:2),'iPoint',iPoint,'jPoint',i + (j-1-1)*Nx
       write(18,*)
       endif
       
       if(wrt_data .eq. 1) then
       write(20,*)iPoint
       write(20,*)R(1,iPoint)
       write(20,*)R(2,iPoint)
       write(20,*)
       endif
 
    enddo
   enddo
   
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
     R(1,iPoint) = R(1,iPoint) + rho*U_old(1,iPoint)*dy*U_old(1,iPoint) ! (rho*U)*U|_e
     R(2,iPoint) = R(2,iPoint) + rho*U_old(1,iPoint)*dy*U_old(2,iPoint) ! (rho*V)*U|_e
     
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

!----------------------------------------------------------------------!
!---------------- Solve Pressure correction equation now --------------!
!----------------------------------------------------------------------!


!--- Compute residual mass flux from vel solution ---!
Mass = 0.d0
P_Correc = 0.d0
   do i=2,Nx-1
    do j=2,Ny-1
    !--- Point definition ---!
    iPoint = i + (j-1)*Nx
    
    !--- Mass flux ---!
    !--- Rhie-Chow interpolation ---!
       !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy + 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i-1,j) - P(i,j))/dx
       F_w(1) = -rho*U_w*dy - 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx + 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j-1) - P(i,j))/dy
       F_s(1) = -rho*V_s*dx - 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx
              
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
       !print*,iPoint,F_e,F_w,F_n,F_s,Mass(iPoint)
    enddo
   enddo
   

   !--- Boundary elements ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=2,Nx-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*Vol/(alfa*dt_m)
       !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy/2.d0 + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.d0
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i-1,j) - P(i,j))/dx
       F_w(1) = -rho*U_w*dy/2.d0 - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.d0
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx
       !South
       F_s(1) = 0.d0
               
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
   enddo
   
   !--- Left inlet (i=1) ---!
   i=1
   do j=2,Ny-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*Vol/(alfa*dt_m)
      !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       !West
       F_w(1) = -rho*U(1,iPoint)*dy
       ! The mass flux is known at the inlet.
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx/2.d0 + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx/2.0
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j-1) - P(i,j))/dy
       F_s(1) = -rho*V_s*dx/2.d0 - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx/2.0
              
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
   enddo
   
   !--- Right outlet (i=Nx) ---!
   ! pressure is specified at the outlet. 
   ! No need to solve the equations.
   i=Nx
   do j=2,Ny-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*Vol/(alfa*dt_m)
      !East 
       F_e(1) = 0.d0
       
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i-1,j) - P(i,j))/dx
       F_w(1) = -rho*U_w*dy - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx/2.d0 + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx/2.0
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j-1) - P(i,j))/dy
       F_s(1) = -rho*V_s*dx/2.d0 - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx/2.0
              
     Mass(iPoint) = 0.d0!F_e(1) + F_w(1) + F_n(1) + F_s(1)
   enddo
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*Vol/(alfa*dt_m)
      !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy/2.d0 + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
       
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i-1,j) - P(i,j))/dx
       F_w(1) = -rho*U_w*dy/2.d0 - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
       
       !North 
       F_n(1) = 0.d0
       
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j-1) - P(i,j))/dy
       F_s(1) = -rho*V_s*dx - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx
              
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
   enddo
   
   !--- Corner elements ---!
   
   i=1
   j=1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.25*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.25*Vol/(alfa*dt_m)
   !East
   jPoint = i+1 + (j-1)*Nx
   U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
   GradPav = (P(i+1,j) - P(i,j))/dx
   F_e(1) = rho*U_e*dy/2.d0 + 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !North 
   jPoint = i + (j+1-1)*Nx
   V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
   GradPav = (P(i,j+1) - P(i,j))/dy
   F_n(1) = rho*V_n*dx/2.d0 + 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dx/2.0
    
    Mass(iPoint) = F_e(1) + F_n(1) - rho*U(1,iPoint)*dy/2 !Contribution from inlet
    
   i=1
   j=Ny
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.25*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.25*Vol/(alfa*dt_m)
   !East
   jPoint = i+1 + (j-1)*Nx
   U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
   GradPav = (P(i+1,j) - P(i,j))/dx
   F_e(1) = rho*U_e*dy/2.d0 + 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !South 
   jPoint = i + (j-1-1)*Nx
   V_s = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j-1) + GradU(3,2,i,j))
   GradPav = (P(i,j-1) - P(i,j))/dy
   F_s(1) = -rho*V_s*dx/2.d0 - 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dx/2.0
    
    Mass(iPoint) = F_e(1) + F_s(1) - rho*U(1,iPoint)*dy/2 !Contribution from inlet
    
   i=Nx
   j=Ny
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.25*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.25*Vol/(alfa*dt_m)
   !West
   jPoint = i-1 + (j-1)*Nx
   U_w = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i-1,j) + GradU(3,1,i,j))
   GradPav = (P(i-1,j) - P(i,j))/dx
   F_w(1) = -rho*U_w*dy/2.d0 - 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !South 
   jPoint = i + (j-1-1)*Nx
   V_s = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j-1) + GradU(3,2,i,j))
   GradPav = (P(i,j-1) - P(i,j))/dy
   F_s(1) = -rho*V_s*dx/2.d0 - 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dx/2.0
    
   Mass(iPoint) = 0.d0!F_w(1) + F_s(1)
   !The outlet element has specified pressure, no need to solve equation.
    
    
   i=Nx
   j=1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.25*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.25*Vol/(alfa*dt_m)
   !West
   jPoint = i-1 + (j-1)*Nx
   U_w = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i-1,j) + GradU(3,1,i,j))
   GradPav = (P(i-1,j) - P(i,j))/dx
   F_w(1) = -rho*U_w*dy/2.d0 - 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !North 
   jPoint = i + (j+1-1)*Nx
   V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
   GradPav = (P(i,j+1) - P(i,j))/dy
   F_n(1) = rho*V_n*dx/2.d0 + 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dx/2.0
    
   Mass(iPoint) = 0.d0!F_w(1) + F_n(1)
   !The outlet element has specified pressure, no need to solve equation.
   
   mass_l2 = 0.d0

   do iPoint=1,nPoint
      mass_l2 = mass_l2 + Mass(iPoint)**2.0
   enddo
   mass_l2 = sqrt(mass_l2/nPoint)
   print*,'Mass: ',mass_l2,ExtIter

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
   
   open(unit=14,file='out/Centerline_channel.txt',status='unknown')
     i = (Nx+1)/2
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(14,*) y(i,j),U(1,iPoint),U(2,iPoint),P(i,j),iPoint
     enddo
    close(14)
    open(unit=24,file='out/Start_channel.txt',status='unknown')
     i = k+1
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(24,*) y(i,j),U(1,iPoint),U(2,iPoint),P(i,j),iPoint
     enddo
    close(24)
    open(unit=34,file='out/Outlet_channel.txt',status='unknown')
     i = Nx
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(34,*) y(i,j),U(1,iPoint),U(2,iPoint),P(i,j),iPoint
     enddo
    close(34)
    open(unit=34,file='out/Interior_channel.txt',status='unknown')
     i = Nx-4
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(34,*) y(i,j),U(1,iPoint),U(2,iPoint),P(i,j),iPoint
     enddo
    close(34)
   endif
enddo !ExtIter


open(unit=13,file='out/Solution.txt',status='unknown')
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







