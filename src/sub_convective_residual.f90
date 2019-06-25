!> \file sub_convective_residual.f90
!! \brief Subroutine to compute residual contribution from convection terms of the momentum equations.

subroutine convective_residual(upwind, muscl)

use global_vars

implicit none

integer           :: i,j,iPoint,jpoint
real              :: FaceFlux
real              :: lambda_i,lambda_j
real              :: lambda_mean
real              :: E_0
real              :: Phi_i
real              :: Phi_j
logical           :: upwind, muscl


do i=1,Nx
    do j=1,Ny
    !--- Node definition ---!
    iPoint = i + (j-1)*Nx
    F_e = 0.0
    F_w = 0.0
    F_n = 0.0
    F_s = 0.0
    if (i.ne.Nx) then
       !--- Compute flux and update Jacobians ---!
       !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       V_e = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       FaceFlux = rho*U_e*dy
       if ((j.eq.1).or.(j.eq.Ny)) FaceFlux = rho*U_e*dy/2.0
       if (upwind) then
         if (FaceFlux .gt. 0.0) then
           U_up = U_old(1,iPoint)
           V_up = U_old(2,iPoint)
           if (muscl) then
             U_up = U_up + GradU(1,1,i,j)*dx/2.0
             V_up = V_up + GradU(2,1,i,j)*dx/2.0
           endif
         else
           U_up = U_old(1,jPoint)
           V_up = U_old(2,jPoint)
           if (muscl) then
             U_up = U_up - GradU(1,1,i+1,j)*dx/2.0
             V_up = V_up - GradU(2,1,i+1,j)*dx/2.0
           endif
         endif
       F_e(1) = FaceFlux*U_up ! (rho*U)*U|_e
       F_e(2) = FaceFlux*V_up ! (rho*V)*U|_e
       else
       F_e(1) = rho*U_e*U_e*dy          ! (rho*U)*U|_e
       F_e(2) = rho*V_e*U_e*dy          ! (rho*V)*U|_e
       endif
       lambda_i = abs(2.0*U_old(1,iPoint)*dy)
       lambda_j = abs(2.0*U_old(1,jPoint)*dy)
       
       lambda_mean = 0.5*(lambda_i + lambda_j)
       if (lambda_mean .eq. 0.0) then
        lambda_i = abs(U_old(1,iPoint)*dy) + sqrt(U_old(1,iPoint)*U_old(1,iPoint)*dy*dy + (artvisc/rho)*dy*dy) 
        lambda_j = abs(U_old(1,jPoint)*dy) + sqrt(U_old(1,jPoint)*U_old(1,jPoint)*dy*dy + (artvisc/rho)*dy*dy) 
        lambda_mean = 0.5*(lambda_i + lambda_j)
       endif
       Phi_i = (lambda_i/(4.0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.0*lambda_mean))**Param_p
       if ((Phi_i + Phi_j).eq.0.0) then
           SF = 0.0
       else 
           SF = 4.0*Phi_i*Phi_j/(Phi_i + Phi_j)
       endif
       sc0 = 3.0*(4.0+4.0)/(4.0*4.0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.0/3.0
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
       FaceFlux = rho*U_w*dy
       if ((j.eq.1).or.(j.eq.Ny)) FaceFlux = rho*U_w*dy/2.0
       if (FaceFlux .gt. 0.0) then
          U_up = U_old(1,jPoint)
          V_up = U_old(2,jPoint)
          if (muscl) then
             U_up = U_up + GradU(1,1,i-1,j)*dx/2.0
             V_up = V_up + GradU(2,1,i-1,j)*dx/2.0
          endif
       else
          U_up = U_old(1,iPoint)
          V_up = U_old(2,iPoint)
          if (muscl) then
             U_up = U_up - GradU(1,1,i,j)*dx/2.0
             V_up = V_up - GradU(2,1,i,j)*dx/2.0
          endif
       endif
       F_w(1) = -FaceFlux*U_up ! (rho*U)*U|_w
       F_w(2) = -FaceFlux*V_up ! (rho*U)*V|_w
       else
       F_w(1) = -rho*U_w*U_w*dy          ! (rho*U)*U|_w
       F_w(2) = -rho*V_w*U_w*dy          ! (rho*U)*V|_w
       endif
       lambda_i = abs(2.0*U_old(1,iPoint)*dy)
       lambda_j = abs(2.0*U_old(1,jPoint)*dy)
       
       lambda_mean = 0.5*(lambda_i + lambda_j)
       if (lambda_mean .eq. 0.0) then
        lambda_i = abs(U_old(1,iPoint)*dy) + sqrt(U_old(1,iPoint)*U_old(1,iPoint)*dy*dy + (artvisc/rho)*dy*dy) 
        lambda_j = abs(U_old(1,jPoint)*dy) + sqrt(U_old(1,jPoint)*U_old(1,jPoint)*dy*dy + (artvisc/rho)*dy*dy) 
        lambda_mean = 0.5*(lambda_i + lambda_j)
       endif
       Phi_i = (lambda_i/(4.0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.0*lambda_mean))**Param_p
       if ((Phi_i + Phi_j).eq.0.0) then 
           SF = 0.0
       else 
            SF = 4.0*Phi_i*Phi_j/(Phi_i + Phi_j)
       endif
       sc0 = 3.0*(4.0+4.0)/(4.0*4.0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.0/3.0
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
       FaceFlux = rho*V_n*dx
       if ((i.eq.1).or.(i.eq.Nx)) FaceFlux = rho*V_n*dx/2.0
       if (FaceFlux .gt. 0.0) then
         U_up = U_old(1,iPoint)
         V_up = U_old(2,iPoint)
         if (muscl) then
             U_up = U_up + GradU(1,2,i,j)*dy/2.0
             V_up = V_up + GradU(2,2,i,j)*dy/2.0
         endif
       else
         U_up = U_old(1,jPoint)
         V_up = U_old(2,jPoint)
         if (muscl) then
             U_up = U_up - GradU(1,2,i,j+1)*dy/2.0
             V_up = V_up - GradU(2,2,i,j+1)*dy/2.0
         endif
       endif
       F_n(1) = FaceFlux*U_up ! (rho*V)*U|_n
       F_n(2) = FaceFlux*V_up ! (rho*V)*V|_n
       else
       F_n(1) = rho*V_n*U_n*dx          ! (rho*V)*U|_n
       F_n(2) = rho*V_n*V_n*dx          ! (rho*V)*V|_n
       endif
       lambda_i = abs(2.0*U_old(2,iPoint)*dx)
       lambda_j = abs(2.0*U_old(2,jPoint)*dx)
       
       lambda_mean = 0.5*(lambda_i + lambda_j)
       if (lambda_mean .eq. 0.0) then
        lambda_i = abs(U_old(2,iPoint)*dx) + sqrt(U_old(2,iPoint)*U_old(2,iPoint)*dx*dx + (artvisc/rho)*dx*dx) 
        lambda_j = abs(U_old(2,jPoint)*dx) + sqrt(U_old(2,jPoint)*U_old(2,jPoint)*dx*dx + (artvisc/rho)*dx*dx)
        lambda_mean = 0.5*(lambda_i + lambda_j)
       endif
       Phi_i = (lambda_i/(4.0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.0*lambda_mean))**Param_p
       if ((Phi_i + Phi_j).eq.0.0) then 
           SF = 0.0
       else 
            SF = 4.0*Phi_i*Phi_j/(Phi_i + Phi_j)
       endif
       sc0 = 3.0*(4.0+4.0)/(4.0*4.0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.0/3.0
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
       FaceFlux = rho*V_s*dx
       if ((i.eq.1).or.(i.eq.Nx)) FaceFlux = rho*V_s*dx/2.0
       if (FaceFlux .gt. 0.0) then
         U_up = U_old(1,jPoint)
         V_up = U_old(2,jPoint)
         if (muscl) then
             U_up = U_up + GradU(1,2,i,j-1)*dy/2.0
             V_up = V_up + GradU(2,2,i,j-1)*dy/2.0
         endif
       else
         U_up = U_old(1,iPoint)
         V_up = U_old(2,iPoint)
         if (muscl) then
             U_up = U_up - GradU(1,2,i,j)*dy/2.0
             V_up = V_up - GradU(2,2,i,j)*dy/2.0
         endif
       endif
       F_s(1) = -FaceFlux*U_up ! (rho*V)*U|_s
       F_s(2) = -FaceFlux*V_up ! (rho*V)*V|_s
       else
       F_s(1) = -rho*V_s*U_s*dx          ! (rho*V)*U|_s
       F_s(2) = -rho*V_s*V_s*dx          ! (rho*V)*V|_s
       endif
       lambda_i = abs(2.0*U_old(2,iPoint)*dx)
       lambda_j = abs(2.0*U_old(2,jPoint)*dx)
       
       lambda_mean = 0.5*(lambda_i + lambda_j)
       if (lambda_mean .eq. 0.0) then
        lambda_i = abs(U_old(2,iPoint)*dx) + sqrt(U_old(2,iPoint)*U_old(2,iPoint)*dx*dx + (artvisc/rho)*dx*dx) 
        lambda_j = abs(U_old(2,jPoint)*dx) + sqrt(U_old(2,jPoint)*U_old(2,jPoint)*dx*dx + (artvisc/rho)*dx*dx)
        lambda_mean = 0.5*(lambda_i + lambda_j)
       endif
       Phi_i = (lambda_i/(4.0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.0*lambda_mean))**Param_p
       if ((Phi_i + Phi_j).eq.0.0) then 
           SF = 0.0
       else 
            SF = 4.0*Phi_i*Phi_j/(Phi_i + Phi_j)
       endif
       sc0 = 3.0*(4.0+4.0)/(4.0*4.0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.0/3.0
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
       
       R(1:2,iPoint) = Fc(1:2,i,j)

     enddo
   enddo



end subroutine convective_residual
