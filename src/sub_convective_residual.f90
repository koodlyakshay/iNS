subroutine convective_residual

use global_vars

do i=1,Nx
    do j=1,Ny
    !--- Node definition ---!
    iPoint = i + (j-1)*Nx
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
       
       R(1:2,iPoint) = Fc(1:2,i,j)
       
       if (wrt_data .eq. 1) then
       write(16,*)F_e(1:2),'iPoint',iPoint,'jPoint',i+1 + (j-1)*Nx,U_e,V_e
       write(16,*)F_w(1:2),'iPoint',iPoint,'jPoint',i-1 + (j-1)*Nx,U_w,V_w
       write(16,*)F_n(1:2),'iPoint',iPoint,'jPoint',i + (j+1-1)*Nx,U_n,V_n
       write(16,*)F_s(1:2),'iPoint',iPoint,'jPoint',i + (j-1-1)*Nx,U_s,V_s
       write(16,*)
       endif
     enddo
   enddo



end subroutine convective_residual