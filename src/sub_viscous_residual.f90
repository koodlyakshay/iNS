!> \file sub_viscous_residual.f90
!! \brief Subroutine to compute the viscous contribution to the residual for the momentum equations.
subroutine viscous_residual

use global_vars


  do i=1,Nx
     do j=1,Ny
    iPoint = i + (j-1)*Nx
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
    
     !--- Update residuals ---!
      Fv(1,i,j) = F_e(1) + F_w(1) + F_n(1) + F_s(1) ! du/dx|_e - du/dx|_w + du/dy|_n - du/dy|_s
      Fv(2,i,j) = F_e(2) + F_w(2) + F_n(2) + F_s(2) ! dv/dx|_e - dv/dx|_w + dv/dy|_n - dv/dy|_s
       
      R(1:2,iPoint) = R(1:2,iPoint) - Fv(1:2,i,j)
       
    enddo
   enddo


end subroutine viscous_residual
