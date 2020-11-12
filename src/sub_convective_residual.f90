!> \file sub_convective_residual.f90
!! \brief Subroutine to compute residual contribution from convection terms of the momentum equations.

subroutine convective_residual(upwind, muscl,quick)

use global_vars
use flow_parmaters


implicit none

logical,intent(in) :: upwind, muscl, quick
integer            :: i,j,iPoint,jpoint, kPoint
real               :: FaceFlux
real               :: lambda_i,lambda_j,lambda_mean
real               :: E_0, sc0, SF
real               :: Phi_i, Phi_j
real               :: artvisc = 4.0
real               :: V_e, V_w, V_n, V_s
real               :: U_e, U_w, U_n, U_s
real               :: U_up, V_up
real               :: F_e(2), F_w(2), F_n(2), F_s(2)


E_0 = 0.0
lambda_mean = 0.0
SF = 0.0
do i=1,Nx
    do j=1,Ny
    !--- Node definition ---!
    iPoint = i + (j-1)*Nx
    F_e = 0.0 
    F_w = 0.0
    F_n = 0.0
    F_s = 0.0
    !East
    if (i.ne.Nx) then
      !--- Compute flux and update Jacobians ---!
      jPoint = i+1 + (j-1)*Nx
      U_e = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
      V_e = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
      if (upwind) then
        FaceFlux = rho*U_e*dy
        if ((j.eq.1).or.(j.eq.Ny)) FaceFlux = rho*U_e*dy/2.0

        if (FaceFlux .gt. 0.0) then
          if (quick) then
            if (i .eq. 1) then 
              call upwind_res(.true., FaceFlux, U_old(:,iPoint),GradU(:,1,i,j),dx/2.0,1,F_e)
            else 
              kPoint = i-1 + (j-1)*Nx
              call quick_res(FaceFlux,U_old(:,iPoint),U_old(:,jPoint),U_old(:,kPoint),1,F_e)
            endif
          else
            call upwind_res(muscl, FaceFlux, U_old(:,iPoint),GradU(:,1,i,j),dx/2.0,1,F_e)
          endif

        else
          if (quick) then
            if (i .eq. Nx-1) then
              call upwind_res(.true., FaceFlux, U_old(:,jPoint),GradU(:,1,i+1,j),-dx/2.0,1,F_e)
            else 
              kPoint = i+2 + (j-1)*Nx
              call quick_res(FaceFlux,U_old(:,jPoint),U_old(:,iPoint),U_old(:,kPoint),1,F_e)
            endif
          else
            call upwind_res(muscl, FaceFlux, U_old(:,jPoint),GradU(:,1,i+1,j),-dx/2.0,1,F_e)
          endif
        endif
      else
        call centered_res(U_old(:,iPoint),U_old(:,jPoint),dy,1,1,F_e)
      endif
              
       Convec_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = 1.0*U_e*dy - E_0*SF*lambda_mean
       Convec_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = 0.0
       Convec_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = 0.5*V_e*dy
       Convec_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = 0.5*U_e*dy - E_0*SF*lambda_mean
       
       Convec_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 1.0*U_e*dy + E_0*SF*lambda_mean
       Convec_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = 0.0
       Convec_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = 0.5*V_e*dy
       Convec_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*U_e*dy + E_0*SF*lambda_mean
    endif  

    !West
    if (i.ne.1) then
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       V_w = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       if (upwind) then
         FaceFlux = rho*U_w*dy
         if ((j.eq.1).or.(j.eq.Ny)) FaceFlux = rho*U_w*dy/2.0
         if (FaceFlux .gt. 0.0) then
           if (quick) then
             if (i.eq.2) then
               call upwind_res(.true., FaceFlux, U_old(:,jPoint),GradU(:,1,i-1,j),dx/2.0,-1,F_w)
             else
               kPoint = i-2 + (j-1)*Nx
               call quick_res(FaceFlux,U_old(:,jPoint),U_old(:,iPoint),U_old(:,kPoint),-1,F_w)
             endif
           else
             call upwind_res(muscl, FaceFlux, U_old(:,jPoint),GradU(:,1,i-1,j),dx/2.0,-1,F_w)
           endif
         else
           if (quick) then
             if(i .eq. Nx) then
               call upwind_res(muscl, FaceFlux, U_old(:,iPoint),GradU(:,1,i,j),-dx/2.0,-1,F_w)
             else
               kPoint = i+1 + (j-1)*Nx
               call quick_res(FaceFlux,U_old(:,iPoint),U_old(:,jPoint),U_old(:,kPoint),-1,F_w)
             endif
           else
             call upwind_res(muscl, FaceFlux, U_old(:,iPoint),GradU(:,1,i,j),-dx/2.0,-1,F_w)
           endif
         endif
       else
         call centered_res(U_old(:,iPoint),U_old(:,jPoint),dy,1,-1,F_w)
       endif
       
       Convec_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = -1.0*U_w*dy - E_0*SF*lambda_mean
       Convec_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = -0.0
       Convec_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = -0.5*V_w*dy
       Convec_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = -0.5*U_w*dy - E_0*SF*lambda_mean
       
       Convec_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = -1.0*U_w*dy + E_0*SF*lambda_mean
       Convec_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = -0.0
       Convec_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = -0.5*V_w*dy
       Convec_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = -0.5*U_w*dy + E_0*SF*lambda_mean
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
           if (quick) then
             if (j.eq.1) then 
               call upwind_res(.true., FaceFlux, U_old(:,iPoint),GradU(:,2,i,j),dy/2.0,1,F_n)
             else
               kPoint = i + (j-1-1)*Nx
               call quick_res(FaceFlux,U_old(:,iPoint),U_old(:,jPoint),U_old(:,kPoint),1,F_n)
             endif
           else
             call upwind_res(muscl, FaceFlux, U_old(:,iPoint),GradU(:,2,i,j),dy/2.0,1,F_n)
           endif
         else
           if (quick) then
             if (j.eq.Ny-1) then
               call upwind_res(.true., FaceFlux, U_old(:,jPoint),GradU(:,2,i,j+1),-dy/2.0,1,F_n)
             else
               kPoint = i + (j-1+2)*Nx
               call quick_res(FaceFlux,U_old(:,jPoint),U_old(:,iPoint),U_old(:,kPoint),1,F_n)
             endif
           else
             call upwind_res(muscl, FaceFlux, U_old(:,jPoint),GradU(:,2,i,j+1),-dy/2.0,1,F_n)
           endif
         endif
       else
         call centered_res(U_old(:,iPoint),U_old(:,jPoint),dx,2,1,F_n)
       endif
              
       Convec_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = 0.5*V_n*dx - E_0*SF*lambda_mean
       Convec_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = 0.5*U_n*dx
       Convec_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = 0.0
       Convec_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = 1.0*V_n*dx - E_0*SF*lambda_mean
       
       Convec_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*V_n*dx + E_0*SF*lambda_mean
       Convec_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = 0.5*U_n*dx
       Convec_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = 0.0
       Convec_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 1.0*V_n*dx + E_0*SF*lambda_mean
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
           if (quick) then
             if (j.eq.2) then
               call upwind_res(.true., FaceFlux, U_old(:,jPoint),GradU(:,2,i,j-1),dy/2.0,-1,F_s)
             else
               kPoint = i + (j-1-2)*Nx
               call quick_res(FaceFlux,U_old(:,jPoint),U_old(:,iPoint),U_old(:,kPoint),-1,F_s)
             endif
           else
             call upwind_res(muscl, FaceFlux, U_old(:,jPoint),GradU(:,2,i,j-1),dy/2.0,-1,F_s)
           endif
         else
           if (quick) then
             if (j.eq.Ny) then 
               call upwind_res(.true., FaceFlux, U_old(:,iPoint),GradU(:,2,i,j),-dy/2.0,-1,F_s)
             else  
               kPoint = i + (j-1+1)*Nx
               call quick_res(FaceFlux, U_old(:,iPoint), U_old(:,jPoint),U_old(:,kPoint),-1,F_s)
             endif
           else
             call upwind_res(muscl, FaceFlux, U_old(:,iPoint),GradU(:,2,i,j),-dy/2.0,-1,F_s)
           endif
         endif
       else
         call centered_res(U_old(:,iPoint),U_old(:,jPoint),dx,2,-1,F_s)
       endif
       
       Convec_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) =- 0.5*V_s*dx - E_0*SF*lambda_mean
       Convec_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) =- 0.5*U_s*dx
       Convec_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) =- 0.0
       Convec_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) =- 1.0*V_s*dx - E_0*SF*lambda_mean
       
       Convec_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) =- 0.5*V_s*dx + E_0*SF*lambda_mean
       Convec_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) =- 0.5*U_s*dx
       Convec_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) =- 0.0
       Convec_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) =- 1.0*V_s*dx + E_0*SF*lambda_mean
    endif    
       !--- Update residuals ---!
       
       Fc(1,i,j) = F_e(1) + F_w(1) + F_n(1) + F_s(1) ! (rho*U)*U|_e - (rho*U)*U|_w + (rho*U)*V|_n - (rho*U)*V|_s
       Fc(2,i,j) = F_e(2) + F_w(2) + F_n(2) + F_s(2) ! (rho*V)*U|_e - (rho*V)*U|_w + (rho*V)*V|_n - (rho*V)*V|_s
       
       R(1:2,iPoint) = Fc(1:2,i,j)
     enddo
   enddo



end subroutine convective_residual

subroutine centered_res(U_i,U_j,area,nrml,dir,cenRes)

! Subroutine to calculate centered residual

use global_vars, only: rho

implicit none

integer,intent(in):: nrml,dir
real,intent(in)   :: U_i(2),U_j(2),area
real,intent(out)  :: cenRes(2)
real              :: diss(2),FF,U_f

U_f = 0.5*(U_i(nrml) + U_j(nrml))
FF = rho*U_f*area

!Get artificial dissipation
call getdissipation(U_i,U_j,area,nrml,diss)

!Calculate Residual
cenRes(:) = FF*0.5*(U_i(:) + U_j(:))*dir + diss(:)


end subroutine centered_res

subroutine getdissipation(U_i, U_j, area,nrml,diss)

! Subroutine to get artificial dissipation for central scheme

use flow_parmaters
use global_vars, only: rho

implicit none

real,intent(in)   :: U_i(2),U_j(2),area
real,intent(out)  :: diss(2)
real              :: artvisc = 4.0,area2
real              :: lambda_i,lambda_j, lambda_mean
real              :: Phi_i, Phi_j
real              :: sc0, SF, E_0
integer,intent(in):: nrml

area2 = area*area

!Dissipation
lambda_i = abs(2.0*U_i(nrml)*area)
lambda_j = abs(2.0*U_j(nrml)*area)
       
lambda_mean = 0.5*(lambda_i + lambda_j)
if (lambda_mean .eq. 0.0) then
  lambda_i = abs(U_i(nrml)*area) + sqrt(U_i(nrml)*U_i(nrml)*area2 + (artvisc/rho)*area2) 
  lambda_j = abs(U_j(nrml)*area) + sqrt(U_j(nrml)*U_j(nrml)*area2 + (artvisc/rho)*area2) 
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

diss(1) = E_0*(U_i(1) - U_j(1))*SF*lambda_mean
diss(2) = E_0*(U_i(2) - U_j(2))*SF*lambda_mean

end subroutine getdissipation

subroutine upwind_res(muscl, FF, U_up, GradU_up, dist, dir,UpRes)

!Subroutine to compute upwind residual either with first order scheme or muscl reconstruction

implicit none

real,intent(in)    :: FF,U_up(2), GradU_up(2), dist
logical,intent(in) :: muscl
real,intent(out)   :: UpRes(2)
integer,intent(in) :: dir
real               :: U_recon(2)

U_recon = U_up

! Reconstruct if necessary
if (muscl) call recon(U_up,GradU_up, dist, U_recon)

! Calculate the flux
UpRes(:) = FF*U_recon(:)*dir

end subroutine upwind_res

subroutine recon(U_var,GradUvar, dist, U_recon)

!Subroutine to reconstruct variable

implicit none

real,intent(in)  :: U_var(2), GradUvar(2),dist
real,intent(out) :: U_recon(2)

U_recon(:) = U_var(:) + GradUvar(:)*dist

end subroutine recon

subroutine quick_res(FF, U_C, U_E, U_UP1, dir, UpRes)

! Subroutine to calculate residual using QUICK scheme

implicit none
real,intent(in)     :: FF,U_C(2),U_E(2),U_UP1(2)
real,intent(out)    :: UpRes(2)
integer,intent(in)  :: dir
real                :: coeff_c, coeff_e, coeff_up1,U_recon(2)

!Assuming equal spacing in each dir
coeff_c = 6.0/8.0
coeff_e = 3.0/8.0
coeff_up1 = -1.0/8.0

! Face Velocity
U_recon(:) = coeff_c*U_C(:) + coeff_e*U_E(:) + coeff_up1*U_UP1(:)

! Residual
UpRes(:) = FF*U_recon(:)*dir

end subroutine quick_res
