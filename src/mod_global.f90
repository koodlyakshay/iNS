!> \file mod_global.f90
!! \brief File with global variables.
module global_vars

implicit none
integer           :: i,j,k
integer           :: Nx,Ny,iPoint,jpoint,nPoint,nDim,jDim,nVar,iVar
integer           :: p_screen1,p_screen2,p_out
integer           :: ExtIter,PIter,MIter,nExtIter,nPIter,nMIter,convergence,liniter,wrt_data
real              :: dx,dy,Area,Normal,Lx,Ly,xmin,ymin,Res_mass,alfa,kappa
real              :: dt,CFL,lambda_inv,lambda_visc,dt_i,dt_v,lambda_i,lambda_j,dt_m,dt_p
real              :: Param_p,sc0,SF,Phi_i,Phi_j,lambda_mean,E_0,mass_l2,Res_l2
real              :: Re, U_inf,rho,mu,Re_l,P_inf,GradPi,GradPav,artvisc,P_outlet
real              :: U_e,U_w,U_n,U_s,F_e(2),F_w(2),F_n(2),F_s(2),Gr_e,Gr_w,Gr_n,Gr_s
real              :: V_e,V_w,V_n,V_s,U_up,V_up,FaceFlux,mass_in,mass_out
real,allocatable  :: x(:,:),y(:,:),Mass(:),D(:,:),P_Correc(:),Vol(:), Vel_Corr(:,:)
real,allocatable  :: U(:,:),U_old(:,:),Sol(:,:),Fc(:,:,:),Fv(:,:,:)
real,allocatable  :: R(:,:),Jac(:,:),P(:,:),V(:,:,:),Res(:),Mat(:,:)
real,allocatable  :: GradU(:,:,:,:),Solp(:),Tot_Jac(:,:),Tot_R(:),Tot_Sol(:),GradPc(:,:,:)
character(len=100)  :: file_out
logical             :: implicit_time,upwind,muscl


end module global_vars
