module global_vars

implicit none
integer             :: i,j,k,Nx,Ny,iPoint,jpoint,nPoint,nDim,jDim,nVar,iVar,p_screen1,p_screen2,p_out
integer             :: ExtIter,PIter,MIter,nExtIter,nPIter,nMIter,convergence,liniter,wrt_data
real*8              :: dx,dy,vol,area,Lx,Ly,xmin,ymin,Res_mass,alfa,kappa
real*8              :: dt,CFL,lambda_inv,lambda_visc,dt_i,dt_v,lambda_i,lambda_j,dt_m,dt_p
real*8              :: Param_p,sc0,SF,Phi_i,Phi_j,lambda_mean,E_0,mass_l2,Res_l2
real*8              :: Re, U_inf,rho,mu,Re_l,P_inf,GradPi,GradPav,artvisc,P_outlet
real*8              :: U_e,U_w,U_n,U_s,F_e(2),F_w(2),F_n(2),F_s(2),Gr_e,Gr_w,Gr_n,Gr_s
real*8              :: V_e,V_w,V_n,V_s,U_up,V_up,FaceFlux
real*8,allocatable  :: x(:,:),y(:,:),Mass(:),D(:,:),P_Correc(:)
real*8,allocatable  :: U(:,:),U_old(:,:),Sol(:,:),Fc(:,:,:),Fv(:,:,:)
real*8,allocatable  :: R(:,:),Jac(:,:),P(:,:),V(:,:,:),Res(:),Mat(:,:)
real*8,allocatable  :: GradU(:,:,:,:),Solp(:),Tot_Jac(:,:),Tot_R(:),Tot_Sol(:)
logical             :: implicit_time,upwind,muscl
namelist            /grid/ Lx,Ly,Nx,Ny,xmin,ymin
namelist            /solver/ wrt_data,implicit_time,upwind,muscl,nExtIter,nPIter,nMIter,p_screen1,p_screen2,p_out
namelist            /numerics/ kappa,Param_p,dt_m,dt_p

end module global_vars
