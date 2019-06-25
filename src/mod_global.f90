!> \file mod_global.f90
!! \brief File with global variables.
module global_vars

implicit none

!integer           :: i,j,k,iPoint,jpoint

integer           :: Nx
integer           :: Ny
integer           :: nvar
real              :: dx
real              :: dy
real              :: Area
real              :: Lx
real              :: Ly
real              :: sc0
real              :: SF


logical           :: restart
real              :: dt_m
real              :: dt_p
real              :: Param_p
real              :: kappa
real              :: xmin
real              :: ymin
real              :: alfa


real              :: U_inf
real              :: rho
real              :: mu
real              :: Re_l
real              :: P_inf
real              :: artvisc
real              :: P_outlet
real              :: U_e
real              :: U_w
real              :: U_n
real              :: U_s
real              :: F_e(2)
real              :: F_w(2)
real              :: F_n(2)
real              :: F_s(2)
real              :: Gr_e
real              :: Gr_w
real              :: Gr_n
real              :: Gr_s
real              :: V_e
real              :: V_w
real              :: V_n
real              :: V_s
real              :: U_up
real              :: V_up
real,allocatable  :: x(:,:)
real,allocatable  :: y(:,:)
real,allocatable  :: Mass(:)
real,allocatable  :: D(:,:)
real,allocatable  :: P_Correc(:)
real,allocatable  :: Vol(:)
real,allocatable  :: Vel_Corr(:,:)
real,allocatable  :: U(:,:)
real,allocatable  :: U_old(:,:)
real,allocatable  :: Sol(:,:)
real,allocatable  :: Fc(:,:,:)
real,allocatable  :: Fv(:,:,:)
real,allocatable  :: R(:,:)
real,allocatable  :: Jac(:,:)
real,allocatable  :: P(:,:)
real,allocatable  :: V(:,:,:)
real,allocatable  :: Res(:)
real,allocatable  :: Mat(:,:)
real,allocatable  :: GradU(:,:,:,:)
real,allocatable  :: Solp(:)
real,allocatable  :: Tot_Jac(:,:)
real,allocatable  :: GradPc(:,:,:)

end module global_vars
