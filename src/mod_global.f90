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
real              :: U_inf
real              :: rho
real              :: mu
real              :: Re_l
real              :: P_inf
real              :: P_outlet




real,allocatable  :: x(:,:)
real,allocatable  :: y(:,:)
real,allocatable  :: Mass(:)
real,allocatable  :: D(:,:)
real,allocatable  :: P_Correc(:)
real,allocatable  :: Vol(:)
real,allocatable  :: dt_m(:)
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
real,allocatable  :: Convec_Jac(:,:)
real,allocatable  :: Visc_Jac(:,:)
real,allocatable  :: GradPc(:,:,:)

end module global_vars
