!> \file sub_read_input.f90
!! \brief Subroutine to read input data from ../inp/grid.inp, ../inp/solver.inp, ../inp/numerics.inp and ../inp/flow.inp files.
!!
!--- Grid definition ---!
! *---*---*---*---*   j=Ny
! *---*---*---*---*    .
! *---*---*---*---*    .
! *---*---*---*---*    .
! *---*---*---*---*   j=1
!i=1,...........,i=Nx
subroutine read_input

use global_vars
use output_vars
use parameters

implicit none


namelist            /grid/ Lx,Ly,Nx,Ny,xmin,ymin
namelist            /solver/ implicit_time,upwind,muscl,nExtIter,nPIter,nMIter,p_screen1,p_screen2,p_out,file_out,restart,p_rest
namelist            /numerics/ kappa,Param_p,dt_m,dt_p
namelist            /flow/ U_inf,rho,mu,Re_l,P_inf

character(len=100)   :: filename

  CALL get_command_argument(1, filename)
  IF (len_TRIM(filename) == 0) then
    print*, "Please provide the input file name!"
    stop
  endif
 
  print*, "input file: ", trim(filename)

  open(unit=16, file=trim(filename), status='old')
  
  read(16, nml=grid)
  
  read(16, nml=solver)
  
  read(16, nml=numerics)
  
  read(16, nml=flow)
  
  
  dx = Lx/(Nx-1)
  dy = Ly/(Ny-1)
  nVar = 2
  
  close(16)
  
end subroutine read_input
