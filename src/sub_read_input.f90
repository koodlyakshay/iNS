!> \file sub_read_input.f90
!! \brief Subroutine to read input data from ../inp/grid.inp, ../inp/solver.inp, ../inp/numerics.inp and ../inp/flow.inp files.
subroutine read_input

use global_vars

!--- Grid definition ---!
! *---*---*---*---*   j=Ny
! *---*---*---*---*    .
! *---*---*---*---*    .
! *---*---*---*---*    .
! *---*---*---*---*   j=1
!i=1,...........,i=Nx
open(unit=16,file='../inp/grid.inp',status='unknown')
open(unit=17,file='../inp/solver.inp',status='unknown')
open(unit=18,file='../inp/numerics.inp',status='unknown')
open(unit=19,file='../inp/flow.inp',status='unknown')

read(16,grid)

read(17,solver)

read(18,numerics)

read(19,flow)


nPoint = Nx*Ny
nDim = 2
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)
nVar = nDim

close(16)
close(17)
close(18)
close(19)
end subroutine read_input
