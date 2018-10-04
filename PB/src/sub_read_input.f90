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

read(16,grid)

read(17,solver)

read(18,numerics)


nPoint = Nx*Ny
nDim = 2
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)
Vol = dx*dy
nVar = nDim

close(16)
close(17)
close(18)
end subroutine read_input
