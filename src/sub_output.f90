subroutine sub_output

use global_vars

real              :: tmp(165)
character         :: fmt1(30),fmt2

open(unit=100,file='../out/test/flow.vtk',status='unknown')
tmp = 1.0

!--- These three lines are compulsory. ---!
write(100,'(a)')'# vtk DataFile Version 3.0'   ! File version and identifier
write(100,'(a)')'vtk output'                   ! Header
write(100,'(a)')'ASCII'                        ! Type of data, only ascii now

!--- Write grid information. ---!
write(100,'(a)')'DATASET STRUCTURED_GRID'      ! Dataset is a keyword
                                               ! We deal with only structured grids

write(100,'(A10,A,I3,A,I3,A,I3)')'DIMENSIONS',achar(9),Nx,achar(9),Ny,achar(9),1
write(100,'(A6,A,I5,A,A6)')'POINTS',achar(9), Nx*Ny,achar(9), 'float'

do j=1,Ny
 do i=1,Nx
  write(100,'(F8.6,F8.6,F8.6)') x(i,j), y(i,j), 0.0
 enddo
enddo

write(100,'(A10,A,I5)') 'POINT_DATA',achar(9), Nx*Ny

write(100,'(a)') 'SCALARS Pressure double 1'
write(100,'(a)') 'LOOKUP_TABLE default'

do j=1,Ny
 do i=1,Nx
  iPoint = i + (j-1)*Nx
  write(100,'(F6.2)') P(i,j)
 enddo
enddo
close(100)

end subroutine sub_output
