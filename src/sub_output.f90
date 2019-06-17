!> \file sub_output.f90
!! \brief Subroutine to write solution output in .vtk format.
subroutine sub_output

use global_vars
use output_vars
 
real            :: dummy
character         :: fmt1(30),fmt2



open(unit=100,file='../out/'//trim(file_out)//'.vtk',status='unknown')
dummy = 0.0

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
  write(100,'(F9.6,F9.6,F9.6)') x(i,j), y(i,j), 0.0
 enddo
enddo

write(100,'(A10,A,I5)') 'POINT_DATA',achar(9), Nx*Ny

write(100,'(a)') 'SCALARS Pressure double 1'
write(100,'(a)') 'LOOKUP_TABLE default'

do j=1,Ny
 do i=1,Nx
  iPoint = i + (j-1)*Nx
  write(100,'(F10.6)') P(i,j)
 enddo
enddo

write(100,'(a)') 'SCALARS Mass double 1'
write(100,'(a)') 'LOOKUP_TABLE default'

do j=1,Ny
 do i=1,Nx
  iPoint = i + (j-1)*Nx
  write(100,'(F10.6)') Mass(iPoint)
 enddo
enddo

write(100,'(a)') 'VECTORS Velocity double'

do j=1,Ny
 do i=1,Nx
  iPoint = i + (j-1)*Nx
  write(100,'(F10.6,A,F10.6,A,F10.6)') U(1,iPoint),achar(9),U(2,iPoint),achar(9),dummy
 enddo
enddo


write(100,'(a)') 'VECTORS VelocityCorrection double'

do j=1,Ny
 do i=1,Nx
  iPoint = i + (j-1)*Nx
  write(100,'(F10.6,A,F10.6,A,F10.6)') Vel_Corr(1,iPoint),achar(9),Vel_Corr(2,iPoint),achar(9),dummy
 enddo
enddo
close(100)

end subroutine sub_output
