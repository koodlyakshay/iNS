!> \file sub_variable_operations.f90
!! \brief File with various subroutines to handle allocation, initialization and deletion of variables.
subroutine allocate_vars

use global_vars

implicit none

integer     :: nPoint
real        :: CFL
real        :: Re

  nPoint = Nx*Ny


  allocate(x(Nx,Ny))
  allocate(y(Nx,Ny))
  allocate(P(Nx,Ny))
  allocate(Vol(nPoint))
  allocate(U(2,Nx*Ny))
  allocate(U_old(2,Nx*Ny))
  allocate(Sol(3,Nx*Ny))
  allocate(V(3,Nx,Ny))
  allocate(Fc(2,Nx,Ny))
  allocate(Fv(2,Nx,Ny))
  allocate(P_Correc(nPoint))
  allocate(Vel_Corr(2,Nx*Ny))
  allocate(R(3,Nx*Ny))
  allocate(Jac(nPoint,nPoint))
  allocate(GradU(3,2,Nx,Ny))
  allocate(Res(Nx*Ny*nVar))
  allocate(Mass(Nx*Ny),D(2,Nx*Ny))
  allocate(Tot_Jac(Nx*Ny*nVar,Nx*Ny*nVar))
  allocate(GradPc(2,Nx,Ny))

end subroutine allocate_vars

subroutine initialize_vars

use global_vars
use flow_parmaters

logical           :: old

  inquire(file = '../out/',exist=old)
  
  if (.not. old) then
    call system('mkdir ../out/')
  endif  
  
  do i=1,Nx
   do j=1,Ny
      x(i,j) = xmin + (i-1)*dx
      y(i,j) = ymin + (j-1)*dy
      if (x(i,j).lt.0.0) k = i
   enddo
  enddo

!First set everything to 0.5dx * 0.5 dy (the lowest)
  Vol(:) = 0.25*dx*dy
!Change the interior points so only boundary points have the old value
  do i=2,Nx-1
   do j=2,Ny-1
     iPoint = i + (j-1)*Nx
     Vol(iPoint) = dx*dy 
   enddo
  enddo

!Change boundary points not at the corner so only corner points have the lowest value 
!and others have the correct value
  i=1
  do j=2,Ny-1
   iPoint = i + (j-1)*Nx
   Vol(iPoint) = 0.5*dx*dy 
  enddo
  i=Nx
  do j=2,Ny-1
   iPoint = i + (j-1)*Nx
   Vol(iPoint) = 0.5*dx*dy 
  enddo
  j=1
  do i=2,Nx-1
   iPoint = i + (j-1)*Nx
   Vol(iPoint) = 0.5*dx*dy 
  enddo
  j=Ny
  do i=2,Nx-1
   iPoint = i + (j-1)*Nx
   Vol(iPoint) = 0.5*dx*dy 
  enddo

!--- Flow defintion ---!
  Re_l = 0.25
  P_outlet = 0.0
  Re = (Re_l*U_inf*rho)/mu
  
  artvisc = 4.0
  CFL = 0.1
  alfa = 0.9

!--- Initialize variables ---!

  do i=1,Nx
   do j=1,Ny
      V(1:3,i,j) = 0.0
      Fc(1:2,i,j) = 0.0
      Fv(1:2,i,j) = 0.0
      P(i,j) = 0.0
   enddo
  enddo

  do iPoint = 1,nPoint
   U(1,iPoint) = 0.0
   U_old(1,iPoint) = 0.0
   U(2,iPoint) = 0.0
   U_old(2,iPoint) = 0.0
  enddo
  
  if (restart) call read_restart
  
  open(unit=10,file='../out/Solver_details.txt',status='unknown')
  write(10,*)'CFL_m: ',U_inf*dt_m/min(dx,dy)
  write(10,*)'CFL_p: ',U_inf*dt_p/min(dx,dy)
  write(10,*)'Re: ',Re
  write(10,*)'dx: ',dx
  write(10,*)'dy: ',dy
  write(10,*)'rho: ',rho
  write(10,*)'xmin: ',xmin
  write(10,*)'xmax: ',x(Nx,Ny)
  write(10,*)'kappa: ',kappa
  write(10,*)'U_inf: ',U_inf
  write(10,*)'mu: ',mu
  write(10,*)'Re_l: ',Re_l
  write(10,*)'P_inf: ',P_inf
  close(10)


end subroutine initialize_vars


subroutine delete_vars

use global_vars

  deallocate(x)
  deallocate(y)
  deallocate(P)
  deallocate(Vol)
  deallocate(U)
  deallocate(U_old)
  deallocate(Sol)
  deallocate(V)
  deallocate(Fc)
  deallocate(Fv)
  deallocate(P_Correc)
  deallocate(Vel_Corr)
  deallocate(R)
  deallocate(Jac)
  deallocate(GradU)
  deallocate(Res)
  deallocate(Mass,D)
  deallocate(Tot_Jac)
  deallocate(GradPc)


end subroutine delete_vars
