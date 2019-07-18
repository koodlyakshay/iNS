subroutine compute_gradient

use global_vars

!--- Compute gradients at nodes ---!
   do i=2,Nx-1
    do j=1,Ny
       !--- (var,dim,i,j) ---!
       GradU(1,1,i,j) = (V(1,i+1,j) - V(1,i-1,j))/(2.0*dx)
       GradU(2,1,i,j) = (V(2,i+1,j) - V(2,i-1,j))/(2.0*dx)
       GradU(3,1,i,j) = (V(3,i+1,j) - V(3,i-1,j))/(2.0*dx)
    enddo
   enddo
   
   do i=1,Nx
    do j=2,Ny-1
       !--- (var,dim,i,j) ---!
       GradU(1,2,i,j) = (V(1,i,j+1) - V(1,i,j-1))/(2.0*dy)
       GradU(2,2,i,j) = (V(2,i,j+1) - V(2,i,j-1))/(2.0*dy)
       GradU(3,2,i,j) = (V(3,i,j+1) - V(3,i,j-1))/(2.0*dy)
    enddo
   enddo
   
   i=1
   do j=2,Ny-1
     !--- (var,dim,i,j) ---!
     GradU(1,1,i,j) = (V(1,i+1,j) - V(1,i,j))/(dx)
     GradU(2,1,i,j) = (V(2,i+1,j) - V(2,i,j))/(dx)
     GradU(3,1,i,j) = (V(3,i+1,j) - V(3,i,j))/(dx)
     !print*,i + (j-1)*Nx,'BC left d/dx'
   enddo
   
   i=Nx
   do j=2,Ny-1
     !--- (var,dim,i,j) ---!
     GradU(1,1,i,j) = (V(1,i,j) - V(1,i-1,j))/(dx)
     GradU(2,1,i,j) = (V(2,i,j) - V(2,i-1,j))/(dx)
     GradU(3,1,i,j) = (V(3,i,j) - V(3,i-1,j))/(dx)
     !print*,i + (j-1)*Nx,'BC right d/dx'
   enddo
   
   j=Ny
   do i=2,Nx-1
     !--- (var,dim,i,j) ---!
     GradU(1,2,i,j) = (V(1,i,j) - V(1,i,j-1))/(dy)
     GradU(2,2,i,j) = (V(2,i,j) - V(2,i,j-1))/(dy)
     GradU(3,2,i,j) = (V(3,i,j) - V(3,i,j-1))/(dy)
     !print*,i + (j-1)*Nx,'BC top d/dy'
   enddo
   
   j=1
   do i=2,Nx-1
     !--- (var,dim,i,j) ---!
     GradU(1,2,i,j) = (V(1,i,j+1) - V(1,i,j))/(dy)
     GradU(2,2,i,j) = (V(2,i,j+1) - V(2,i,j))/(dy)
     GradU(3,2,i,j) = (V(3,i,j+1) - V(3,i,j))/(dy)
     !print*,i + (j-1)*Nx,'BC bottom d/dy'
   enddo
   
   i=Nx
   j=Ny
   !--- (var,dim,i,j) ---!
   GradU(1,1,i,j) = (V(1,i,j) - V(1,i-1,j))/dx
   GradU(1,2,i,j) = (V(1,i,j) - V(1,i,j-1))/dy
   GradU(2,1,i,j) = (V(2,i,j) - V(2,i-1,j))/dx
   GradU(2,2,i,j) = (V(2,i,j) - V(2,i,j-1))/dy
   GradU(3,1,i,j) = (V(3,i,j) - V(3,i-1,j))/dx
   GradU(3,2,i,j) = (V(3,i,j) - V(3,i,j-1))/dy
   
   i=1
   j=1
   !--- (var,dim,i,j) ---!
   GradU(1,1,i,j) = (V(1,i+1,j) - V(1,i,j))/dx
   GradU(1,2,i,j) = (V(1,i,j+1) - V(1,i,j))/dy
   GradU(2,1,i,j) = (V(2,i+1,j) - V(2,i,j))/dx
   GradU(2,2,i,j) = (V(2,i,j+1) - V(2,i,j))/dy
   GradU(3,1,i,j) = (V(3,i+1,j) - V(3,i,j))/dx
   GradU(3,2,i,j) = (V(3,i,j+1) - V(3,i,j))/dy

   i=1
   j=Ny
   !--- (var,dim,i,j) ---!
   GradU(1,1,i,j) = (V(1,i+1,j) - V(1,i,j))/dx
   GradU(1,2,i,j) = (V(1,i,j) - V(1,i,j-1))/dy
   GradU(2,1,i,j) = (V(2,i+1,j) - V(2,i,j))/dx
   GradU(2,2,i,j) = (V(2,i,j) - V(2,i,j-1))/dy
   GradU(3,1,i,j) = (V(3,i+1,j) - V(3,i,j))/dx
   GradU(3,2,i,j) = (V(3,i,j) - V(3,i,j-1))/dy
      
   i=Nx
   j=1
   !--- (var,dim,i,j) ---!
   GradU(1,1,i,j) = (V(1,i,j) - V(1,i-1,j))/dx
   GradU(1,2,i,j) = (V(1,i,j+1) - V(1,i,j))/dy
   GradU(2,1,i,j) = (V(2,i,j) - V(2,i-1,j))/dx
   GradU(2,2,i,j) = (V(2,i,j+1) - V(2,i,j))/dy
   GradU(3,1,i,j) = (V(3,i,j) - V(3,i-1,j))/dx
   GradU(3,2,i,j) = (V(3,i,j+1) - V(3,i,j))/dy

end subroutine compute_gradient
