subroutine compute_pcgradient

use global_vars

!--- Compute gradients at nodes ---!
   do i=2,Nx-1
    do j=1,Ny
       !--- (var,dim,i,j) ---!
       GradPc(1,i,j) = (P_correc(i+1 + (j-1)*Nx) - P_correc(i-1 + (j-1)*Nx))/(2.0*dx)
    enddo
   enddo
   
   do i=1,Nx
    do j=2,Ny-1
       !--- (var,dim,i,j) ---!
       GradPc(2,i,j) = (P_correc(i + (j+1-1)*Nx) - P_correc(i + (j-1-1)*Nx))/(2.0*dy)
    enddo
   enddo
   
   i=1
   do j=2,Ny-1
     !--- (var,dim,i,j) ---!
!     GradPc(1,i,j) = (P_correc(i+1 + (j-1)*Nx) - P_correc(i + (j-1)*Nx))/(dx)
     GradPc(1,i,j) = (-1.5*P_Correc(i + (j-1)*Nx) + 2.0*P_Correc(i+1 + (j-1)*Nx) - 0.5*P_Correc(i+2 + (j-1)*Nx))/dx
   enddo
   
   i=Nx
   do j=2,Ny-1
     !--- (var,dim,i,j) ---!
!     GradPc(1,i,j) = (P_correc(i + (j-1)*Nx) - P_correc(i-1 + (j-1)*Nx))/(dx)
     GradPc(1,i,j) = (1.5*P_Correc(i + (j-1)*Nx) - 2.0*P_Correc(i-1 + (j-1)*Nx) + 0.5*P_Correc(i-2 + (j-1)*Nx))/dx
   enddo
   
   j=Ny
   do i=2,Nx-1
     !--- (var,dim,i,j) ---!
!     GradPc(2,i,j) = (P_correc(i + (j-1)*Nx) - P_correc(i + (j-1-1)*Nx))/(dy)
     GradPc(2,i,j) = (1.5*P_Correc(i + (j-1)*Nx) - 2.0*P_Correc(i + (j-1-1)*Nx) + 0.5*P_Correc(i + (j-1-2)*Nx))/dy
   enddo
   
   j=1
   do i=2,Nx-1
     !--- (var,dim,i,j) ---!
!     GradPc(2,i,j) = (P_correc(i + (j+1-1)*Nx) - P_correc(i + (j-1)*Nx))/(dy)
     GradPc(2,i,j) = (-1.5*P_Correc(i + (j-1)*Nx) + 2.0*P_Correc(i + (j-1+1)*Nx) - 0.5*P_Correc(i + (j-1+2)*Nx))/dy
   enddo
   
   i=Nx
   j=Ny
   !--- (var,dim,i,j) ---!
!   GradPc(1,i,j) = (P_correc(i + (j-1)*Nx) - P_correc(i-1 + (j-1)*Nx))/dx
   GradPc(1,i,j) = (1.5*P_Correc(i + (j-1)*Nx) - 2.0*P_Correc(i-1 + (j-1)*Nx) + 0.5*P_Correc(i-2 + (j-1)*Nx))/dx
!   GradPc(2,i,j) = (P_correc(i + (j-1)*Nx) - P_correc(i + (j-1-1)*Nx))/dy
   GradPc(2,i,j) = (1.5*P_Correc(i + (j-1)*Nx) - 2.0*P_Correc(i + (j-1-1)*Nx) + 0.5*P_Correc(i + (j-1-2)*Nx))/dy
   
   i=1
   j=1
   !--- (var,dim,i,j) ---!
!   GradPc(1,i,j) = (P_correc(i+1 + (j-1)*Nx) - P_correc(i + (j-1)*Nx))/dx
   GradPc(1,i,j) = (-1.5*P_Correc(i + (j-1)*Nx) + 2.0*P_Correc(i+1 + (j-1)*Nx) - 0.5*P_Correc(i+2 + (j-1)*Nx))/dx
!   GradPc(2,i,j) = (P_correc(i + (j+1-1)*Nx) - P_correc(i + (j-1)*Nx))/dy
   GradPc(2,i,j) = (-1.5*P_Correc(i + (j-1)*Nx) + 2.0*P_Correc(i + (j-1+1)*Nx) - 0.5*P_Correc(i + (j-1+2)*Nx))/dy

   i=1
   j=Ny
   !--- (var,dim,i,j) ---!
!   GradPc(1,i,j) = (P_correc(i+1 + (j-1)*Nx) - P_correc(i + (j-1)*Nx))/dx
   GradPc(1,i,j) = (-1.5*P_Correc(i + (j-1)*Nx) + 2.0*P_Correc(i+1 + (j-1)*Nx) - 0.5*P_Correc(i+2 + (j-1)*Nx))/dx
!   GradPc(2,i,j) = (P_correc(i + (j-1)*Nx) - P_correc(i + (j-1-1)*Nx))/dy
   GradPc(2,i,j) = (1.5*P_Correc(i + (j-1)*Nx) - 2.0*P_Correc(i + (j-1-1)*Nx) + 0.5*P_Correc(i + (j-1-2)*Nx))/dy
      
   i=Nx
   j=1
   !--- (var,dim,i,j) ---!
!   GradPc(1,i,j) = (P_correc(i + (j-1)*Nx) - P_correc(i-1 + (j-1)*Nx))/dx
   GradPc(1,i,j) = (1.5*P_Correc(i + (j-1)*Nx) - 2.0*P_Correc(i-1 + (j-1)*Nx) + 0.5*P_Correc(i-2 + (j-1)*Nx))/dx
!   GradPc(2,i,j) = (P_correc(i + (j+1-1)*Nx) - P_correc(i + (j-1)*Nx))/dy
   GradPc(2,i,j) = (-1.5*P_Correc(i + (j-1)*Nx) + 2.0*P_Correc(i + (j-1+1)*Nx) - 0.5*P_Correc(i + (j-1+2)*Nx))/dy

end subroutine compute_pcgradient

subroutine compute_pcgradientgg

use global_vars
real*8 part_grad


do i=1,Nx
 do j=1,Ny
 
 GradPc(1:2,i,j) = 0.0
  iPoint = i + (j-1)*Nx

  !East
  if (i .ne. Nx) then
    jPoint = i+1 + (j-1)*Nx
    if ((j.eq.1).or.(j.eq.Ny)) then
       part_grad = 5.d-1*(P_Correc(iPoint) + P_Correc(jPoint))*dy/2.d0
    else 
       part_grad = 5.d-1*(P_Correc(iPoint) + P_Correc(jPoint))*dy
    endif
    GradPc(1,i,j) = GradPc(1,i,j) + part_grad
  endif
  !West
  if (i .ne. 1) then
    jPoint = i-1 + (j-1)*Nx
     if ((j.eq.1).or.(j.eq.Ny)) then
       part_grad = -5.d-1*(P_Correc(iPoint) + P_Correc(jPoint))*dy/2.d0
    else 
       part_grad = -5.d-1*(P_Correc(iPoint) + P_Correc(jPoint))*dy
    endif
    GradPc(1,i,j) = GradPc(1,i,j) + part_grad
  endif
  !North
  if (j.ne.Ny) then
    jPoint = i + (j+1-1)*Nx
    if ((i.eq.1).or.(i.eq.Nx)) then
       part_grad = 5.d-1*(P_Correc(iPoint) + P_Correc(jPoint))*dx/2.d0
    else
       part_grad = 5.d-1*(P_Correc(iPoint) + P_Correc(jPoint))*dx
    endif
    GradPc(2,i,j) = GradPc(2,i,j) + part_grad
  endif
  !South
  if (j.ne.1) then
    jPoint = i + (j-1-1)*Nx
    if ((i.eq.1).or.(i.eq.Nx)) then
       part_grad = -5.d-1*(P_Correc(iPoint) + P_Correc(jPoint))*dx/2.d0
    else
       part_grad = -5.d-1*(P_Correc(iPoint) + P_Correc(jPoint))*dx
    endif
    GradPc(2,i,j) = GradPc(2,i,j) + part_grad
  endif
 enddo
enddo

!--- Boundary elements ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    part_grad = -(P_Correc(iPoint))*dx
    GradPc(2,i,j) = GradPc(2,i,j) + part_grad
   enddo
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    part_grad = (P_Correc(iPoint))*dx
    GradPc(2,i,j) = GradPc(2,i,j) + part_grad
   enddo
   
   !--- Left inlet (i=1) ---!
   i=1
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    part_grad = -(P_Correc(iPoint))*dy
    GradPc(1,i,j) = GradPc(1,i,j) + part_grad
   enddo
   
   
   !--- Right outlet (i=Nx) ---!
   i=Nx
   do j=2,Ny-1
    iPoint = i + (j-1)*Nx
    part_grad = (P_Correc(iPoint))*dy
    GradPc(1,i,j) = GradPc(1,i,j) + part_grad
   enddo
   
   
   i=1
   j=1
   iPoint = i + (j-1)*Nx
 
   !Inlet face(West)
   part_grad = -(P_Correc(iPoint))*dy/2.d0
   GradPc(1,i,j) = GradPc(1,i,j) + part_grad
   
   !South
   part_grad = -(P_Correc(iPoint))*dx/2.d0
   GradPc(2,i,j) = GradPc(2,i,j) + part_grad
   
   

   i=1
   j=Ny
   iPoint = i + (j-1)*Nx

   !Inlet face(West)
   part_grad = -(P_Correc(iPoint))*dy/2.d0
   GradPc(1,i,j) = GradPc(1,i,j) + part_grad
       
   !North
   part_grad = (P_Correc(iPoint))*dx/2.d0
   GradPc(2,i,j) = GradPc(2,i,j) + part_grad
   
   
   i=Nx
   j=Ny
   iPoint = i + (j-1)*Nx

   !East (outlet)
   part_grad = (P_Correc(iPoint))*dy/2.d0
   GradPc(1,i,j) = GradPc(1,i,j) + part_grad
   
   !North
   part_grad = (P_Correc(iPoint))*dx/2.d0
   GradPc(2,i,j) = GradPc(2,i,j) + part_grad
   
    
   i=Nx
   j=1
   iPoint = i + (j-1)*Nx

   !East
   part_grad = (P_Correc(iPoint))*dy/2.d0
   GradPc(1,i,j) = GradPc(1,i,j) + part_grad
   
   !South
   part_grad = -(P_Correc(iPoint))*dx/2.d0
   GradPc(2,i,j) = GradPc(2,i,j) + part_grad
   
   
   do i=1,Nx
    do j=1,Ny
     iPoint = i + (j-1)*Nx
     GradPc(1,i,j) = GradPc(1,i,j)/Vol(iPoint)
     GradPc(2,i,j) = GradPc(2,i,j)/Vol(iPoint)
    enddo
   enddo

end subroutine compute_pcgradientgg
