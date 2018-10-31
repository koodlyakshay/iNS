subroutine pressure_residual

 use global_vars

  do i=1,Nx
    do j=1,Ny   
    
    iPoint = i + (j-1)*Nx
    if (i.ne.Nx) then
       !East
       F_e(1) = 0.5*(P(i+1,j) + P(i,j))*dy ! p*Area|_e
       F_e(2) = 0.0
       !print*,iPoint,(i+1 + (j-1)*Nx),F_e,'p'
    endif
    if (i.ne.1) then
       !West
       F_w(1) = -0.5*(P(i,j) + P(i-1,j))*dy ! p*Area|_w
       F_s(2) = 0.0
       !print*,iPoint,(i-1 + (j-1)*Nx),F_w,'p'
    endif
    if (j.ne.Ny) then
       !North
       F_n(1) = 0.0
       F_n(2) = 0.5*(P(i,j+1) + P(i,j))*dx ! p*Area|_n
       !print*,iPoint,(i + (j+1-1)*Nx),F_n,'p'
    endif
    if (j.ne.1) then
       !South
       F_s(1) = 0.0
       F_s(2) = -0.5*(P(i,j) + P(i,j-1))*dx ! p*Area|_s
       !print*,iPoint,(i + (j-1-1)*Nx),F_s,'p'
    endif
       R(1,iPoint) = R(1,iPoint) + GradU(3,1,i,j)*Vol(iPoint)!F_e(1) + F_w(1) + F_n(1) + F_s(1) ! p*Area|_e - p*Area|_w
       R(2,iPoint) = R(2,iPoint) + GradU(3,2,i,j)*Vol(iPoint)!F_e(2) + F_w(2) + F_n(2) + F_s(2) ! p*Area|_n - p*Area|_s
       
       if (wrt_data .eq. 1) then
       write(18,*)x(i,j),y(i,j),GradU(3,1,i,j)*Vol(iPoint),GradU(3,1,i,j),GradU(3,2,i,j)*Vol(iPoint),GradU(3,2,i,j)
       endif
 
    enddo
   enddo


end subroutine pressure_residual
