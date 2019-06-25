!> \file sub_flow_correction.f90
!! \brief Subroutine to correct velocity and pressure based on the pressure correction equation.

subroutine sub_flow_correction


use global_vars

implicit none

integer     :: i, j, iPoint, jPoint


   do i=2,Nx-1
    do j=2,Ny-1
      iPoint = i + (j-1)*Nx
      P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
      !East
      jPoint = i+1 + (j-1)*Nx
      F_e(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy
      !West
      jPoint = i-1 + (j-1)*Nx
      F_w(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy
      !North
      jPoint = i + (j+1-1)*Nx
      F_n(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx
      !South
      jPoint = i + (j-1-1)*Nx
      F_s(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx
      
      Vel_Corr(1,iPoint) = (F_e(1) - F_w(1))/&
                               Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
      Vel_Corr(2,iPoint) = (F_n(1) - F_s(1))/&
                               Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)

      U(1,iPoint) = U(1,iPoint) - Vel_Corr(1,iPoint)
      
      U(2,iPoint) = U(2,iPoint) - Vel_Corr(2,iPoint)
      
    enddo
   enddo
   !--- Apply BC ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
    !East
    jPoint = i+1 + (j-1)*Nx
    F_e(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy/2.0
    !West
    jPoint = i-1 + (j-1)*Nx
    F_w(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy/2.0
    !North
    jPoint = i + (j+1-1)*Nx
     F_n(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx
    !South
    jPoint = i + (j-1-1)*Nx
    F_s(1) = (P_Correc(iPoint))*dx
      
    Vel_Corr(1,iPoint) = (F_e(1) - F_w(1))/&
                               Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
    Vel_Corr(2,iPoint) = (F_n(1) - F_s(1))/&
                               Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)                            
    
    U(1,iPoint) = U(1,iPoint) - Vel_Corr(1,iPoint)
      
    U(2,iPoint) = U(2,iPoint) - Vel_Corr(2,iPoint)  
    
   enddo
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
    iPoint = i + (j-1)*Nx
    P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)

    !East
    jPoint = i+1 + (j-1)*Nx
    F_e(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy/2.0
    !West
    jPoint = i-1 + (j-1)*Nx
    F_w(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy/2.0
    !North
    jPoint = i + (j+1-1)*Nx
    F_n(1) = (P_Correc(iPoint))*dx
    !South
    jPoint = i + (j-1-1)*Nx
    F_s(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx
      
    Vel_Corr(1,iPoint) = (F_e(1) - F_w(1))/&
                               Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
    Vel_Corr(2,iPoint) = (F_n(1) - F_s(1))/&
                               Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)
     
    U(1,iPoint) = U(1,iPoint) - Vel_Corr(1,iPoint)
      
    U(2,iPoint) = U(2,iPoint) - Vel_Corr(2,iPoint)
  
   enddo
   !--- Left inlet (i=1) ---!
   i=1
   do j=1,Ny
    iPoint = i + (j-1)*Nx
    P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
   
    !East
    jPoint = i+1 + (j-1)*Nx
    F_e(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy
    if ((j.eq.Ny) .or. (j.eq.1)) F_e(1) = F_e(1)/2.0
    !West
    jPoint = i-1 + (j-1)*Nx
    F_w(1) = P_Correc(iPoint)*dy
    if ((j.eq.Ny) .or. (j.eq.1)) F_w(1) = F_w(1)/2.0
    !North
    jPoint = i + (j+1-1)*Nx
    if (j.ne.Ny) then
      F_n(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx/2.0
    else 
      F_n(1) = (P_Correc(iPoint))*dx/2.0
    endif    
    !South
    jPoint = i + (j-1-1)*Nx
    if (j .ne. 1) then 
        F_s(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx/2.0
    else
        F_s(1) = (P_Correc(iPoint))*dx/2.0
    endif
      
    Vel_Corr(1,iPoint) = 0.0
    
    Vel_Corr(2,iPoint) = 0.0
    
    U(1,iPoint) = U(1,iPoint) - Vel_Corr(1,iPoint)
      
    U(2,iPoint) = U(2,iPoint) - Vel_Corr(2,iPoint)
    
    enddo
   !--- Right outlet (i=Nx) ---!
   i=Nx
   do j=1,Ny
    F_n = 0.0
    F_s = 0.0
    F_e = 0.0
    F_w = 0.0
    iPoint = i + (j-1)*Nx
    !P(i,j) = P(i,j) + (1.0-alfa)*P_Correc(iPoint)
    P(i,j) = P_outlet
    !if ((j.ne.Ny) .or. (j.ne.1)) then
    !East
    jPoint = i+1 + (j-1)*Nx
    F_e(1) = (P_Correc(iPoint))*dy
    if ((j.eq.Ny) .or. (j.eq.1)) F_e(1) = F_e(1)/2.0
    !West
    jPoint = i-1 + (j-1)*Nx
    F_w(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dy    
    if ((j.eq.Ny) .or. (j.eq.1)) F_w(1) = F_w(1)/2.0
    !North
    jPoint = i + (j+1-1)*Nx
    if (j.ne.Ny) then
      F_n(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx/2.0
    else 
      F_n(1) = (P_Correc(iPoint))*dx/2.0
    endif
    !South
    jPoint = i + (j-1-1)*Nx
    if (j.ne.1) then
      F_s(1) = 0.5*(P_Correc(iPoint) + P_Correc(jPoint))*dx/2.0
    else 
      F_s(1) = (P_Correc(iPoint))*dx/2.0
    endif
      
    Vel_Corr(1,iPoint) = (F_e(1) - F_w(1))/&
                               Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
    Vel_Corr(2,iPoint) = (F_n(1) - F_s(1))/&
                               Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)
    
    U(1,iPoint) = U(1,iPoint) - Vel_Corr(1,iPoint)
      
    U(2,iPoint) = U(2,iPoint) - Vel_Corr(2,iPoint)
   enddo

end subroutine sub_flow_correction
