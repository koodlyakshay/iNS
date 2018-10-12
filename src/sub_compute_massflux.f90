subroutine compute_massflux

use global_vars


Mass = 0.d0
mass_in = 0.0
mass_out = 0.0
P_Correc = 0.d0
   do i=2,Nx-1
    do j=2,Ny-1
    !--- Point definition ---!
    iPoint = i + (j-1)*Nx
    
    !--- Mass flux ---!
    !--- Rhie-Chow interpolation ---!
       !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy + 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i,j) - P(i-1,j))/dx
       F_w(1) = -rho*U_w*dy - 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx + 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j) - P(i,j-1))/dy
       F_s(1) = -rho*V_s*dx - 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx
              
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
       !print*,iPoint,F_e,F_w,F_n,F_s,Mass(iPoint)
    enddo
   enddo
   

   !--- Boundary elements ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=2,Nx-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*Vol/(alfa*dt_m)
       !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy/2.d0 + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.d0
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i,j) - P(i-1,j))/dx
       F_w(1) = -rho*U_w*dy/2.d0 - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.d0
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx
       !South
       F_s(1) = 0.d0
               
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
   enddo
   
   !--- Left inlet (i=1) ---!
   i=1
   do j=2,Ny-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*Vol/(alfa*dt_m)
      !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       !West
       F_w(1) = -rho*U(1,iPoint)*dy
       ! The mass flux is known at the inlet.
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx/2.d0 + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx/2.0
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j) - P(i,j-1))/dy
       F_s(1) = -rho*V_s*dx/2.d0 - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx/2.0
              
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
     mass_in = mass_in + F_w(1)
   enddo
   
   !--- Right outlet (i=Nx) ---!
   ! pressure is specified at the outlet. 
   ! No need to solve the equations.
   i=Nx
   do j=2,Ny-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*Vol/(alfa*dt_m)
      !East 
       F_e(1) = rho*U(1,iPoint)*dy
       
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i,j) - P(i-1,j))/dx
       F_w(1) = -rho*U_w*dy - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx/2.d0 + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx/2.0
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j) - P(i,j-1))/dy
       F_s(1) = -rho*V_s*dx/2.d0 - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx/2.0
              
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
     mass_out = mass_out + F_e(1)
   enddo
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.5*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.5*Vol/(alfa*dt_m)
      !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy/2.d0 + 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
       
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i,j) - P(i-1,j))/dx
       F_w(1) = -rho*U_w*dy/2.d0 - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
       
       !North 
       F_n(1) = 0.d0
       
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j) - P(i,j-1))/dy
       F_s(1) = -rho*V_s*dx - 0.5*(0.5*D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dx
              
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
   enddo
   
   !--- Corner elements ---!
   
   i=1
   j=1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.25*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.25*Vol/(alfa*dt_m)
   !East
   jPoint = i+1 + (j-1)*Nx
   U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
   GradPav = (P(i+1,j) - P(i,j))/dx
   F_e(1) = rho*U_e*dy/2.d0 + 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !North 
   jPoint = i + (j+1-1)*Nx
   V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
   GradPav = (P(i,j+1) - P(i,j))/dy
   F_n(1) = rho*V_n*dx/2.d0 + 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dx/2.0
    
    Mass(iPoint) = F_e(1) + F_n(1) - rho*U(1,iPoint)*dy/2.0 !Contribution from inlet
    mass_in = mass_in  - rho*U(1,iPoint)*dy/2.0
    
   i=1
   j=Ny
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.25*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.25*Vol/(alfa*dt_m)
   !East
   jPoint = i+1 + (j-1)*Nx
   U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
   GradPav = (P(i+1,j) - P(i,j))/dx
   F_e(1) = rho*U_e*dy/2.d0 + 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !South 
   jPoint = i + (j-1-1)*Nx
   V_s = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j-1) + GradU(3,2,i,j))
   GradPav = (P(i,j) - P(i,j-1))/dy
   F_s(1) = -rho*V_s*dx/2.d0 - 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dx/2.0
    
    Mass(iPoint) = F_e(1) + F_s(1) - rho*U(1,iPoint)*dy/2.0 !Contribution from inlet
    mass_in = mass_in  - rho*U(1,iPoint)*dy/2.0
    
   i=Nx
   j=Ny
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.25*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.25*Vol/(alfa*dt_m)
   !West
   jPoint = i-1 + (j-1)*Nx
   U_w = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i-1,j) + GradU(3,1,i,j))
   GradPav = (P(i,j) - P(i-1,j))/dx
   F_w(1) = -rho*U_w*dy/2.d0 - 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !South 
   jPoint = i + (j-1-1)*Nx
   V_s = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j-1) + GradU(3,2,i,j))
   GradPav = (P(i,j) - P(i,j-1))/dy
   F_s(1) = -rho*V_s*dx/2.d0 - 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dx/2.0
    
   Mass(iPoint) = F_w(1) + F_s(1) + rho*U(1,iPoint)*dy/2.0
   mass_out = mass_out  + rho*U(1,iPoint)*dy/2.0
    
    
   i=Nx
   j=1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 0.25*Vol/(alfa*dt_m)
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 0.25*Vol/(alfa*dt_m)
   !West
   jPoint = i-1 + (j-1)*Nx
   U_w = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i-1,j) + GradU(3,1,i,j))
   GradPav = (P(i,j) - P(i-1,j))/dx
   F_w(1) = -rho*U_w*dy/2.d0 - 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !North 
   jPoint = i + (j+1-1)*Nx
   V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
   GradPav = (P(i,j+1) - P(i,j))/dy
   F_n(1) = rho*V_n*dx/2.d0 + 0.5*(0.25*D(1,iPoint)+0.5*D(1,jPoint))*(GradPi - GradPav)*dx/2.0
    
   Mass(iPoint) = F_w(1) + F_n(1) + rho*U(1,iPoint)*dy/2.0
   mass_out = mass_out  + rho*U(1,iPoint)*dy/2.0
   
   mass_l2 = 0.d0

   do iPoint=1,nPoint
      mass_l2 = mass_l2 + Mass(iPoint)**2.0
   enddo
   mass_l2 = sqrt(mass_l2/nPoint)
   if (modulo(ExtIter,p_out) .eq. 0) then 
     print*,'Mass: ',mass_l2,ExtIter
     print*,'Net mass: ',mass_in + mass_out
   endif

end subroutine compute_massflux
