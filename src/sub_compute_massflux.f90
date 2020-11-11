!> \file sub_compute_massflux.f90
!! \brief Subroutine to compute massflux using Rhie-Chow momentum interpolation method.

subroutine compute_massflux(extiter)

use global_vars
use output_vars
use flow_parmaters

implicit none

integer           :: i,j,iPoint,jpoint
integer           :: nPoint
integer           :: ExtIter
real              :: GradPi,  GradPav
real              :: mass_in, mass_out, mass_l2
real              :: V_e, V_w, V_n, V_s
real              :: U_e, U_w, U_n, U_s
real              :: U_up, V_up
real              :: F_e(2), F_w(2), F_n(2), F_s(2)

  nPoint = Nx*Ny
  Mass = 0.0
  mass_in = 0.0
  mass_out = 0.0
  P_Correc = 0.0
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
    enddo
   enddo

   !--- Boundary elements ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=2,Nx-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(alfa*dt_m(iPoint))
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(alfa*dt_m(iPoint))
       !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy/2.0 + 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
       
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i,j) - P(i-1,j))/dx
       F_w(1) = -rho*U_w*dy/2.0 - 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
       
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx + 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx
       
       !South
       F_s(1) = 0.0
               
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
   enddo
   
   !--- Left inlet (i=1) ---!
   i=1
   do j=2,Ny-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(alfa*dt_m(iPoint))
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(alfa*dt_m(iPoint))
      !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy + 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy
       
       !West
       F_w(1) = -rho*U(1,iPoint)*dy
       ! The mass flux is known at the inlet.
       !North 
       jPoint = i + (j+1-1)*Nx
       V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
       GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
       GradPav = (P(i,j+1) - P(i,j))/dy
       F_n(1) = rho*V_n*dx/2.0 + 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx/2.0
       
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j) - P(i,j-1))/dy
       F_s(1) = -rho*V_s*dx/2.0 - 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx/2.0
       
              
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
     mass_in = mass_in + F_w(1)
   enddo
   
   !--- Right outlet (i=Nx) ---!
   i=Nx
   do j=2,Ny-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(alfa*dt_m(iPoint))
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(alfa*dt_m(iPoint))
       !East(Outlet face)
       jPoint = i+1 + (j-1)*Nx
       U_e =  U(1,iPoint)
       GradPi = (GradU(3,1,i,j))
       GradPav = (P(i,j) - P(i-1,j))/dx
       F_e(1) = rho*U_e*dy + D(1,iPoint)*(GradPi - GradPav)*dy
       
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
       F_n(1) = rho*V_n*dx/2.0 + 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx/2.0
       
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j) - P(i,j-1))/dy
       F_s(1) = -rho*V_s*dx/2.0 - 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx/2.0
                     
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
     mass_out = mass_out + F_e(1)
   enddo
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=2,Nx-1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(alfa*dt_m(iPoint))
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(alfa*dt_m(iPoint))
      !East 
       jPoint = i+1 + (j-1)*Nx
       U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
       GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
       GradPav = (P(i+1,j) - P(i,j))/dx
       F_e(1) = rho*U_e*dy/2.0 + 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
              
       !West
       jPoint = i-1 + (j-1)*Nx
       U_w = 0.5*( U(1,jPoint) + U(1,iPoint))
       GradPi = 0.5*(GradU(3,1,i,j) + GradU(3,1,i-1,j))
       GradPav = (P(i,j) - P(i-1,j))/dx
       F_w(1) = -rho*U_w*dy/2.0 - 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
       
       !North 
       F_n(1) = 0.0
       
       !South
       jPoint = i + (j-1-1)*Nx
       V_s = 0.5*( U(2,jPoint) + U(2,iPoint))
       GradPi = 0.5*(GradU(3,2,i,j) + GradU(3,2,i,j-1))
       GradPav = (P(i,j) - P(i,j-1))/dy
       F_s(1) = -rho*V_s*dx - 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx
                     
     Mass(iPoint) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
   enddo
   
   !--- Corner elements ---!
   
   i=1
   j=1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(alfa*dt_m(iPoint))
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(alfa*dt_m(iPoint))
   !East
   jPoint = i+1 + (j-1)*Nx
   U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
   GradPav = (P(i+1,j) - P(i,j))/dx
   F_e(1) = rho*U_e*dy/2.0 + 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
      
   !North 
   jPoint = i + (j+1-1)*Nx
   V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
   GradPav = (P(i,j+1) - P(i,j))/dy
   F_n(1) = rho*V_n*dx/2.0 + 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx/2.0
       
    Mass(iPoint) = F_e(1) + F_n(1) - rho*U(1,iPoint)*dy/2.0 !Contribution from inlet
    mass_in = mass_in  - rho*U(1,iPoint)*dy/2.0
    
   i=1
   j=Ny
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(alfa*dt_m(iPoint))
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(alfa*dt_m(iPoint))
   !East
   jPoint = i+1 + (j-1)*Nx
   U_e = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j))
   GradPav = (P(i+1,j) - P(i,j))/dx
   F_e(1) = rho*U_e*dy/2.0 + 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
      
   !South 
   jPoint = i + (j-1-1)*Nx
   V_s = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j-1) + GradU(3,2,i,j))
   GradPav = (P(i,j) - P(i,j-1))/dy
   F_s(1) = -rho*V_s*dx/2.0 - 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx/2.0

    Mass(iPoint) = F_e(1) + F_s(1) - rho*U(1,iPoint)*dy/2.0 !Contribution from inlet
    mass_in = mass_in  - rho*U(1,iPoint)*dy/2.0
    
    
   i=Nx
   j=Ny
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(alfa*dt_m(iPoint))
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(alfa*dt_m(iPoint))
   !East (Outlet face)
   jPoint = i+1 + (j-1)*Nx
   U_e =  U(1,iPoint)
   GradPi = (GradU(3,1,i,j))
   GradPav = (P(i,j) - P(i-1,j))/dx
   F_e(1) = rho*U_e*dy/2.0 + D(1,iPoint)*(GradPi - GradPav)*dy/2.0
   
   !West
   jPoint = i-1 + (j-1)*Nx
   U_w = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i-1,j) + GradU(3,1,i,j))
   GradPav = (P(i,j) - P(i-1,j))/dx
   F_w(1) = -rho*U_w*dy/2.0 - 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !South 
   jPoint = i + (j-1-1)*Nx
   V_s = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j-1) + GradU(3,2,i,j))
   GradPav = (P(i,j) - P(i,j-1))/dy
   F_s(1) = -rho*V_s*dx/2.0 - 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx/2.0
    
   Mass(iPoint) = F_w(1) + F_s(1) + F_e(1)
   mass_out = mass_out  + F_e(1)
    
   i=Nx
   j=1
   iPoint = i + (j-1)*Nx
   Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(alfa*dt_m(iPoint))
   Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(alfa*dt_m(iPoint))
   !East(Outlet face)
   jPoint = i+1 + (j-1)*Nx
   U_e =  U(1,iPoint)
   GradPi = (GradU(3,1,i,j))
   GradPav = (P(i,j) - P(i-1,j))/dx
   F_e(1) = rho*U_e*dy/2.0 + D(1,iPoint)*(GradPi - GradPav)*dy/2.0

   !West
   jPoint = i-1 + (j-1)*Nx
   U_w = 0.5*( U(1,iPoint) + U(1,jPoint))
   GradPi = 0.5*(GradU(3,1,i-1,j) + GradU(3,1,i,j))
   GradPav = (P(i,j) - P(i-1,j))/dx
   F_w(1) = -rho*U_w*dy/2.0 - 0.5*(D(1,iPoint)+D(1,jPoint))*(GradPi - GradPav)*dy/2.0
   
   !North 
   jPoint = i + (j+1-1)*Nx
   V_n = 0.5*( U(2,iPoint) + U(2,jPoint))
   GradPi = 0.5*(GradU(3,2,i,j+1) + GradU(3,2,i,j))
   GradPav = (P(i,j+1) - P(i,j))/dy
   F_n(1) = rho*V_n*dx/2.0 + 0.5*(D(2,iPoint)+D(2,jPoint))*(GradPi - GradPav)*dx/2.0

    
   Mass(iPoint) = F_w(1) + F_n(1) + F_e(1)
   mass_out = mass_out  + F_e(1)


   mass_l2 = 0.0

   !open(unit=11,file='test.txt',status='unknown')
   do i=1,Nx
     do j=1,Ny
      iPoint = i+ (j-1)*Nx
      mass_l2 = mass_l2 + Mass(iPoint)**2.0
      !write(11,*) iPoint, Mass(iPoint)
     enddo
   enddo
   
   
   mass_l2 = sqrt(mass_l2/nPoint)
   if (modulo(ExtIter,p_out) .eq. 0) then 
     print*,'Mass: ',log10(mass_l2),ExtIter,'Net mass: ',mass_in + mass_out
     print*
   endif
   !close(11)
   !stop

end subroutine compute_massflux
