subroutine initialize_vars

use global_vars

do i=1,Nx
 do j=1,Ny
    x(i,j) = xmin + (i-1)*dx
    y(i,j) = ymin + (j-1)*dy
    if (x(i,j).lt.0.0) k = i
 enddo
enddo

!--- Flow defintion ---!

U_inf = 1.0
P_inf = 0.0
P_outlet = 0.0
rho = 1.0
Re_l = 1.0
mu = 1.d0/400.d0!0.798e-3
Re = (Re_l*U_inf*rho)/mu
artvisc = 4.0
CFL = 0.1
nExtIter = 200
nMIter = 500
nPIter = 100
alfa = 0.9
p_screen1 = 250
p_screen2 = 50
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

end subroutine initialize_vars
