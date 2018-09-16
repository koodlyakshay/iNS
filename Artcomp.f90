!Program to solve Lid driven cavity using SIMPLE 
!Implementation based on SU2 
!Governing Eqn:
!             \partial_t U + \nabla F_c - \nabla F_v = Q
!U -> Solution vector, F_c -> Convective flux, F_v -> Viscous flux, Q -> Source term
!Discretized Eqn:
!(\frac{\Omega_i}{\Delat t^n} + \frac{\partial R_i(U^n)}{\partial U_j}).\Delata U_j^n = -R_i(U^n)
!Residual, R_i(U) = \Sigma_{j\in N(i)} (F_{c_{ij}} - F_{v_{ij}})\Delta S_{ij} - Q_i|\Omega_i|
!Jacobian, J_{ij} = \frac{\partial R_i(U^n)}{\partial U_j}
program iNS

implicit none

integer             :: i,j,k,Nx,Ny,iPoint,jpoint,nPoint,nDim,jDim,nVar,iVar
integer             :: ExtIter,IniIter,nExtIter,nIniIter,convergence,liniter
real*8              :: dx,dy,vol,area,Lx,Ly,xmin,ymin,Res_l2,alfa,beta2,kappa
real*8              :: dt,CFL,lambda_inv,lambda_visc,dt_i,dt_v,lambda_i,lambda_j
real*8              :: Param_p,sc0,SF,Phi_i,Phi_j,lambda_mean,E_0,U_max
real*8              :: Re, U_inf,rho,mu,Re_l,P_inf,GradPi,GradPav
real*8              :: U_e,U_w,U_n,U_s,F_e(3),F_w(3),F_n(3),F_s(3),Gr_e,Gr_w,Gr_n,Gr_s
real*8              :: V_e,V_w,V_n,V_s,P_e,P_w,P_n,P_s
real*8,allocatable  :: x(:,:),y(:,:),Mass(:)
real*8,allocatable  :: U(:,:),U_old(:,:),Sol(:,:),Fc(:,:,:),Fv(:,:,:)
real*8,allocatable  :: R(:,:),Jac(:,:,:),P(:,:),V(:,:,:),Res(:),Mat(:,:)
real*8,allocatable  :: GradU(:,:,:,:),Solp(:),Jacp(:,:),Rp(:),Tot_Jac(:,:),Tot_R(:),Tot_Sol(:)
logical             :: implicit_time
!--- Output Files ---!
open(unit=10,file='Solver_details_AC.txt',status='unknown')
open(unit=11,file='Jacobian.txt',status='unknown')

open(unit=15,file='Gradients.txt',status='unknown')
open(unit=16,file='Viscbudgets.txt',status='unknown')
open(unit=17,file='Convbudgets.txt',status='unknown')
!--- Grid definition ---!
! *---*---*---*---*   j=Ny
! *---*---*---*---*    .
! *---*---*---*---*    .
! *---*---*---*---*    .
! *---*---*---*---*   j=1
!i=1,...........,i=Nx

Lx = 1.0
Ly = 1.0
Nx = 17
Ny = 17
nPoint = Nx*Ny
nDim = 2
xmin = 0.0
ymin = 0.0
dx = Lx/(Nx-1)
dy = Ly/(Ny-1)
area = dx
Vol = dx*dy
nVar = nDim + 1
allocate(x(Nx,Ny),y(Nx,Ny))
allocate(U(3,Nx*Ny),U_old(3,Nx*Ny),V(3,Nx,Ny),Fc(3,Nx,Ny),Fv(3,Nx,Ny))
allocate(R(3,Nx*Ny),Jac(3,Nx*Ny,Nx*Ny),GradU(3,2,Nx,Ny),Res(Nx*Ny*nVar),Mass(Nx*Ny))
allocate(Tot_Jac(Nx*Ny*nVar,Nx*Ny*nVar),Tot_R(Nx*Ny*nVar),Tot_Sol(Nx*Ny*nVar))

do i=1,Nx
 do j=1,Ny
    x(i,j) = xmin + (i-1)*dx
    y(i,j) = ymin + (j-1)*dy
 enddo
enddo

!--- Flow defintion ---!

U_inf = 1.0
P_inf = 0.0
rho = 1.0
Re_l = Lx
mu = 1.d0/400.d0
Re = (Re_l*U_inf*rho)/mu
beta2 = 1.0
CFL = 1.0
nExtIter = 50000
alfa = 1.0
!--- Initialize variables ---!

do i=1,Nx
 do j=1,Ny
    V(1:3,i,j) = 0.0
    Fc(1:3,i,j) = 0.0
    Fv(1:3,i,j) = 0.0
 enddo
enddo

!--- 1 -> P, 2 -> u, 3 -> v ---!
do iPoint = 1,nPoint 
 U(1,iPoint) = P_inf
 U_old(1,iPoint) = P_inf
 U(2,iPoint) = 0.0
 U_old(2,iPoint) = 0.0
 U(3,iPoint) = 0.0
 U_old(3,iPoint) = 0.0
enddo
i=(Nx+1)/2
j=1
iPoint = i + (j-1)*Nx
U(1,iPoint) = 0.0
!--- Set Time step ---!
lambda_inv = U_inf*dy + sqrt(U_inf*U_inf*dy*dy + (beta2/rho)*area*area) 
lambda_visc = mu*area*area/rho
dt_i = CFL*Vol/lambda_inv
dt_v = CFL*0.25*Vol*Vol/lambda_visc
!dt = min(dt_i,dt_v)
dt = 1.0e-4
implicit_time = .false.
Param_p = 0.1
kappa = 0.03
write(10,*)'Lambda_inv: ',lambda_inv
write(10,*)'Lambda_visc: ',lambda_visc
write(10,*)'dt_i: ',dt_i
write(10,*)'dt_v: ',dt_v
write(10,*)'dt: ',dt
write(10,*)'Area',area
write(10,*)'Vol: ',Vol
write(10,*)'Re: ',Re
write(10,*)'Artificial Comprssibility: ',beta2
close(10)
!--- Begin Solver ---!

do ExtIter = 1,nExtIter
   
   R = 0.0
   Tot_Jac = 0.0
   !--- Set old solution ---!
  
   do iPoint=1,nPoint
     U_old(1:3,iPoint) = U(1:3,iPoint)
   enddo
   
   
   !--- Set Primitive ---!
   !--- 1 -> P, 2 -> u, 3 ->v ---!
   do i=1,Nx
    do j=1,Ny
    iPoint = i + (j-1)*Nx
      V(1:3,i,j) = U_old(1:3,iPoint)
    enddo
   enddo
   
   !--- Compute gradients at nodes, to be used for viscous flux ---!
   do i=2,Nx-1
    do j=1,Ny
       !--- (var,dim,i,j) ---!
       GradU(1,1,i,j) = (V(1,i+1,j) - V(1,i,j))/(1.0*dx)
       GradU(2,1,i,j) = (V(2,i+1,j) - V(2,i,j))/(1.0*dx)
       GradU(3,1,i,j) = (V(3,i+1,j) - V(3,i,j))/(1.0*dx)
       !print*,i + (j-1)*Nx,'loop 1 d/dx'
    enddo
   enddo
   
   do i=1,Nx
    do j=2,Ny-1
       !--- (var,dim,i,j) ---!
       GradU(1,2,i,j) = (V(1,i,j+1) - V(1,i,j))/(1.0*dy)
       GradU(2,2,i,j) = (V(2,i,j+1) - V(2,i,j))/(1.0*dy)
       GradU(3,2,i,j) = (V(3,i,j+1) - V(3,i,j))/(1.0*dy)
       !print*,i + (j-1)*Nx,'loop 2 d/dy'
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
   !print*,i + (j-1)*Nx,'corner nx,ny'
   
   i=1
   j=1
   !--- (var,dim,i,j) ---!
   GradU(1,1,i,j) = (V(1,i+1,j) - V(1,i,j))/dx
   GradU(1,2,i,j) = (V(1,i,j+1) - V(1,i,j))/dy
   GradU(2,1,i,j) = (V(2,i+1,j) - V(2,i,j))/dx
   GradU(2,2,i,j) = (V(2,i,j+1) - V(2,i,j))/dy
   GradU(3,1,i,j) = (V(3,i+1,j) - V(3,i,j))/dx
   GradU(3,2,i,j) = (V(3,i,j+1) - V(3,i,j))/dy
   !print*,i + (j-1)*Nx,'corner nx,ny'

   i=1
   j=Ny
   !--- (var,dim,i,j) ---!
   GradU(1,1,i,j) = (V(1,i+1,j) - V(1,i,j))/dx
   GradU(1,2,i,j) = (V(1,i,j) - V(1,i,j-1))/dy
   GradU(2,1,i,j) = (V(2,i+1,j) - V(2,i,j))/dx
   GradU(2,2,i,j) = (V(2,i,j) - V(2,i,j-1))/dy
   GradU(3,1,i,j) = (V(3,i+1,j) - V(3,i,j))/dx
   GradU(3,2,i,j) = (V(3,i,j) - V(3,i,j-1))/dy
   !print*,i + (j-1)*Nx,'corner nx,ny'
      
   i=Nx
   j=1
   !--- (var,dim,i,j) ---!
   GradU(1,1,i,j) = (V(1,i,j) - V(1,i-1,j))/dx
   GradU(1,2,i,j) = (V(1,i,j+1) - V(1,i,j))/dy
   GradU(2,1,i,j) = (V(2,i,j) - V(2,i-1,j))/dx
   GradU(2,2,i,j) = (V(2,i,j+1) - V(2,i,j))/dy
   GradU(3,1,i,j) = (V(3,i,j) - V(3,i-1,j))/dx
   GradU(3,2,i,j) = (V(3,i,j+1) - V(3,i,j))/dy
   !print*,i + (j-1)*Nx,'corner nx,ny'
 
   
   !--- Compute spatial discretization ---!
   !write(16,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   !write(17,*)'-------------------------------------------Iteration-------------------------------------------------',ExtIter
   !--- Interior Nodes ---!
   do i=2,Nx-1
    do j=2,Ny-1
    !--- Node definition ---!
    iPoint = i + (j-1)*Nx
    
    !--- Convective terms (central) ---!    
       !--- Compute flux and update Jacobians ---!
       !East 
       jPoint = i+1 + (j-1)*Nx
       P_e = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       U_e = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       V_e = 0.5*( U_old(3,iPoint) + U_old(3,jPoint))
       
       F_e(1) = beta2*rho*U_e*dy              ! (rho*U)|_e
       F_e(2) = rho*U_e*U_e*dy + P_e*dy ! (rho*U)*U|_e + P
       F_e(3) = rho*V_e*U_e*dy          ! (rho*V)*U|_e
       
       lambda_i = abs(U_old(2,iPoint)*dy) + sqrt(U_old(2,iPoint)*U_old(2,iPoint)*dy*dy + (beta2/rho)*area*area) 
       lambda_j = abs(U_old(2,jPoint)*dy) + sqrt(U_old(2,jPoint)*U_old(2,jPoint)*dy*dy + (beta2/rho)*area*area) 
       lambda_mean = 0.5*(lambda_i + lambda_j)
       !print*,'lambda_i ',lambda_i,', lambda_j ',lambda_j,' E'
       Phi_i = (lambda_i/(4.d0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.d0*lambda_mean))**Param_p
       SF = 4.d0*Phi_i*Phi_j/(Phi_i + Phi_j)
       sc0 = 3.d0*(4.d0+4.d0)/(4.d0*4.d0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.d0/3.d0
       F_e(1) = F_e(1) + E_0*(U_old(1,iPoint) - U_old(1,jPoint))*SF*lambda_mean
       F_e(2) = F_e(2) + E_0*(U_old(2,iPoint) - U_old(2,jPoint))*SF*lambda_mean
       F_e(3) = F_e(3) + E_0*(U_old(3,iPoint) - U_old(3,jPoint))*SF*lambda_mean
       
       !--- iPoint, jPoint ---!
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) + 0.0 - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) + 0.5*dy*beta2
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) + 0.5*dy 
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) + 1.0*U_e*dy - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) + 0.5*V_e*dy
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) + 0.5*U_e*dy - E_0*SF*lambda_mean
       
       !--- iPoint, iPoint ---!
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.0 + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.5*dy*beta2
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.5*dy 
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 1.0*U_e*dy + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) + 0.5*V_e*dy
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) + 0.5*U_e*dy + E_0*SF*lambda_mean
       
       !West
       jPoint = i-1 + (j-1)*Nx
       P_w = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       U_w = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       V_w = 0.5*( U_old(3,iPoint) + U_old(3,jPoint))
       F_w(1) = -beta2*rho*U_w*dy              ! (rho*U)|_w
       F_w(2) = -rho*U_w*U_w*dy - P_w*dy ! (rho*U)*U|_w + P
       F_w(3) = -rho*V_w*U_w*dy          ! (rho*V)*U|_w
       
       lambda_i = abs(U_old(2,iPoint)*dy) + sqrt(U_old(2,iPoint)*U_old(2,iPoint)*dy*dy + (beta2/rho)*area*area) 
       lambda_j = abs(U_old(2,jPoint)*dy) + sqrt(U_old(2,jPoint)*U_old(2,jPoint)*dy*dy + (beta2/rho)*area*area) 
       lambda_mean = 0.5*(lambda_i + lambda_j)
       !print*,'lambda_i ',lambda_i,', lambda_j ',lambda_j,' W'
       Phi_i = (lambda_i/(4.d0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.d0*lambda_mean))**Param_p
       SF = 4.d0*Phi_i*Phi_j/(Phi_i + Phi_j)
       sc0 = 3.d0*(4.d0+4.d0)/(4.d0*4.d0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.d0/3.d0
       F_w(1) = F_w(1) + E_0*(U_old(1,iPoint) - U_old(1,jPoint))*SF*lambda_mean
       F_w(2) = F_w(2) + E_0*(U_old(2,iPoint) - U_old(2,jPoint))*SF*lambda_mean
       F_w(3) = F_w(3) + E_0*(U_old(3,iPoint) - U_old(3,jPoint))*SF*lambda_mean
       
       !--- iPoint, jPoint ---!
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.0 - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.5*dy*beta2
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.5*dy 
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - 1.0*U_w*dy - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) - 0.5*V_w*dy
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) - 0.5*U_w*dy - E_0*SF*lambda_mean
       
       !--- iPoint, iPoint ---!
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) - 0.0 + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) - 0.5*dy*beta2
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) - 0.5*dy
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) - 1.0*U_w*dy + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) - 0.5*V_w*dy
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) - 0.5*U_w*dy + E_0*SF*lambda_mean
       
       !North 
       jPoint = i + (j+1-1)*Nx
       P_n = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       U_n = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       V_n = 0.5*( U_old(3,iPoint) + U_old(3,jPoint))
       F_n(1) = beta2*rho*V_n*dx              ! (rho*V)|_n
       F_n(2) = rho*U_n*V_n*dx          ! (rho*U)*V|_n
       F_n(3) = rho*V_n*V_n*dx + P_n*dx ! (rho*V)*V|_n + P
       
       lambda_i = abs(U_old(3,iPoint)*dx) + sqrt(U_old(3,iPoint)*U_old(3,iPoint)*dx*dx + (beta2/rho)*area*area) 
       lambda_j = abs(U_old(3,jPoint)*dx) + sqrt(U_old(3,jPoint)*U_old(3,jPoint)*dx*dx + (beta2/rho)*area*area) 
       lambda_mean = 0.5*(lambda_i + lambda_j)
       !print*,'lambda_i ',lambda_i,', lambda_j ',lambda_j,' N'
       Phi_i = (lambda_i/(4.d0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.d0*lambda_mean))**Param_p
       SF = 4.d0*Phi_i*Phi_j/(Phi_i + Phi_j)
       sc0 = 3.d0*(4.d0+4.d0)/(4.d0*4.d0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.d0/3.d0
       F_n(1) = F_n(1) + E_0*(U_old(1,iPoint) - U_old(1,jPoint))*SF*lambda_mean
       F_n(2) = F_n(2) + E_0*(U_old(2,iPoint) - U_old(2,jPoint))*SF*lambda_mean
       F_n(3) = F_n(3) + E_0*(U_old(3,iPoint) - U_old(3,jPoint))*SF*lambda_mean
       
       !--- iPoint, jPoint ---!
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) + 0.0 - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) + 0.5*dx*beta2
       
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) + 0.5*V_n*dx - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) + 0.5*U_n*dy
       
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) + 0.5*dx 
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) + 1.0*V_n*dx - E_0*SF*lambda_mean
       
       !--- iPoint, iPoint ---!
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.0 + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) + 0.5*dx*beta2
       
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + 0.5*V_n*dx + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) + 0.5*U_n*dy
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) + 0.5*dx
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) + 1.0*V_n*dy + E_0*SF*lambda_mean
       
       !South
       jPoint = i + (j-1-1)*Nx
       P_s = 0.5*( U_old(1,iPoint) + U_old(1,jPoint))
       U_s = 0.5*( U_old(2,iPoint) + U_old(2,jPoint))
       V_s = 0.5*( U_old(3,iPoint) + U_old(3,jPoint))
       F_s(1) = -beta2*rho*V_s*dx              ! (rho*V)|_s
       F_s(2) = -rho*U_s*V_s*dx          ! (rho*U)*V|_s
       F_s(3) = -rho*V_s*V_s*dx - P_s*dx ! (rho*V)*V|_s + P
       
       lambda_i = abs(U_old(3,iPoint)*dx) + sqrt(U_old(3,iPoint)*U_old(3,iPoint)*dx*dx + (beta2/rho)*area*area) 
       lambda_j = abs(U_old(3,jPoint)*dx) + sqrt(U_old(3,jPoint)*U_old(3,jPoint)*dx*dx + (beta2/rho)*area*area) 
       lambda_mean = 0.5*(lambda_i + lambda_j)
       !print*,'lambda_i ',lambda_i,', lambda_j ',lambda_j,' S'
       Phi_i = (lambda_i/(4.d0*lambda_mean))**Param_p
       Phi_j = (lambda_j/(4.d0*lambda_mean))**Param_p
       SF = 4.d0*Phi_i*Phi_j/(Phi_i + Phi_j)
       sc0 = 3.d0*(4.d0+4.d0)/(4.d0*4.d0)  ! (Neighbor_i + Neighbor_j)/(Neighbor_i*Neighbor_j)
       E_0 = kappa*sc0*2.d0/3.d0
       F_s(1) = F_s(1) + E_0*(U_old(1,iPoint) - U_old(1,jPoint))*SF*lambda_mean
       F_s(2) = F_s(2) + E_0*(U_old(2,iPoint) - U_old(2,jPoint))*SF*lambda_mean
       F_s(3) = F_s(3) + E_0*(U_old(3,iPoint) - U_old(3,jPoint))*SF*lambda_mean
       
       !--- iPoint, jPoint ---!
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.0 - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) - 0.5*dx*beta2
       
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - 0.5*V_s*dx - E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) - 0.5*U_s*dx
       
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) - 0.5*dx
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) - 1.0*V_s*dy - E_0*SF*lambda_mean
       
       !--- iPoint, iPoint ---!
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) - 0.0 + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) - 0.5*dx*beta2
       
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) - 0.5*V_s*dx + E_0*SF*lambda_mean
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) - 0.5*U_s*dx
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) - 0.5*dx
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) - 1.0*V_s*dy + E_0*SF*lambda_mean
       
       
       !--- Note ---!
       ! diagonal elements of jac are 0.5*2.0*V_s but simplified as 1.0 
       !--- Update residuals ---!
       
       Fc(1,i,j) = F_e(1) + F_w(1) + F_n(1) + F_s(1) ! (rho*U)|_e - (rho*U)|_w + (rho*V)|_n - (rho*V)|_s
       Fc(2,i,j) = F_e(2) + F_w(2) + F_n(2) + F_s(2) ! (rho*U)*U|_e - (rho*U)*U|_w + (rho*U)*V|_n - (rho*U)*V|_s
       Fc(3,i,j) = F_e(3) + F_w(3) + F_n(3) + F_s(3) ! (rho*V)*U|_e - (rho*V)*U|_w + (rho*V)*V|_n - (rho*V)*V|_s
       
       !write(17,*)F_e(1:3),'iPoint',iPoint,'jPoint',i+1 + (j-1)*Nx,i+1,j,U_e
       !write(17,*)F_w(1:3),'iPoint',iPoint,'jPoint',i-1 + (j-1)*Nx,i+1,j,U_w
       !write(17,*)F_n(1:3),'iPoint',iPoint,'jPoint',i + (j+1-1)*Nx,i+1,j,U_n
       !write(17,*)F_s(1:3),'iPoint',iPoint,'jPoint',i + (j-1-1)*Nx,i+1,j,U_s
       !write(17,*)
    !--- Viscous terms ---!
       !East
       F_e(1) = 0.0
       Gr_e = 1.0*0.5*(GradU(2,1,i+1,j) + GradU(2,1,i,j))
       F_e(2) = mu*Gr_e*dy ! du/dx|_e
       Gr_e = 1.0*0.5*(GradU(3,1,i+1,j) + GradU(3,1,i,j)) !+ 1.0*0.5*(GradU(2,2,i+1,j) + GradU(2,2,i,j))
       F_e(3) = mu*Gr_e*dy ! dv/dx|_e + du/dy|_e
       !West
       F_w(1) = 0.0
       Gr_w = 1.0*0.5*(GradU(2,1,i-1,j) + GradU(2,1,i,j)) 
       F_w(2) = -mu*Gr_w*dy ! du/dx|_w
       Gr_w = 1.0*0.5*(GradU(3,1,i-1,j) + GradU(3,1,i,j)) !+ 1.0*0.5*(GradU(2,2,i-1,j) + GradU(2,2,i,j)) 
       F_w(3) = -mu*Gr_w*dy ! dv/dx|_w + du/dy|_w
       !North
       F_n(1) = 0.0
       Gr_n = 1.0*0.5*( GradU(2,2,i,j+1) + GradU(2,2,i,j)) !+ 1.0*0.5*(GradU(3,1,i,j+1) + GradU(3,1,i,j)) 
       F_n(2) = mu*Gr_n*dx ! dv/dx|_n + du/dy|_n
       Gr_n = 1.0*0.5*( GradU(3,2,i,j+1) + GradU(3,2,i,j)) 
       F_n(3) = mu*Gr_n*dx ! dv/dy|_n
       !South
       F_s(1) = 0.0
       Gr_s = 1.0*0.5*( GradU(2,2,i,j-1) + GradU(2,2,i,j)) !+ 1.0*0.5*(GradU(3,1,i,j-1) + GradU(3,1,i,j))
       F_s(2) = -mu*Gr_s*dx ! dv/dx|_s + du/dy|_s
       Gr_s = 1.0*0.5*( GradU(3,2,i,j-1) + GradU(3,2,i,j)) 
       F_s(3) = -mu*Gr_s*dx ! dv/dy|_s
              
       Fv(1,i,j) = F_e(1) + F_w(1) + F_n(1) + F_s(1)
       Fv(2,i,j) = F_e(2) + F_w(2) + F_n(2) + F_s(2) ! du/dx|_e - du/dx|_w + du/dy|_n - du/dy|_s
       Fv(3,i,j) = F_e(3) + F_w(3) + F_n(3) + F_s(3) ! dv/dx|_e - dv/dx|_w + dv/dy|_n - dv/dy|_s
       
       !write(16,*)F_e(1:3),'iPoint',iPoint,i,j,'jPoint',i+1 + (j-1)*Nx,i+1,j
       !write(16,*)F_w(1:3),'iPoint',iPoint,i,j,'jPoint',i-1 + (j-1)*Nx,i-1,j
       !write(16,*)F_n(1:3),'iPoint',iPoint,i,j,'jPoint',i + (j+1-1)*Nx,i,j+1
       !write(16,*)F_s(1:3),'iPoint',iPoint,i,j,'jPoint',i + (j-1-1)*Nx,i,j-1
       !write(16,*)
       !--- Update Jacobians ---!
       !East
       jPoint = i+1 + (j-1)*Nx
       
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) - mu*area/dx
       
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) + mu*area/dx
       
       !West
       jPoint = i-1 + (j-1)*Nx
       
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) - mu*area/dx
       
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) + mu*area/dx
       
       !North
       jPoint = i + (j+1-1)*Nx
       
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(jPoint-1)*nVar+3) - mu*area/dx
       
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) + mu*area/dx
       
       !South
       jPoint = i + (j-1-1)*Nx
       
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+2) - mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(jPoint-1)*nVar+3) - 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) - 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) - mu*area/dx
       
       
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + mu*area/dx
       Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+3) + 0.0
       
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+1) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+2) + 0.0
       Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) + mu*area/dx
       
       !--- Update residual ---!
       
       R(1,iPoint) = Fc(1,i,j) - Fv(1,i,j)
       R(2,iPoint) = Fc(2,i,j) - Fv(2,i,j)
       R(3,iPoint) = Fc(3,i,j) - Fv(3,i,j)
    enddo
   enddo
   
   !--- Apply BC ---!
   !--- Lower wall (j=1) ---!
   j=1
   do i=1,Nx
   iPoint = i + (j-1)*Nx
    !--- Zero velocity ---!
     U_old(2:3,iPoint) = 0.0
     jPoint = i + (j+1-1)*Nx 
     !U_old(1,iPoint) = U_old(1,jPoint) 
     R(1,iPoint) = beta2*rho*0.5*(U_old(3,jPoint))*dx
     jPoint = i+1 + (j-1)*Nx 
     R(1,iPoint) = R(1,iPoint) + beta2*rho*0.5*(U_old(2,jPoint))*dy
     jPoint = i-1 + (j-1)*Nx 
     R(1,iPoint) = R(1,iPoint) - beta2*rho*0.5*(U_old(2,jPoint))*dy
     
     !--- Update residual ---!
     R(1:3,iPoint) = 0.0
     
     Tot_Jac((iPoint-1)*nVar+1,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+2,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+3,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 1.0
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 1.0
     Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = 1.0

   enddo
   
   !--- Left wall (i=1) ---!
   i=1
   do j=1,Ny
   iPoint = i + (j-1)*Nx
    !--- Zero velocity ---!
     U_old(2:3,iPoint) = 0.0
     jPoint = i+1 + (j-1)*Nx 
     !U_old(1,iPoint) = U_old(1,jPoint)
     R(1,iPoint) = beta2*rho*0.5*(U_old(2,jPoint))*dy
     jPoint = i + (j+1-1)*Nx 
     R(1,iPoint) = R(1,iPoint) + beta2*rho*0.5*(U_old(3,jPoint))*dx
     jPoint = i + (j-1-1)*Nx 
     R(1,iPoint) = R(1,iPoint) - beta2*rho*0.5*(U_old(3,jPoint))*dx
     
     !--- Update residual ---!
     R(1:3,iPoint) = 0.0
     
     Tot_Jac((iPoint-1)*nVar+1,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+2,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+3,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 1.0
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 1.0
     Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = 1.0
   enddo
   
   !--- Right wall (i=Nx) ---!
   i=Nx
   do j=1,Ny
   iPoint = i + (j-1)*Nx
    !--- Zero velocity ---!
     U_old(2:3,iPoint) = 0.0
     jPoint = i-1 + (j-1)*Nx 
     !U_old(1,iPoint) = U_old(1,jPoint)
     R(1,iPoint) = -beta2*rho*0.5*(U_old(2,jPoint))*dy
     jPoint = i + (j+1-1)*Nx 
     R(1,iPoint) = R(1,iPoint) + beta2*rho*0.5*(U_old(3,jPoint))*dx
     jPoint = i + (j-1-1)*Nx 
     R(1,iPoint) = R(1,iPoint) - beta2*rho*0.5*(U_old(3,jPoint))*dx
     
     !--- Update residual ---!
     R(1:3,iPoint) = 0.0
     
     Tot_Jac((iPoint-1)*nVar+1,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+2,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+3,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 1.0
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 1.0
     Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = 1.0
   enddo
   
   !--- Upper wall (j=Ny) ---!
   j=Ny
   do i=1,Nx
   iPoint = i + (j-1)*Nx
    !--- Moving wall ---!
     U_old(2,iPoint) = U_inf
     U_old(3,iPoint) = 0.0
     jPoint = i + (j-1 -1)*Nx ! Normal gradient of pressure is zero
     !U_old(1,iPoint) = U_old(1,jPoint)
     R(1,iPoint) = -beta2*rho*0.5*(U_old(3,jPoint)+U_old(3,jPoint))*dx
     jPoint = i+1 + (j-1)*Nx 
     R(1,iPoint) = R(1,iPoint) + beta2*rho*0.5*(U_old(2,jPoint)+U_old(2,jPoint))*dy
     jPoint = i-1 + (j-1)*Nx 
     R(1,iPoint) = R(1,iPoint) - beta2*rho*0.5*(U_old(2,jPoint)+U_old(2,jPoint))*dy
     
     !--- Update residual ---!
     R(1:3,iPoint) = 0.0
     
     Tot_Jac((iPoint-1)*nVar+1,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+2,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+3,:) = 0.0
     Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = 1.0
     Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = 1.0
     Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = 1.0
   enddo
   
   Res_l2 = 0.d0
   if (implicit_time) then
   !--- Time Integration (Implicit) ---!
   !write(11,*)'---------------------Iteration----------------------------------',ExtIter
     do iPoint=1,nPoint
        Tot_R((iPoint-1)*nVar+1) = -R(1,iPoint)
        Tot_R((iPoint-1)*nVar+2) = -R(2,iPoint)
        Tot_R((iPoint-1)*nVar+3) = -R(3,iPoint)
        Res_l2 = Res_l2 + (R(2,iPoint)**2.0)
        U_max = abs(U_old(2,iPoint)*dy) + abs(U_old(3,iPoint)*dx)!max(abs(U_old(2,iPoint)),abs(U_old(3,iPoint)))
        lambda_inv = U_max + sqrt(U_max*dy*U_max*dy + (beta2/rho)*area*area) 
        if (abs(lambda_inv) .lt. 1.0e-6) lambda_inv = 2.d0*U_inf
        !dt = CFL*Vol/lambda_inv
                
        Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + Vol/dt
        Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + Vol/dt
        Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) = Tot_Jac((iPoint-1)*nVar+3,(iPoint-1)*nVar+3) + Vol/dt 
        
        Tot_Sol((iPoint-1)*nVar+1) = 0.0
        Tot_Sol((iPoint-1)*nVar+2) = 0.0
        Tot_Sol((iPoint-1)*nVar+3) = 0.0
        
        !write(11,*)Tot_Jac((iPoint-1)*nVar+1,:)
        !write(11,*)Tot_Jac((iPoint-1)*nVar+1,:)
        !write(11,*)Tot_Jac((iPoint-1)*nVar+1,:)
        !write(11,*)
     enddo
     Res_l2 = sqrt(Res_l2/(Nx*Ny))
     !--- Solve ---!
     if (modulo(ExtIter,1000) .eq. 0) print*,'Starting Solve.....',Extiter
     convergence = -2
     liniter = 20
          
     Mat = Tot_Jac
     call seidel(0,nPoint*nVar,Mat,Tot_R(:),1.d0,Tot_Sol(:),Res(:),liniter,convergence)
     if (convergence .ne. 0) print*, 'Error in sys',convergence,ExtIter
     if (modulo(ExtIter,1000) .eq. 0) print*,'Finished sys, Res(u): ',Res_l2
     do iPoint =1,nPoint
        U(1,iPoint) = U_old(1,iPoint) + Tot_Sol((iPoint-1)*nVar+1)
        U(2,iPoint) = U_old(2,iPoint) + alfa*Tot_Sol((iPoint-1)*nVar+2)
        U(3,iPoint) = U_old(3,iPoint) + alfa*Tot_Sol((iPoint-1)*nVar+3)
     enddo
   
   else   
   !--- Time Integration (Explicit) ---!
     do iPoint=1,nPoint
     
        U_max = max(abs(U_old(2,iPoint)),abs(U_old(3,iPoint)))
        lambda_inv = U_max*dy + sqrt(U_max*dy*U_max*dy + (beta2/rho)*area*area) 
        if (abs(lambda_inv) .lt. 1.0e-6) lambda_inv = U_inf + (beta2/rho)*area
        !dt = CFL*Vol/lambda_inv
     
        U(1,iPoint) = U_old(1,iPoint) - R(1,iPoint)*dt/Vol
        U(2,iPoint) = U_old(2,iPoint) - R(2,iPoint)*dt/Vol
        U(3,iPoint) = U_old(3,iPoint) - R(3,iPoint)*dt/Vol
        
        Res_l2 = Res_l2 + (R(2,iPoint)**2.0)
     enddo
     Res_l2 = sqrt(Res_l2/(Nx*Ny))
     if (modulo(ExtIter,1000) .eq. 0) print*,'Res(u): ',Res_l2,' ExtIter: ',ExtIter
     
   endif
   
   !--- Convergence monitoring ---!
   
   
   !--- Output solution ---!
   if (modulo(ExtIter,100) .eq. 0) then
   
   open(unit=14,file='Centerline_test.txt',status='unknown')
     i = (Nx+1)/2
     do j=1,Ny
       iPoint = i + (j-1)*Nx
       write(14,*) y(i,j),U(1,iPoint),U(2,iPoint),U(3,iPoint),iPoint
     enddo
    close(14)
   endif

   if (modulo(ExtIter,1000) .eq. 0) then 
   open(unit=13,file='Residual_ac.txt',status='unknown')
     do iPoint = 1,nPoint
       write(13,*) iPoint
       write(13,*) U(1,iPoint),U(2,iPoint),U(3,iPoint),U_old(1,iPoint),U_old(2,iPoint),U_old(3,iPoint)
       write(13,*) Tot_Sol((iPoint-1)*nVar+1),Tot_Sol((iPoint-1)*nVar+2),Tot_Sol((iPoint-1)*nVar+3)
     enddo    
    close(13)
   endif
   
   
   !write(15,*)'---------------------Iteration----------------------------------',ExtIter
   do i=1,Nx
    do j=1,Ny
     !iPoint = i+(j-1)*Nx
     !write(15,*) GradU(1,1,i,j),GradU(1,2,i,j),Fc(1,i,j),Fv(1,i,j),U(1,i+(j-1)*Nx),R(1,iPoint),Fc(1,i,j)-Fv(1,i,j),i+(j-1)*Nx
     !write(15,*) GradU(2,1,i,j),GradU(2,2,i,j),Fc(2,i,j),Fv(2,i,j),U(2,i+(j-1)*Nx),R(2,iPoint),Fc(2,i,j)-Fv(2,i,j),i+(j-1)*Nx
     !write(15,*) GradU(3,1,i,j),GradU(3,2,i,j),Fc(3,i,j),Fv(3,i,j),U(3,i+(j-1)*Nx),R(3,iPoint),Fc(3,i,j)-Fv(3,i,j),i+(j-1)*Nx
     !write(15,*)
    enddo
   enddo
enddo !ExtIter

!--- Output solution ---!
open(unit=13,file='Residual_ac.txt',status='unknown')
open(unit=14,file='Centerline_test.txt',status='unknown')
i = (Nx+1)/2
do j=1,Ny
  iPoint = i + (j-1)*Nx
  write(14,*) y(i,j),U(1,iPoint),U(2,iPoint),U(3,iPoint),iPoint
enddo

do iPoint = 1,nPoint
  write(13,*) iPoint
  write(13,*) U(1,iPoint),U(2,iPoint),U(3,iPoint)
  write(13,*) Tot_Sol((iPoint-1)*nVar+1),Tot_Sol((iPoint-1)*nVar+2),Tot_Sol((iPoint-1)*nVar+3)
enddo



close(11)
close(12)
close(13)
close(14)
close(15)
close(16)
close(17)
end program iNS





subroutine seidel(crit,n,mat,b,omega,x,residu,iter,rc)
parameter(ITERMAX=1000)            ! Maximal number of iterations
parameter(ONE=1.d0,TWO=2.d0,ZERO=0.d0)
  integer crit, n, iter, rc
  REAL*8 mat(n,n),b(n),omega
  REAL*8 x(n),residu(n)
!*====================================================================*
!*                                                                    *
!*  seidel solves the linear system  mat * x = b  iteratively.        *
!*  Here  mat  is a nonsingular  n x n  matrix, b is the right hand   *
!*  side for the linear system and x is the solution.                 *
!*                                                                    *
!*  seidel uses the Gauss Seidel Method with relaxation for a given   *
!*  relaxation coefficient 0 < omega < 2.                             *
!*  If  omega = 1, the standard Gauss Seidel method (without          *
!*  relaxation) is performed.                                         *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Applications:                                                    *
!*   =============                                                    *
!*      Solve linear systems with nonsingular system matrices that    *
!*      satisfy one of the following criteria: row sum criterion,     *
!*      column sum criterion or the criterion of Schmidt and v. Mises.*
!*      Only if at least one of these criteria is satisfied for mat,  *
!*      convergence of the scheme is guaranteed [See BIBLI 11].       *
!*                                                                    *
!*====================================================================*
!*                                                                    *
!*   Input parameters:                                                *
!*   ================                                                 *
!*      crit     integer crit                                         *
!*               select criterion                                     *
!*               =1 : row sum criterion                               *
!*               =2 : column sum criterion                            *
!*               =3 : criterion of Schmidt-v.Mises                    *
!*               other : no check                                     *
!*      n        integer n ( n > 0 )                                  *
!*               size of mat, b and x                                 *
!*      mat      REAL*8   mat(n,n)                                    *
!*               Matrix of the liear system                           *
!*      b        REAL*8 b(n)                                          *
!*               Right hand side                                      *
!*      omega    REAL*8 omega; (0 < omega < 2)                        *
!*               Relaxation coefficient.                              *
!*      x        REAL*8  x(n)                                         *
!*               Starting vector for iteration                        *
!*                                                                    *
!*   Output parameters:                                               *
!*   ==================                                               *
!*      x        REAL*8  x(n)                                         *
!*               solution vector                                      *
!*      residu   REAL*8   residu(n)                                   *
!*               residual vector  b - mat * x; close to zero vector   *
!*      iter     integer iter                                         *
!*               Number of iterations performed                       *
!*      rc       integer return code                                  *
!*               =  0     solution has been found                     *
!*               =  1     n < 1  or omega <= 0 or omega >= 2          *
!*               =  2     improper mat or b or x (not used here)      *
!*               =  3     one diagonal element of mat vanishes        *
!*               =  4     Iteration number exceeded                   *
!*               = 11     column sum criterion violated               *
!*               = 12     row sum criterion violated                  *
!*               = 13     Schmidt-v.Mises criterion violated          *
!*                                                                    *
!*====================================================================*
  REAL*8 tmp, eps;

   rc = 0                       
   iter = 0 !Initialize iteration counter
  if (n<1.or.omega<=ZERO.or.omega>=TWO) then
    rc=1
    return
  end if

  eps = 1.d-10

  do i=1, n                       !transform mat so that all
                                          !diagonals equal 1
    if (mat(i,i) == ZERO) then
      rc=3
      return
    end if
    tmp = ONE / mat(i,i)
    do j=1, n
      mat(i,j)= mat(i,j)*tmp
    end do
    b(i) = b(i)*tmp               !adjust right hand side b
  
  end do


  !check convergence criteria
  if (crit==1) then
     do i = 1, n                  !row sum criterion
       tmp=ZERO
       do j=1,n
         tmp = tmp + dabs(mat(i,j))
       end do
       if (tmp >= TWO) then
         rc=11
         return
       end if 
     end do
  else if (crit==2) then  
     do j=1, n                    !column sum criterion
	   tmp=ZERO
       do i=1,n
         tmp = tmp + dabs(mat(i,j))
       end do
       if (tmp >= TWO) then
         rc=12
	 return
       end if
     end do
  else if (crit==3) then
     tmp=ZERO
     do i=1, n
       do j=1, n                  !criterion of Schmidt
         tmp = tmp + mat(i,j)**2  !von Mises
       end do
     end do
     tmp = DSQRT(tmp - ONE)
     if (tmp >= ONE) then
       rc=13
       return
     end if
  end if

  do i=1, n 
    residu(i) = x(i)              !store x in residu
  end do
  

  do while (iter <= ITERMAX)      !Begin iteration
  
    iter=iter+1

    do i=1, n
      tmp=b(i)
      do j=1, n
        tmp =  tmp - mat(i,j) * residu(j)
      end do 
      residu(i) = residu(i) + omega * tmp
    end do

    do i=1, n                     !check break-off criterion
      tmp = x(i) - residu(i)
      if (DABS (tmp) <= eps) then
        x(i) = residu(i)          !If rc = 0 at end of loop
        rc = 0                    !  -> stop iteration
      else
        do j=1, n 
          x(j) = residu(j)
        end do
        rc = 4
        goto 10
      end if
    end do
    if (rc == 0) goto 20          !solution found
10 end do                         !End iteration

20 do i=1, n                      !find residual vector
     tmp=b(i)
     do j=1, n
       tmp = tmp - mat(i,j) * x(j)
     end do
     residu(i) = tmp
   end do

  return

end

