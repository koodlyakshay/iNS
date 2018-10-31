subroutine allocate_vars

use global_vars

allocate(x(Nx,Ny),y(Nx,Ny),P(Nx,Ny),Vol(nPoint))
allocate(U(2,Nx*Ny),U_old(2,Nx*Ny),Sol(3,Nx*Ny),V(3,Nx,Ny),Fc(2,Nx,Ny),Fv(2,Nx,Ny),P_Correc(nPoint))
allocate(R(3,Nx*Ny),Jac(nPoint,nPoint),GradU(3,2,Nx,Ny),Res(Nx*Ny*nVar),Mass(Nx*Ny),D(2,Nx*Ny))
allocate(Tot_Jac(Nx*Ny*nVar,Nx*Ny*nVar),Tot_R(Nx*Ny*nVar),Tot_Sol(Nx*Ny*nVar),GradPc(2,Nx,Ny))

end subroutine allocate_vars
