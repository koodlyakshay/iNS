!> \file sub_time_integration.f90
!! \brief File with subroutines to perform Euler explicit and implicit time integrations.
subroutine set_time_step

use global_vars
use output_vars
use flow_parmaters

implicit none

integer           :: iPoint, nPoint
real              :: U_ref,lambda,dt_i,dt_v


do iPoint=1,Nx*Ny
  U_ref = U_inf + sum(abs(U_old(:,iPoint)))
  lambda = U_ref*sqrt(Vol(iPoint))
  dt_i = CFL*Vol(iPoint)/lambda
  lambda = (4.0/3.0)*mu*Vol(iPoint)/rho
  dt_v = CFL*0.25*Vol(iPoint)*Vol(iPoint)/lambda
  dt_m(iPoint) = min(dt_v,dt_i)
enddo

endsubroutine set_time_step

subroutine explicit_euler(extiter)

use global_vars
use output_vars
use flow_parmaters

implicit none

integer           :: iPoint, nPoint
integer           :: ExtIter, miter, nmiter
real              :: Res_l2
real              :: Tot_Sol(Nx*Ny*nVar), Tot_R(Nx*Ny*nVar)

  nPoint = Nx*Ny
  Sol(1:2,:) = 0.0
  Res_l2 = 0.0
  
  call set_time_step


  do iPoint=1,nPoint
    U(1,iPoint) = U_old(1,iPoint) - alfa*R(1,iPoint)*dt_m(iPoint)/Vol(iPoint)
    U(2,iPoint) = U_old(2,iPoint) - alfa*R(2,iPoint)*dt_m(iPoint)/Vol(iPoint)
    
    Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(dt_m(iPoint)*alfa)
    Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(dt_m(iPoint)*alfa)
    
    Res_l2 = Res_l2 + R(1,iPoint)**2.0
    D(1,iPoint) = Vol(iPoint)/Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
    D(2,iPoint) = Vol(iPoint)/Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)
  enddo

  Res_l2 = sqrt(Res_l2/nPoint)
  if ((modulo(ExtIter,p_out) .eq. 0) .and. (MIter .eq. nMIter)) print*,'Res: ',log10(Res_l2),MIter,ExtIter

end subroutine explicit_euler

subroutine implicit_euler(extiter)

use global_vars
use output_vars
use flow_parmaters

implicit none

integer           :: iPoint, nPoint, nsize
integer           :: ExtIter, miter, nmiter,ipiv(Nx*Ny*nVar)
integer           :: convergence, liniter
real              :: Res_l2
real              :: Tot_Sol(Nx*Ny*nVar), Tot_R(Nx*Ny*nVar)

  nPoint = Nx*Ny
  Res_l2 = 0.0
  
  Tot_Jac = Convec_Jac + Visc_Jac
  
  call set_time_step
  
!  !$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(iPoint), SHARED(nVar,Tot_Sol,Tot_R,R,Tot_Jac,D,Vol,dt_m),REDUCTION(+:Res_l2)
  do iPoint=1,nPoint
    Tot_R((iPoint-1)*nVar+1) = -R(1,iPoint)
    Tot_R((iPoint-1)*nVar+2) = -R(2,iPoint)

    Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + Vol(iPoint)/dt_m(iPoint)
    Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + Vol(iPoint)/dt_m(iPoint)

    Tot_Sol((iPoint-1)*nVar+1) = 0.0
    Tot_Sol((iPoint-1)*nVar+2) = 0.0

    D(1,iPoint) = Vol(iPoint)/Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
    D(2,iPoint) = Vol(iPoint)/Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)
    
    Res_l2 = Res_l2 + R(1,iPoint)**2.0
   
  enddo
!  !$OMP END PARALLEL DO
     
   convergence = -2
   liniter = 1000
   !--- Solve ---!
   Mat = Tot_Jac
   call seidel(0,nPoint*nVar,Mat(1:npoint*nvar,1:npoint*nvar),Tot_R(1:npoint*nvar),1.0,Tot_Sol(1:npoint*nvar), &
                        Res(1:npoint*nvar),liniter,convergence)
   if ((convergence .ne. 0) .and. (convergence .ne. 4)) print*, 'Error in mom',convergence,ExtIter
   
   !nSize = nPoint*nVar
   !call dgesv(nSize,1,Mat(1:nSize,1:nSize),nSize,ipiv,Tot_R(1:nSize),nSize,convergence)
   
   !if ((convergence .ne. 0)) print*, 'Error in mom',convergence,ExtIter
   !Tot_Sol = Tot_R
   
   !if (modulo(ExtIter,5000) .eq. 0) print*,'Finished mom'
   Res_l2 = sqrt(Res_l2/nPoint)
   if ((modulo(ExtIter,p_out) .eq. 0)) print*,'Res: ',log10(Res_l2),ExtIter
!   !$OMP PARALLEL DO DEFAULT(NONE), PRIVATE(iPoint), SHARED(nVar,U,U_old,Tot_Sol)
   do iPoint =1,nPoint
      U(1,iPoint) = U_old(1,iPoint) + Tot_Sol((iPoint-1)*nVar+1)
      U(2,iPoint) = U_old(2,iPoint) + Tot_Sol((iPoint-1)*nVar+2)
   enddo
!   !$OMP END PARALLEL DO
   
end subroutine implicit_euler
