!> \file sub_time_integration.f90
!! \brief File with subroutines to perform Euler explicit and implicit time integrations.
subroutine explicit_euler

use global_vars

Sol(1:2,:) = 0.0
     do iPoint=1,nPoint
        U(1,iPoint) = U_old(1,iPoint) - alfa*R(1,iPoint)*dt_m/Vol(iPoint)
        U(2,iPoint) = U_old(2,iPoint) - alfa*R(2,iPoint)*dt_m/Vol(iPoint)
        
        Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Vol(iPoint)/(dt_m*alfa)
        Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Vol(iPoint)/(dt_m*alfa)
        
        Res_l2 = Res_l2 + R(1,iPoint)**2.0
        D(1,iPoint) = Vol(iPoint)/Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1)
        D(2,iPoint) = Vol(iPoint)/Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2)
     enddo
     Res_l2 = sqrt(Res_l2/nPoint)
     if (modulo(MIter,p_screen1) .eq. 0) print*,'Res: ',log10(Res_l2),MIter,ExtIter
end subroutine explicit_euler

subroutine implicit_euler

use global_vars

do iPoint=1,nPoint
        Tot_R((iPoint-1)*nVar+1) = -R(1,iPoint)
        Tot_R((iPoint-1)*nVar+2) = -R(2,iPoint)
       
        Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) = Tot_Jac((iPoint-1)*nVar+1,(iPoint-1)*nVar+1) + Vol(iPoint)/dt
        Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) = Tot_Jac((iPoint-1)*nVar+2,(iPoint-1)*nVar+2) + Vol(iPoint)/dt 
        
        Tot_Sol((iPoint-1)*nVar+1) = 0.0
        Tot_Sol((iPoint-1)*nVar+2) = 0.0
     enddo
     
     if (modulo(ExtIter,2) .eq. 0) print*,'Starting Solve.....',Extiter
     convergence = -2
     liniter = 1000
     !--- Solve ---!
     Mat = Tot_Jac
     call seidel(0,nPoint*nVar,Mat,Tot_R(:),1.d0,Tot_Sol(:),Res(:),liniter,convergence)
     if (convergence .ne. 0) print*, 'Error in mom',convergence,ExtIter
     if (modulo(ExtIter,5000) .eq. 0) print*,'Finished mom'
     do iPoint =1,nPoint
        U(1,iPoint) = U_old(1,iPoint) + Tot_Sol((iPoint-1)*nVar+1)
        U(2,iPoint) = U_old(2,iPoint) + Tot_Sol((iPoint-1)*nVar+2)
     enddo


end subroutine implicit_euler
