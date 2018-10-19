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
