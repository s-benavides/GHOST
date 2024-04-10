! Step 2 of Runge-Kutta for the MHD equations
! Computes the nonlinear terms and evolves the equations in dt/o
         CALL laplak2(ps,C3) ! Vorticity -W = -k^2 * psi
         CALL laplak2(az,C4) ! Current -J = -k^2 * a
         CALL poisson(ps,az,C5) ! -curl( u x a ) = [psi,a]
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL vectrans(az,C5,ext) ! a [ psi, a ]
         CALL poisson(ps,C3,ps) ! [psi, W]
         CALL poisson(az,C4,az) ! [a, J]
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL entrans(C1,ps,ext) ! psi [psi,W]. uuu
            CALL entrans_bbu(C1,az,ext) ! psi [a , J]. bbu
            CALL entrans_bub(C4,C5,ext) ! -J [psi,a]. bub+ubb

         rmp = 1.0_GP/real(o,kind=GP)
         DO i = ista,iend
            DO j = 1,n
            IF ((ka2(j,i).le.kmax).and.(ka2(j,i).ge.tiny)) THEN
               ps(j,i) = C1(j,i)+dt*((nu*ka2(j,i)**(hek-1) + hnu*ka2(j,i)**(-hok-1))*C3(j,i)+(az(j,i)-ps(j,i))    &
              /ka2(j,i)+fk(j,i))*rmp
               az(j,i) = C2(j,i)+dt*((mu*ka2(j,i)**(hem-1) + hmu*ka2(j,i)**(-hom-1))*C4(j,i)+C5(j,i)+mk(j,i))*rmp
            ELSE
               ps(j,i) = 0.0_GP
               az(j,i) = 0.0_GP
            ENDIF
            
            END DO
         END DO

