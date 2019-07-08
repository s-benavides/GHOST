! Step 2 of Runge-Kutta for the HD equations in a rotating frame
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL prodre3(vx,vy,vz,C4,C5,C6)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend               ! Coriolis force
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C4(k,j,i) = C4(k,j,i)-2*omegaz*vy(k,j,i)
                  C5(k,j,i) = C5(k,j,i)+2*omegaz*vx(k,j,i)
               END DO
            END DO
         END DO
         CALL nonlhd3(C4,C5,C6,C7,1)
         CALL nonlhd3(C4,C5,C6,C8,2)
         CALL nonlhd3(C4,C5,C6,C4,3)
         CALL diss(vx,vx,hek,hok,nu,hnu)   ! vx =-(k^(2*hek)+k^(-2*hok))*vx
         CALL diss(vy,vy,hek,hok,nu,hnu)
         CALL diss(vz,vz,hek,hok,nu,hnu)

         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            CALL entrans(C1,C2,C3,C7,C8,C4,ext,1,odir)
            CALL entpara(C1,C2,C3,C7,C8,C4,ext,1,odir)
            CALL entperp(C1,C2,C3,C7,C8,C4,ext,1,odir)
            CALL heltrans(C1,C2,C3,C7,C8,C4,ext,1,odir)
            CALL heltpara(C1,C2,C3,C7,C8,C4,ext,1,odir)
            CALL heltperp(C1,C2,C3,C7,C8,C4,ext,1,odir)
         ENDIF

         rmp = 1.0_GP/(real(o,kind=GP))
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.lt.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz
            IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1(k,j,i)+dt*(vx(k,j,i)+C7(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i)+dt*(vy(k,j,i)+C8(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i)+dt*(vz(k,j,i)+C4(k,j,i) &
              +fz(k,j,i))*rmp
            ELSE
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
            ENDIF
         END DO
         END DO
         END DO
