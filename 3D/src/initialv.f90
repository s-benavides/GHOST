! Initial condition for the vector potential.
! This file contains the expression used for the initial 
! vector potential. You can use temporary real arrays R1-R3 
! of size (1:nx,1:ny,ksta:kend) and temporary complex arrays 
! C1-C8 of size (1:nz,1:ny,ista:iend) to do intermediate 
! computations. The variable u0 should control the global 
! amplitude of the forcing, and variables vparam0-9 can be  
! used to control the amplitudes of individual terms. At the 
! end, the three components of the potential in spectral 
! space should be stored in the arrays vx, ay, and az.

! Superposition of ABC flows at large scales, and k^-m noise at 
! small scales. These initial conditions require an isotropic box.
!     kdn   : minimum wave number (rounded to next integer)
!     kup   : maximum wave number (rounded to next integer)
!     vparam0: A amplitude
!     vparam1: B amplitude
!     vparam2: C amplitude
!     vparam3: index of spectral fall off of noise, fixed here
!     vparam4: noise level as fraction of mag. energy

      IF ( (abs(Lx-Ly).gt.tiny).or.(abs(Ly-Lz).gt.tiny) ) THEN
        IF (myrank.eq.0) &
           PRINT *,'ABC initial conditions require Lx=Ly=Lz'
        STOP
      ENDIF

      vparam3 = 4.5

      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
            R1(i,j,k) = 0.
            DO ki = INT(mkdn),INT(mkup)
               R1(i,j,k) = R1(i,j,k)+(vparam1*COS(2*pi*ki*(float(j)-1)/    &
                          float(ny))+vparam2*SIN(2*pi*ki*(float(k)-1)/      &
                          float(nz)))/float(ki**2)
            END DO
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
            R1(i,j,k) = 0.
            DO ki = INT(mkdn),INT(mkup)
               R1(i,j,k) = R1(i,j,k)+(vparam0*SIN(2*pi*ki*(float(i)-1)/    &
                          float(nx))+vparam2*COS(2*pi*ki*(float(k)-1)/      &
                          float(nz)))/float(ki**2)
            END DO
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,R1,vy,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx
            R1(i,j,k) = 0.
            DO ki = INT(mkdn),INT(mkup)
               R1(i,j,k) = R1(i,j,k)+(vparam0*COS(2*pi*ki*(float(i)-1)/    &
                          float(nx))+vparam1*SIN(2*pi*ki*(float(j)-1)/      &
                          float(ny)))/float(ki**2)
            END DO
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,R1,vz,MPI_COMM_WORLD)

      IF (ista.eq.1) THEN
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0.
         DO j = 2,ny/2+1

            IF ((kn2(1,j,1).le.kmax/2).and.(kn2(1,j,1).gt.mkup**2)) THEN
               dump = 3.e5/sqrt(kk2(1,j,1))**vparam3
               phase = 2*pi*randu(seed)
               C1(1,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C1(1,ny-j+2,1) = conjg(C1(1,j,1))
               phase = 2*pi*randu(seed)
               C2(1,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C2(1,ny-j+2,1) = conjg(C2(1,j,1))
               phase = 2*pi*randu(seed)
               C3(1,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C3(1,ny-j+2,1) = conjg(C3(1,j,1))
            ELSE
               C1(1,j,1) = 0.
               C1(1,ny-j+2,1) = 0.
               C2(1,j,1) = 0.
               C2(1,ny-j+2,1) = 0.
               C3(1,j,1) = 0.
               C3(1,ny-j+2,1) = 0.
            ENDIF

         END DO
         DO k = 2,nz/2+1

            IF ((kn2(k,1,1).le.kmax/2).and.(kn2(k,1,1).gt.mkup**2)) THEN
               dump = 3.e5/sqrt(kk2(k,1,1))**vparam3
               phase = 2*pi*randu(seed)
               C1(k,1,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C1(nz-k+2,1,1) = conjg(C1(k,1,1))
               phase = 2*pi*randu(seed)
               C2(k,1,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C2(nz-k+2,1,1) = conjg(C2(k,1,1))
               phase = 2*pi*randu(seed)
               C3(k,1,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C3(nz-k+2,1,1) = conjg(C3(k,1,1))
            ELSE
               C1(k,1,1) = 0.
               C1(nz-k+2,1,1) = 0.
               C2(k,1,1) = 0.
               C2(nz-k+2,1,1) = 0.
               C3(k,1,1) = 0.
               C3(nz-k+2,1,1) = 0.
            ENDIF

         END DO
         DO j = 2,ny
            DO k = 2,nz/2+1
     
            IF ((kn2(k,j,1).le.kmax/2).and.(kn2(k,j,1).gt.mkup**2)) THEN
               dump = 3.e5/sqrt(kk2(k,j,1))**vparam3
               phase = 2*pi*randu(seed)
               C1(k,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
               phase = 2*pi*randu(seed)
               C2(k,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
               phase = 2*pi*randu(seed)
               C3(k,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C3(nz-k+2,ny-j+2,1) = conjg(C3(k,j,1))
            ELSE
               C1(k,j,1) = 0.
               C1(nz-k+2,ny-j+2,1) = 0.
               C2(k,j,1) = 0.
               C2(nz-k+2,ny-j+2,1) = 0.
               C3(k,j,1) = 0.
               C3(nz-k+2,ny-j+2,1) = 0.
            ENDIF

            END DO
         END DO
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kn2(k,j,i).le.kmax/2).and.(kn2(k,j,i).gt.mkup**2)) THEN
                  dump = 3.e5/sqrt(kk2(k,j,i))**vparam3
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kn2(k,j,i).le.kmax/2).and.(kn2(k,j,i).gt.mkup**2)) THEN
                  dump = 3.e5/sqrt(kk2(k,j,i))**vparam3
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
        END DO
      ENDIF

      CALL rotor3(C2,C3,C7,1)
      CALL rotor3(C1,C3,C8,2)
      CALL rotor3(C1,C2,C9,3)
      CALL normalize(C7,C8,C9,sqrt(vparam4)*u0,0,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,ny
            DO k = 1,nz
               vx(k,j,i) = vx(k,j,i)+C7(k,j,i)
               vy(k,j,i) = vy(k,j,i)+C8(k,j,i)
               vz(k,j,i) = vz(k,j,i)+C9(k,j,i)
            END DO
        END DO
      END DO
      CALL normalize(vx,vy,vz,u0,0,MPI_COMM_WORLD)
