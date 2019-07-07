! Initial condition for the velocity.
! This file contains the expression used for the initial 
! velocity field. You can use temporary real arrays R1-R3 
! of size (1:n,1:n,ksta:kend) and temporary complex arrays 
! C1-C8 of size (n,n,ista:iend) to do intermediate 
! computations. The variable u0 should control the global 
! amplitude of the velocity, and variables vparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the velocity in spectral
! space should be stored in the arrays vx, vy, and vz.

! Superposition of ABC flows at large scales, and k^-m noise at 
! small scales
!     kdn    : minimum wave number
!     kup    : maximum wave number
!     vparam0: A amplitude
!     vparam1: B amplitude
!     vparam2: C amplitude
!     vparam3: spectral fall off of noise, fixed here
!     vparam4: noise amplitude as fraction of energy

      phase = pi/2.
      vparam3 = 3.5 !(N/64.)^3  NOT ! resolution times 9
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
            R1(i,j,k) = 0.
            DO ki = INT(kdn),INT(kup)
               R1(i,j,k) = R1(i,j,k)+(vparam1*COS(2*pi*ki*(float(j)-1)/     &
                          float(n)+phase)+vparam2*SIN(2*pi*ki*(float(k)-1)/ &
                          float(n)+phase))/float(ki)
            END DO
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
            R1(i,j,k) = 0.
            DO ki = INT(kdn),INT(kup)
               R1(i,j,k) = R1(i,j,k)+(vparam0*SIN(2*pi*ki*(float(i)-1)/     &
                          float(n)+phase)+vparam2*COS(2*pi*ki*(float(k)-1)/ &
                          float(n)+phase))/float(ki)
            END DO
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,R1,vy,MPI_COMM_WORLD)
      DO k = ksta,kend
         DO j = 1,n
            DO i = 1,n
            R1(i,j,k) = 0.
            DO ki = INT(kdn),INT(kup)
               R1(i,j,k) = R1(i,j,k)+(vparam0*COS(2*pi*ki*(float(i)-1)/     &
                          float(n)+phase)+vparam1*SIN(2*pi*ki*(float(j)-1)/ &
                          float(n)+phase))/float(ki)
            END DO
            END DO
         END DO
      END DO
      CALL fftp3d_real_to_complex(planrc,R1,vz,MPI_COMM_WORLD)

      IF (ista.eq.1) THEN
         C4(1,1,1) = 0.
         C5(1,1,1) = 0.
         C6(1,1,1) = 0.
         DO j = 2,n/2+1

            IF ((kn2(1,j,1).le.kmax/2).and.(kn2(1,j,1).gt.kup**2)) THEN
               dump = 3.e5/sqrt(kk2(1,j,1))**vparam3
               phase = 2*pi*randu(seed)
               C4(1,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C4(1,n-j+2,1) = conjg(C4(1,j,1))
               phase = 2*pi*randu(seed)
               C5(1,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C5(1,n-j+2,1) = conjg(C5(1,j,1))
               phase = 2*pi*randu(seed)
               C6(1,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C6(1,n-j+2,1) = conjg(C6(1,j,1))
            ELSE
               C4(1,j,1) = 0.
               C4(1,n-j+2,1) = 0.
               C5(1,j,1) = 0.
               C5(1,n-j+2,1) = 0.
               C6(1,j,1) = 0.
               C6(1,n-j+2,1) = 0.
            ENDIF

         END DO
         DO k = 2,n/2+1

            IF ((kn2(k,1,1).le.kmax/2).and.(kn2(k,1,1).gt.kup**2)) THEN
               dump = 3.e5/sqrt(kk2(k,1,1))**vparam3
               phase = 2*pi*randu(seed)
               C4(k,1,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C4(n-k+2,1,1) = conjg(C4(k,1,1))
               phase = 2*pi*randu(seed)
               C5(k,1,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C5(n-k+2,1,1) = conjg(C5(k,1,1))
               phase = 2*pi*randu(seed)
               C6(k,1,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C6(n-k+2,1,1) = conjg(C6(k,1,1))
            ELSE
               C4(k,1,1) = 0.
               C4(n-k+2,1,1) = 0.
               C5(k,1,1) = 0.
               C5(n-k+2,1,1) = 0.
               C6(k,1,1) = 0.
               C6(n-k+2,1,1) = 0.
            ENDIF

         END DO
         DO j = 2,n
            DO k = 2,n/2+1
     
            IF ((kn2(k,j,1).le.kmax/2).and.(kn2(k,j,1).gt.kup**2)) THEN
               dump = 3.e5/sqrt(kk2(k,j,1))**vparam3
               phase = 2*pi*randu(seed)
               C4(k,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C4(n-k+2,n-j+2,1) = conjg(C4(k,j,1))
               phase = 2*pi*randu(seed)
               C5(k,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C5(n-k+2,n-j+2,1) = conjg(C5(k,j,1))
               phase = 2*pi*randu(seed)
               C6(k,j,1) = (vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               C6(n-k+2,n-j+2,1) = conjg(C6(k,j,1))
            ELSE
               C4(k,j,1) = 0.
               C4(n-k+2,n-j+2,1) = 0.
               C5(k,j,1) = 0.
               C5(n-k+2,n-j+2,1) = 0.
               C6(k,j,1) = 0.
               C6(n-k+2,n-j+2,1) = 0.
            ENDIF

            END DO
         END DO
         DO i = 2,iend
            DO j = 1,n
               DO k = 1,n

               IF ((kn2(k,j,i).le.kmax/2).and.(kn2(k,j,i).gt.kup**2)) THEN
                  dump = 3.e5/sqrt(kk2(k,j,i))**vparam3
                  phase = 2*pi*randu(seed)
                  C4(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C5(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C6(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               ELSE
                  C4(k,j,i) = 0.
                  C5(k,j,i) = 0.
                  C6(k,j,i) = 0.
               ENDIF

               END DO
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               DO k = 1,n

               IF ((kn2(k,j,i).le.kmax/2).and.(kn2(k,j,i).gt.kup**2)) THEN
                  dump = 3.e5/sqrt(kk2(k,j,i))**vparam3
                  phase = 2*pi*randu(seed)
                  C4(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C5(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C6(k,j,i) = 2*(vparam2*COS(phase)+im*vparam1*SIN(phase))*dump
               ELSE
                  C4(k,j,i) = 0.
                  C5(k,j,i) = 0.
                  C6(k,j,i) = 0.
               ENDIF

               END DO
            END DO
        END DO
      ENDIF

      CALL rotor3(C5,C6,C7,1)
      CALL rotor3(C4,C6,C8,2)
      CALL rotor3(C4,C5,C6,3)
      CALL normalize(C7,C8,C6,sqrt(vparam4)*u0,1,MPI_COMM_WORLD)
      DO i = ista,iend
         DO j = 1,n
            DO k = 1,n
               vx(k,j,i) = vx(k,j,i)+C7(k,j,i)
               vy(k,j,i) = vy(k,j,i)+C8(k,j,i)
               vz(k,j,i) = vz(k,j,i)+C6(k,j,i)
            END DO
        END DO
      END DO
      CALL normalize(vx,vy,vz,u0,1,MPI_COMM_WORLD)

