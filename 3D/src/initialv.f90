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

! Superposition of ABC flows with k^(-4) spectrum
!     kdn    : minimum wave number
!     kup    : maximum wave number
!     vparam0: A amplitude
!     vparam1: B amplitude
!     vparam2: C amplitude

      IF ( (abs(Lx-Ly).gt.tiny).or.(abs(Lx-Lz).gt.tiny) ) THEN
        IF (myrank.eq.0) &
           PRINT *,'ABC initial conditions require Lx=Ly=Lz'
        STOP
      ENDIF

      DO k = ksta,kend
         DO j = 1,ny
            DO i = 1,nx

            R1(i,j,k) = 0.
            R2(i,j,k) = 0.
            R3(i,j,k) = 0.

            DO ki = INT(kdn),INT(kup)
               R1(i,j,k) = R1(i,j,k)+(vparam1*COS(2*pi*ki*(real(j,kind=GP)-1)/ &
                          real(ny,kind=GP))+vparam2*SIN(2*pi*ki*(real(k,kind=GP)-1)/ &
                          real(nz,kind=GP)))/ki**2
               R2(i,j,k) = R2(i,j,k)+(vparam0*SIN(2*pi*ki*(real(i,kind=GP)-1)/ &
                          real(nx,kind=GP))+vparam2*COS(2*pi*ki*(real(k,kind=GP)-1)/ &
                          real(nz,kind=GP)))/ki**2 
               R3(i,j,k) = R3(i,j,k)+(vparam0*COS(2*pi*ki*(real(i,kind=GP)-1)/ &
                          real(nx,kind=GP))+vparam1*SIN(2*pi*ki*(real(j,kind=GP)-1)/ &
                          real(ny,kind=GP)))/ki**2
            END DO

            END DO
         END DO
      END DO

      CALL fftp3d_real_to_complex(planrc,R1,vx,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R2,vy,MPI_COMM_WORLD)
      CALL fftp3d_real_to_complex(planrc,R3,vz,MPI_COMM_WORLD)
      CALL normalize(vx,vy,vz,u0,1,MPI_COMM_WORLD)
