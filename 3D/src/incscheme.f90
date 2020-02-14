! User-defined forcing scheme.
! A forcing scheme can be implemented in this file, which
! is used by the code when 'rand=3' is set in 'parameter.dat'.
! This scheme is executed every 'fstep' time steps, or every
! 'cort' time interval, if using CFL condition. See the
! folder 'examples' for an example. If not needed, this file
! can be left empty.
! Forcing arrays are complex (in Fourier space) of size
! (n,n,ista:iend) and are called:
!       (fx,fy,fz)   for the velocity
!       (mx,my,mz)   for the e.m.f. (magnetic field)
!       (fs,fs1,...) for scalar fields
!       (fre,fim)    for quantum solvers
! You can use temporary real arrays R1-R3 of size
! (1:n,1:n,ksta:kend) and temporary complex arrays C1-C8 of
! size (n,n,ista:iend) to do intermediate computations,
! and two real arrays Faux1 and Faux2 of size (10) to store
! information of the history of the forcing if needed.

! ----------  Santi's Comments  -----------------------------

! This forcing is the same as initialfv.f90_random. If it 
! is included using rand = 3, then the phases are constantly
! being changed, but the amplitude is just f0, no longer
! f0/sqrt(dt), which is used if rand.eq.1.

! ----------------------------------------------------------

! External mechanical forcing.
! This file contains the expression used for the external
! mechanical forcing. You can use temporary real arrays
! R1-R3 of size (1:n,1:n,ksta:kend) and temporary complex
! arrays C1-C8 of size (n,n,ista:iend) to do intermediate
! computations. The variable f0 should control the global
! amplitude of the forcing, and variables fparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the forcing in spectral
! space should be stored in the arrays fx, fy, and fz.

! Superposition of harmonic modes with random phases
!     kdn : minimum wave number
!     kup : maximum wave number

! Notes:
! - curl(C1,C2,C3) = (fx,fy,fz)
! - Forces ONLY 2D modes
! - Have: |C| = 1/kf so that |F| = f0 in the end.
! - If I use this with f0 = f0/sqrt(dt), will this give me an energy injection of 
!   f0^2 as I would guess?

      IF (ista.eq.1) THEN
      ! Sets largest mode forcing to zero
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0. 
         ! cycles through (1,j,1), have to manually set a(-k) = cong(a(k))  
         DO j = 2,ny/2+1

            IF ((kk2(1,j,1).le.kup**2).and.(kk2(1,j,1).ge.kdn**2)) THEN
               dump = 1./sqrt(kk2(1,j,1))
               phase = 2*pi*randu(seed)
               C1(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(1,ny-j+2,1) = conjg(C1(1,j,1))
               phase = 2*pi*randu(seed)
               C2(1,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(1,ny-j+2,1) = conjg(C2(1,j,1))
               phase = 2*pi*randu(seed)
               C3(1,j,1) = (COS(phase)+im*SIN(phase))*dump
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
         ! cycles through (k,1,1)
         DO k = 2,nz/2+1

               C1(k,1,1) = 0. 
               C1(nz-k+2,1,1) = 0.
               C2(k,1,1) = 0. 
               C2(nz-k+2,1,1) = 0.
               C3(k,1,1) = 0.
               C3(nz-k+2,1,1) = 0.

         END DO
         ! cycles through (k,j,1)
         DO j = 2,ny
            DO k = 2,nz/2+1
     
               C1(k,j,1) = 0.
               C1(nz-k+2,ny-j+2,1) = 0.
               C2(k,j,1) = 0.
               C2(nz-k+2,ny-j+2,1) = 0.
               C3(k,j,1) = 0.
               C3(nz-k+2,ny-j+2,1) = 0.

            END DO
         END DO
         DO i=2,iend
            DO j = 1,ny
         ! cycles through (k.ne.1,j,i)
             DO k = 2,nz
               C1(k,j,i) = 0.
               C2(k,j,i) = 0.
               C3(k,j,i) = 0.
             END DO

         ! finally cycles through (1,j,i)

               IF ((kk2(1,j,i).le.kup**2).and.(kk2(1,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(1,j,i))
                  phase = 2*pi*randu(seed)
                  C1(1,j,i) = (COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(1,j,i) = (COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(1,j,i) = (COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(1,j,i) = 0.
                  C2(1,j,i) = 0.
                  C3(1,j,i) = 0.
               ENDIF

            END DO
         END DO
      ELSE
        DO i = ista,iend

            DO j = 1,ny
         ! cycles through (k.ne.1,j,i)
             DO k = 2,nz
               C1(k,j,i) = 0.
               C2(k,j,i) = 0.
               C3(k,j,i) = 0.
             END DO

         ! finally cycles through (1,j,i)

               IF ((kk2(1,j,i).le.kup**2).and.(kk2(1,j,i).ge.kdn**2)) THEN
                  dump = 1./sqrt(kk2(1,j,i))
                  phase = 2*pi*randu(seed)
                  C1(1,j,i) = (COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(1,j,i) = (COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(1,j,i) = (COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(1,j,i) = 0.
                  C2(1,j,i) = 0.
                  C3(1,j,i) = 0.
               ENDIF

            END DO

        END DO
      ENDIF

      CALL rotor3(C2,C3,fx,1)
      CALL rotor3(C1,C3,fy,2)
      CALL rotor3(C1,C2,fz,3)
      IF (rand.eq.0) THEN
              dump = f0
      ELSE 
              dump = sqrt(2.0d0*f0)/sqrt(dt)  ! So that inj = f0
      ENDIF
      CALL normalize(fx,fy,fz,dump,1,MPI_COMM_WORLD)
