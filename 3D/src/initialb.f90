! Initial conditions.
! This file contains the expression used for random initial
! conditions. You can use temporary real arrays
! R1-R3 of size (1:n,1:n,ksta:kend) and temporary complex
! arrays C1-C8 of size (n,n,ista:iend) to do intermediate
! computations. The variable u0 should control the global
! amplitude of the IC, and variables fparam0-9 can be
! used to control the amplitudes of individual terms. At the
! end, the three components of the velocity in spectral
! space should be stored in the arrays ax, ay, and az.

! Superposition of harmonic modes with random phases
!     kdn : minimum wave number
!     mkup : maximum wave number

! Notes:
! - Activates all scales < kf 
! - curl(C1,C2,C3) = (ax,ay,az)
! - Have: |C| = 1/kf**2 so that E_b(k) constant in the range.
! - Isotropic 

      IF (ista.eq.1) THEN
      ! Sets largest mode forcing to zero
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0. 
         ! cycles through (1,j,1), have to manually set a(-k) = cong(a(k))  
         DO j = 2,ny/2+1

            IF ((kk2(1,j,1).le.mkup**2).and.(kk2(1,j,1).ge.tiny)) THEN
               dump = 1./kk2(1,j,1)
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

            IF ((kk2(k,1,1).le.mkup**2).and.(kk2(k,1,1).ge.tiny)) THEN
               dump = 1./kk2(k,1,1)
               phase = 2*pi*randu(seed)
               C1(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,1,1) = conjg(C1(k,1,1))
               phase = 2*pi*randu(seed)
               C2(k,1,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,1,1) = conjg(C2(k,1,1))
               phase = 2*pi*randu(seed)
               C3(k,1,1) = (COS(phase)+im*SIN(phase))*dump
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
         ! cycles through (k,j,1)
         DO j = 2,ny
            DO k = 2,nz/2+1
     
            IF ((kk2(k,j,1).le.mkup**2).and.(kk2(k,j,1).ge.tiny)) THEN
               dump = 1./kk2(k,j,1)
               phase = 2*pi*randu(seed)
               C1(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C1(nz-k+2,ny-j+2,1) = conjg(C1(k,j,1))
               phase = 2*pi*randu(seed)
               C2(k,j,1) = (COS(phase)+im*SIN(phase))*dump
               C2(nz-k+2,ny-j+2,1) = conjg(C2(k,j,1))
               phase = 2*pi*randu(seed)
               C3(k,j,1) = (COS(phase)+im*SIN(phase))*dump
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
         ! finally cycles through (k,j,i) 
         DO i = 2,iend
            DO j = 1,ny
               DO k = 1,nz

               IF ((kk2(k,j,i).le.mkup**2).and.(kk2(k,j,i).ge.tiny)) THEN
                  dump = 1./kk2(k,j,i)
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = (COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = (COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = (COS(phase)+im*SIN(phase))*dump
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

               IF ((kk2(k,j,i).le.mkup**2).and.(kk2(k,j,i).ge.tiny)) THEN
                  dump = 1./kk2(k,j,i)
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = (COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = (COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = (COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF

               END DO
            END DO
        END DO
      ENDIF

      CALL rotor3(C2,C3,ax,1)
      CALL rotor3(C1,C3,ay,2)
      CALL rotor3(C1,C2,az,3)
      CALL normalize(ax,ay,az,a0,0,MPI_COMM_WORLD)
