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

! Random forcing (by Santiago Benavides)
! This function is meant to be used when rand = 1. It is
! made to choose one wavevector randomly, then making fx,
! fy, fz = 0 except at k_f, where it is random. This is in
! contrast to the other random forcing function which puts
! energy at all wavenumbers within the given shell at one
! time, which in certain cases is too strong.
! 
! If rand = 1, then this forcing function will give an
! energy injection rate of f0 at kup. 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

! Pick random wavenumber of length close to kup:
      IF (myrank.eq.0) THEN
        kfx = randu(seed)
        kfy = randu(seed)
        kfz = randu(seed)
        rmp = sqrt(kfx**2+kfy**2+kfz**2)
        IF (rmp>0.0) THEN
                kfx = floor(kup*kfx/rmp + 0.5)
                kfy = floor(kup*kfy/rmp + 0.5)
                kfz = floor(kup*kfz/rmp + 0.5)
        ELSE
                print*,"kfx=kfy=kfz=0!"
                kfx = 1.
                kfy = 1.
                kfz = 1.
                rmp = sqrt(kfx**2+kfy**2+kfz**2)
                kfx = floor(kup*kfx/rmp + 0.5)
                kfy = floor(kup*kfy/rmp + 0.5)
                kfz = floor(kup*kfz/rmp + 0.5)
        ENDIF
      ENDIF
        CALL MPI_BCAST(kfx,1,GC_REAL,0,MPI_COMM_WORLD,ierr)        
        CALL MPI_BCAST(kfy,1,GC_REAL,0,MPI_COMM_WORLD,ierr)        
        CALL MPI_BCAST(kfz,1,GC_REAL,0,MPI_COMM_WORLD,ierr)        

      IF (ista.eq.1) THEN
      ! Sets largest mode forcing to zero
         C1(1,1,1) = 0.
         C2(1,1,1) = 0.
         C3(1,1,1) = 0. 
         ! cycles through (1,j,1), have to manually set a(-k) = cong(a(k))  
         DO j = 2,ny/2+1
             IF ((kx(1).eq.kfx).and.((ky(j).eq.kfy).or.(ky(j)&
                        .eq.(-kfy))).and.(kz(1).eq.kfz)) THEN
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
             IF ((kx(1).eq.kfx).and.(ky(1).eq.kfy).and.&
                ((kz(k).eq.kfz).or.(kz(k).eq.(-kfz)))) THEN
               dump = 1./sqrt(kk2(k,1,1))
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
               IF ((kx(1).eq.kfx).and.((ky(j).eq.kfy).or.(ky(j).eq.&
               (-kfy))).and.((kz(k).eq.kfz).or.(kz(k).eq.(-kfz)))) THEN
               dump = 1./sqrt(kk2(k,j,1))
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

               IF (((kx(i).eq.kfx).or.(kx(i).eq.(-kfx))).and.&
                        (ky(j).eq.kfy).and.(kz(k).eq.kfz)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
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
                 
               IF (((kx(i).eq.kfx).or.(kx(i).eq.(-kfx))).and.&
                        (ky(j).eq.kfy).and.(kz(k).eq.kfz)) THEN
                  dump = 1./sqrt(kk2(k,j,i))
                  phase = 2*pi*randu(seed)
                  C1(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C2(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
                  phase = 2*pi*randu(seed)
                  C3(k,j,i) = 2*(COS(phase)+im*SIN(phase))*dump
               ELSE
                  C1(k,j,i) = 0.
                  C2(k,j,i) = 0.
                  C3(k,j,i) = 0.
               ENDIF
               END DO
            END DO
        END DO
      ENDIF

        
      CALL rotor3(C2,C3,fx,1)
      CALL rotor3(C1,C3,fy,2)
      CALL rotor3(C1,C2,fz,3) 
      IF (rand.eq.1) THEN
              dump = sqrt(2.0d0*f0)/sqrt(dt)  ! So that inj = f0
      ELSE 
              dump = f0
      ENDIF
      CALL normalize(fx,fy,fz,dump,1,MPI_COMM_WORLD)
