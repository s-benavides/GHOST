!=================================================================
! PSEUDOSPECTRAL subroutines
!
! Subroutines for computing spatial derivatives and nonlinear 
! terms in incompressible MHD and Hall-MHD equations in 2D 
! using a pseudo-spectral method. You should use the FFTPLANS 
! and MPIVARS modules (see the file 'fftp2D_mod.f90') in each 
! program that call any of the subroutines in this file. 
!
! NOTATION: index 'i' is 'x' 
!           index 'j' is 'y'
!
! 2003 Pablo D. Mininni.
!      Department of Physics, 
!      Facultad de Ciencias Exactas y Naturales.
!      Universidad de Buenos Aires.
!      e-mail: mininni@df.uba.ar 
!=================================================================

!*****************************************************************
      SUBROUTINE poissonb0(a,b,c,b0)
!-----------------------------------------------------------------
!
! Poisson bracket of the scalar fields A and B 
! in real space.
!
! Parameters
!     a : input matrix
!     b : input matrix
!     c : Poisson bracket {A,B} [output]
!     b0: amplitude of the uniform field in y
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE fft
      IMPLICIT NONE
      
      COMPLEX(KIND=GP), INTENT(IN),  DIMENSION(n,ista:iend) :: a,b
      COMPLEX(KIND=GP), INTENT(OUT), DIMENSION(n,ista:iend) :: c
      COMPLEX(KIND=GP), DIMENSION(n,ista:iend) :: c1,c2
      REAL(KIND=GP), DIMENSION(n,jsta:jend)    :: r1,r2,r3
      REAL(KIND=GP), INTENT(IN) :: b0
      REAL(KIND=GP) :: tmp
      INTEGER       :: i,j

!
! Computes dA/dx.dB/dy
!
      CALL derivk2(a,c1,1)
      CALL derivk2(b,c2,2)
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = r1(i,j)*r2(i,j)
         END DO
      END DO
!
! Computes (dA/dy+b0).dB/dx
!
      CALL derivk2(a,c1,2)
      CALL derivk2(b,c2,1)
      c2(1,1) = -b0*real(n,KIND=GP)**2
      CALL fftp2d_complex_to_real(plancr,c1,r1,MPI_COMM_WORLD)
      CALL fftp2d_complex_to_real(plancr,c2,r2,MPI_COMM_WORLD)
      tmp = 1.0_GP/real(n,kind=GP)**4
      DO j = jsta,jend
         DO i = 1,n
            r3(i,j) = (r3(i,j)-r1(i,j)*r2(i,j))*tmp
         END DO
      END DO

      CALL fftp2d_real_to_complex(planrc,r3,c,MPI_COMM_WORLD)

      RETURN
      END SUBROUTINE poissonb0

!*****************************************************************
      SUBROUTINE mhdcheck(a,b,c,d,t,nu,hnu,mu,hmu,hek,hok,hem,hom,kdn,kup,mkdn,mkup)
!-----------------------------------------------------------------
!
! Consistency check for the conservation of energy in MHD 2D
!
! Parameters
!     a  : streamfunction
!     b  : vector potential
!     c  : external kinetic force
!     d  : external magnetic force
!     t  : time
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      USE ali
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b,c,d
      DOUBLE PRECISION, INTENT(IN)  :: nu,hnu,mu,hmu
      DOUBLE PRECISION :: enk, denk, henk, injk, enkf  ! KE
      DOUBLE PRECISION :: enm, denm, henm, injm, enmf  ! ME
      DOUBLE PRECISION :: asq, dasq, hasq, inja, asqf  ! A^2
      DOUBLE PRECISION :: entot, tmp0, tmp1, tmp2, tmp3, tmp4, tmp5
      REAL(KIND=GP), INTENT(IN) :: t
      REAL(KIND=GP), INTENT(IN) :: kup,kdn
      REAL(KIND=GP), INTENT(IN) :: mkup,mkdn
      REAL(KIND=GP) :: tmq
      INTEGER, INTENT(IN) :: hek,hok ! Hyperviscosity powers
      INTEGER, INTENT(IN) :: hem,hom ! Hyperviscosity powers
      INTEGER       :: i,j

      tmq = 1.0_GP/real(n,kind=GP)**4
!
! Computes the mean energy, enstrophy, square 
! current, and square vector potential.
!
!! ENERGY
!       KINETIC 
      CALL energy(a,enk,1) ! Mean kinetic energy
      CALL energy2(a,denk,1+hek) ! Mean kinetic hyperdissipation
      denk = nu*denk
      CALL energy2(a,henk,1-hok) ! Mean kinetic hypodissipation
      henk = hnu*henk
!       MAGNETIC
      CALL energy(b,enm,1) ! Mean magnetic energy
      CALL energy2(b,denm,1+hem) ! Mean magnetic hyperdissipation
      denm = mu*denm
      CALL energy2(b,henm,1-hom) ! Mean magnetic hypodissipation
      henm = hmu*henm

      entot = enk + enm
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!     SQUARE VECTOR POTENTIAL
      CALL energy(b,asq,2)            ! A^2
      CALL energy(b,dasq, hem)         ! Diss
      dasq = mu*dasq
      CALL energy(b,hasq,-hom)         ! Hypodiss
      hasq = hmu*hasq      

!
! Computes, for KE, ME, and A^2 the (1) energy injection rate and (2) energy at forcing scale
!
      tmp0 = 0.0D0
      tmp1 = 0.0D0
      tmp2 = 0.0D0
      tmp3 = 0.0D0
      tmp4 = 0.0D0
      tmp5 = 0.0D0
      IF (ista.eq.1) THEN
         DO j = 1,n
            tmp0 = tmp0+real(d(j,1)*conjg(b(j,1)))*tmq
            tmp1 = tmp1+ka2(j,1)*real(c(j,1)*conjg(a(j,1)))*tmq
            tmp2 = tmp2+ka2(j,1)*real(d(j,1)*conjg(b(j,1)))*tmq
            ! Energy at forcing scale
            IF ((ka2(j,1).gt.(kdn**2/2.0)).and.(ka2(j,1).le.(kup**2*2.0))) THEN
                tmp4 = tmp4+ka2(j,1)*abs(a(j,1))**2*tmq
            ELSE IF ((ka2(j,1).gt.(mkdn**2/2.0)).and.(ka2(j,1).le.(mkup**2*2.0))) THEN
                tmp3 = tmp3+abs(b(j,1))**2*tmq
                tmp5 = tmp5+ka2(j,1)*abs(b(j,1))**2*tmq
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               tmp0 = tmp0+2*real(d(j,i)*conjg(b(j,i)))*tmq
               tmp1 = tmp1+2*ka2(j,i)*real(c(j,i)*conjg(a(j,i)))*tmq
               tmp2 = tmp2+2*ka2(j,i)*real(d(j,i)*conjg(b(j,i)))*tmq
            ! Energy at forcing scale
            IF ((ka2(j,i).gt.(kdn**2/2.0)).and.(ka2(j,i).le.(kup**2*2.0))) THEN
                tmp4 = tmp4+2*ka2(j,i)*abs(a(j,i))**2*tmq
            ELSE IF ((ka2(j,i).gt.(mkdn**2/2.0)).and.(ka2(j,i).le.(mkup**2*2.0))) THEN
                tmp3 = tmp3+2*abs(b(j,i))**2*tmq
                tmp5 = tmp5+2*ka2(j,i)*abs(b(j,i))**2*tmq
            ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               tmp0 = tmp0+2*real(d(j,i)*conjg(b(j,i)))*tmq
               tmp1 = tmp1+2*ka2(j,i)*real(c(j,i)*conjg(a(j,i)))*tmq
               tmp2 = tmp2+2*ka2(j,i)*real(d(j,i)*conjg(b(j,i)))*tmq
            ! Energy at forcing scale
            IF ((ka2(j,i).gt.(kdn**2/2.0)).and.(ka2(j,i).le.(kup**2*2.0))) THEN
                tmp4 = tmp4+2*ka2(j,i)*abs(a(j,i))**2*tmq
            ELSE IF ((ka2(j,i).gt.(mkdn**2/2.0)).and.(ka2(j,i).le.(mkup**2*2.0))) THEN
                tmp3 = tmp3+2*abs(b(j,i))**2*tmq
                tmp5 = tmp5+2*ka2(j,i)*abs(b(j,i))**2*tmq
            ENDIF
            END DO
         END DO
      ENDIF
      CALL MPI_REDUCE(tmp0,inja,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tmp1,injk,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tmp2,injm,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tmp3,asqf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tmp4,enkf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      CALL MPI_REDUCE(tmp5,enmf,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)

!!
!! Computes the cross correlation between
!! velocity and magnetic fields
!!
!      tmp = 0.
!      IF (ista.eq.1) THEN
!         DO j = 1,n
!            tmp = tmp+ka2(j,1)*real(b(j,1)*conjg(a(j,1)))*tmq
!         END DO
!         DO i = 2,iend
!            DO j = 1,n
!               tmp = tmp+2*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmq
!            END DO
!         END DO
!      ELSE
!         DO i = ista,iend
!            DO j = 1,n
!               tmp = tmp+2*ka2(j,i)*real(b(j,i)*conjg(a(j,i)))*tmq
!            END DO
!         END DO
!      ENDIF
!      CALL MPI_REDUCE(tmp,udb,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
!                      MPI_COMM_WORLD,ierr)
!
! Creates external files to store the results
!
      IF (myrank.eq.0) THEN
         OPEN(1,file='energy_k.txt',position='append')
         WRITE(1,10) t,enk,denk,henk,injk,enkf
   10    FORMAT( E26.18,E26.18,E26.18,E26.18,E26.18,E26.18 )
         CLOSE(1)
         OPEN(1,file='energy_m.txt',position='append')
         WRITE(1,10) t,enm,denm,henm,injm,enmf
         CLOSE(1)
         OPEN(1,file='energy_a.txt',position='append')
         WRITE(1,10) t,asq,dasq,hasq,inja,asqf
         CLOSE(1)
      ENDIF      

      RETURN
      END SUBROUTINE mhdcheck

!*****************************************************************
      SUBROUTINE pmspectrum(ps,a,nmb,kin)
!-----------------------------------------------------------------
!
! Computes the energy power spectrum in 2D for E+/-. 
! The output is written to a file by the first node.
!
! Parameters
!     ps : streamfunction 
!     a  : vector potential
!     nmb: the extension used when writting the file
!     kin: =0 computes the E+ spectrum
!          =1 computes the E- spectrum
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: ps,a
      DOUBLE PRECISION, DIMENSION(n/2+1)  :: Ek,Ektot
      REAL(KIND=GP)       :: q, sgn, tmp
      INTEGER, INTENT(IN) :: kin
      INTEGER             :: kmn
      INTEGER             :: i,j
      CHARACTER(len=*), INTENT(IN) :: nmb

      sgn = 2.0_GP
      IF (kin.eq.1) THEN
        sgn = -2.0_GP
      ENDIF
!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0_GP
      END DO
!
! Computes the energy spectrum
!
      tmp = 1.0_GP/real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.5)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
               q = abs(ps(j,1))**2+abs(a(j,1))**2+     &
                   sgn*real(ps(j,1)*conjg(a(j,1)))
               Ek(kmn) = Ek(kmn)+ka2(j,1)*q*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  q = abs(ps(j,i))**2+abs(a(j,i))**2+  & 
                      sgn*real(ps(j,i)*conjg(a(j,i)))
                  Ek(kmn) = Ek(kmn)+2*ka2(j,i)*q*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
                  q = abs(ps(j,i))**2+abs(a(j,i))**2+  &
                      sgn*real(ps(j,i)*conjg(a(j,i)))
                  Ek(kmn) = Ek(kmn)+2*ka2(j,i)*q*tmp
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         IF (kin.eq.0) THEN
            OPEN(1,file='epspectrum.' // nmb // '.txt')
         ELSE 
            OPEN(1,file='emspectrum.' // nmb // '.txt')
         ENDIF
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE pmspectrum

!*****************************************************************
      SUBROUTINE vectrans(a,b,nmb)
!-----------------------------------------------------------------
!
! Computes the square vector potential transfer in 
! Fourier space in 2D MHD. The output is written 
! to a file by the first node.
!
! Parameters
!     a  : vector potential
!     b  : Poisson bracket of the streamfunction and a
!     nmb: the extension used when writting the file
!
      USE fprecision
      USE commtypes
      USE mpivars
      USE filefmt
      USE grid
      USE kes
      IMPLICIT NONE

      COMPLEX(KIND=GP), INTENT(IN), DIMENSION(n,ista:iend) :: a,b
      DOUBLE PRECISION, DIMENSION(n/2+1) :: Ek,Ektot
      REAL(KIND=GP)       :: tmp
      INTEGER             :: kmn
      INTEGER             :: i,j
      CHARACTER(len=*), INTENT(IN) :: nmb

!
! Sets Ek to zero
!
      DO i = 1,n/2+1
         Ek(i) = 0.0_GP
      END DO
!
! Computes the square vector potential flux
!
      tmp = 1.0_GP/real(n,kind=GP)**4
      IF (ista.eq.1) THEN
         DO j = 1,n
            kmn = int(sqrt(ka2(j,1))+.501)
            IF ((kmn.gt.0).and.(kmn.le.n/2+1)) THEN
              Ek(kmn) = Ek(kmn)+real(b(j,1)*conjg(a(j,1)))*tmp
            ENDIF
         END DO
         DO i = 2,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+2*real(b(j,i)*conjg(a(j,i)))*tmp
               ENDIF
            END DO
         END DO
      ELSE
         DO i = ista,iend
            DO j = 1,n
               kmn = int(sqrt(ka2(j,i))+.5)
               IF (kmn.le.n/2+1) THEN
                  Ek(kmn) = Ek(kmn)+2*real(b(j,i)*conjg(a(j,i)))*tmp
               ENDIF
            END DO
         END DO
      ENDIF
!
! Computes the reduction between nodes
! and exports the result to a file
!
      CALL MPI_REDUCE(Ek,Ektot,n/2+1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
      IF (myrank.eq.0) THEN
         OPEN(1,file='vectransf.' // nmb // '.txt')
         WRITE(1,20) Ektot
   20    FORMAT( E23.15 ) 
         CLOSE(1)
      ENDIF

      RETURN
      END SUBROUTINE vectrans
