!=================================================================
! MODULES for FFTP v3
! Parallel Fast Fourier Transform in 2D
!
! 2007 Pablo D. Mininni.
!      National Center for Atmospheric Research.
!      e-mail: mininni@ucar.uba.ar 
!=================================================================

!=================================================================

  MODULE fftplans
!
! Set the variable ikind to:  4 in 32 bits machines
!                             8 in 64 bits machines
! Set the variable csize to:  8 if L1 cache is <= 64 kb
!                            16 if L1 cache is 128 kb
! The variable nstrip controls strip mining during the 
! transposition. Often set to 1.
!
      USE fprecision
      INCLUDE 'fftw3.f'
 
      INTEGER, PARAMETER  :: ikind = IKIND_
      INTEGER, PARAMETER  :: csize = CSIZE_
      INTEGER, PARAMETER  :: nstrip = NSTRIP_
      INTEGER, PARAMETER  :: FFTW_REAL_TO_COMPLEX = FFTW_FORWARD
      INTEGER, PARAMETER  :: FFTW_COMPLEX_TO_REAL = FFTW_BACKWARD
      DOUBLE PRECISION    :: ffttime = 0.0
      DOUBLE PRECISION    :: tratime = 0.0

      TYPE FFTPLAN
         COMPLEX(KIND=GP), DIMENSION (:,:), POINTER :: ccarr
         COMPLEX(KIND=GP), DIMENSION (:,:), POINTER :: carr
         REAL(KIND=GP), DIMENSION (:,:), POINTER    :: rarr
         INTEGER(kind=ikind) :: planr,planc
         INTEGER :: n
         INTEGER, DIMENSION (:), POINTER :: itype1, itype2
      END TYPE FFTPLAN
      SAVE

  END MODULE fftplans
!=================================================================

  MODULE mpivars
!     INCLUDE 'mpif.h'
      INTEGER, SAVE :: ista,iend
      INTEGER, SAVE :: jsta,jend
      INTEGER, SAVE :: ksta,kend
      INTEGER, SAVE :: nprocs,myrank
      INTEGER, SAVE :: ierr

  END MODULE mpivars
!=================================================================
