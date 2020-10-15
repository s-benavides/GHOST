! Global quantities computed in MHD runs
            CALL mhdcheck(vx,vy,vz,ax,ay,az,hek,hok,hem,hom,dump,dt,1,1,0)
            CALL cross(vx,vy,vz,fx,fy,fz,eps,1)
            CALL cross(ax,ay,az,mx,my,mz,epm,0)
            CALL maxabs(vx,vy,vz,rmp,0)
            CALL maxabs(ax,ay,az,rmq,1)
!
! Calculates the energy perpendicular to B and Omega separately.
!

      tmp = 0.0D0
      tmq = 0.0D0                            !dpe,dpa
      CALL energy_arbdir(vx,vy,vz,bx0,by0,bz0,tmp,tmq)
      tmr = 0.0D0
      tms = 0.0D0
      CALL energy_arbdir(vx,vy,vz,omegax,omegay,omegaz,tmr,tms)

      IF (myrank.eq.0) THEN
               OPEN(1,file='eperp_epara.txt',position='append')
          WRITE(1,FMT='(E25.18,E25.18,E25.18,E25.18,E25.18)') dump,tmp,tmq,tmr,tms
               CLOSE(1)                 ! perp2B,para2B,perp2omega,para2omega
       ENDIF

! Now for magnetic energy
      CALL rotor3(ay,az,c1,1) ! computing bx
      CALL rotor3(ax,az,c2,2)
      CALL rotor3(ax,ay,c3,3)

      tmp = 0.0D0
      tmq = 0.0D0                            !dpe,dpa
      CALL energy_arbdir(c1,c2,c3,bx0,by0,bz0,tmp,tmq)
      tmr = 0.0D0
      tms = 0.0D0
      CALL energy_arbdir(c1,c2,c3,omegax,omegay,omegaz,tmr,tms)

      IF (myrank.eq.0) THEN
               OPEN(1,file='ebperp_ebpara.txt',position='append')
          WRITE(1,FMT='(E25.18,E25.18,E25.18,E25.18,E25.18)') dump,tmp,tmq,tmr,tms
               CLOSE(1)      !magnetic energy: perp2B,para2B,perp2omega,para2omega
       ENDIF


!
! Computes |u|^2 at k_f.
!
      tmp = 0.0D0
      tmq = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

         IF (ista.eq.1) THEN
            DO j = 1,ny
               DO k = 1,nz
                  ki = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((ki.gt.(kdn/sqrt(2.0))).and.(ki.le.(kup*sqrt(2.0)))) THEN
                  tmp = tmp+(abs(vx(k,j,1))**2+abs(vy(k,j,1))**2+ &
                         abs(vz(k,j,1))**2)*tmq
                  ENDIF
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,ny
                  DO k = 1,nz
                  ki = int(sqrt(kk2(k,j,i))/Dkk+.501)
                  IF ((ki.gt.(kdn/sqrt(2.0))).and.(ki.le.(kup*sqrt(2.0)))) THEN
                     tmp = tmp+2*(abs(vx(k,j,i))**2+abs(vy(k,j,i))**2+ &
                            abs(vz(k,j,i))**2)*tmq
                  ENDIF
                  END DO
               END DO
            END DO
          ELSE
            DO i = ista,iend
               DO j = 1,ny
                  DO k = 1,nz
                  ki = int(sqrt(kk2(k,j,i))/Dkk+.501)
                  IF ((ki.gt.(kdn/sqrt(2.0))).and.(ki.le.(kup*sqrt(2.0)))) THEN
                     tmp = tmp + 2*(abs(vx(k,j,i))**2+abs(vy(k,j,i))**2+ &
                            abs(vz(k,j,i))**2)*tmq
                   ENDIF
                END DO
               END DO
            END DO
          ENDIF


     CALL MPI_REDUCE(tmp,tmr,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Computes |b|^2 at k_f.
!
      CALL rotor3(ay,az,c1,1) ! computing bx
      CALL rotor3(ax,az,c2,2)
      CALL rotor3(ax,ay,c3,3)
      tmp = 0.0D0
      tmq = 1.0_GP/ &
            (real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP))**2

         IF (ista.eq.1) THEN
            DO j = 1,ny
               DO k = 1,nz
                  ki = int(sqrt(kk2(k,j,1))/Dkk+.501)
                  IF ((ki.gt.(kdn/sqrt(2.0))).and.(ki.le.(kup*sqrt(2.0)))) THEN
                  tmp = tmp+(abs(c1(k,j,1))**2+abs(c2(k,j,1))**2+ &
                         abs(c3(k,j,1))**2)*tmq
                  ENDIF
               END DO
            END DO
            DO i = 2,iend
               DO j = 1,ny
                  DO k = 1,nz
                  ki = int(sqrt(kk2(k,j,i))/Dkk+.501)
                  IF ((ki.gt.(kdn/sqrt(2.0))).and.(ki.le.(kup*sqrt(2.0)))) THEN
                     tmp = tmp+2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                            abs(c3(k,j,i))**2)*tmq
                  ENDIF
                  END DO
               END DO
            END DO
          ELSE
            DO i = ista,iend
               DO j = 1,ny
                  DO k = 1,nz
                  ki = int(sqrt(kk2(k,j,i))/Dkk+.501)
                  IF ((ki.gt.(kdn/sqrt(2.0))).and.(ki.le.(kup*sqrt(2.0)))) THEN
                     tmp = tmp + 2*(abs(c1(k,j,i))**2+abs(c2(k,j,i))**2+ &
                            abs(c3(k,j,i))**2)*tmq
                   ENDIF
                END DO
               END DO
            END DO
          ENDIF

     CALL MPI_REDUCE(tmp,tms,1,MPI_DOUBLE_PRECISION,MPI_SUM,0, &
                      MPI_COMM_WORLD,ierr)
!
! Computes the hybrid helicity:
!           
                tmp = sqrt((bx0)**2+(by0)**2+(bz0)**2) 
            IF ((tmp>0).and.(omegaz>0)) THEN
                CALL helicity(ax,ay,az,helm)  !computes magnetic helicity
                CALL cross(vx,vy,vz,c1,c2,c3,chel,1)   ! computes cross helicity
                tmp = omegaz/tmp      ! computes 'sigma' = Omega/B0
                tmp = chel - tmp*helm   ! computes the hybrid helicity. 
              IF (myrank.eq.0) THEN
               OPEN(1,file='hybridhel.txt',position='append')
               WRITE(1,FMT='(E25.18,E25.18)') dump,tmp
               CLOSE(1)
              ENDIF
            ENDIF 

!
! Computes b.grad(v) and B_0.grad(v)
!
      CALL rotor3(ay,az,c1,1) ! computing bx
      CALL rotor3(ax,az,c2,2)
      CALL rotor3(ax,ay,c3,3)
      CALL vector3(vx,vy,vz,c1,c2,c3,c4,c5,c6) ! v x b
      CALL cross(ax,ay,az,c4,c5,c6,tmv,0) ! < b . curl(v x b) >
      ! Manually calculate v x B_0 since B_0 is uniform
      c1 = vy*bz0 - vz*by0
      c2 = vz*bx0 - vx*bz0
      c3 = vx*by0 - vy*by0  
      CALL cross(ax,ay,az,c1,c2,c3,tmp,0) ! < b . curl(v x B_0) >

!!!!!!!!!!!!!!!!!!!!!

            IF (myrank.eq.0) THEN
               OPEN(1,file='uf.txt',position='append')
               WRITE(1,FMT='(E25.18,E25.18,E25.18)') dump,tmr,tms
               CLOSE(1)                          ! t,  u^2, b^2 at kf
               OPEN(1,file='b0_bal.txt',position='append')
               WRITE(1,FMT='(E25.18,E26.18,E26.18)') dump,tmv,tmp
               CLOSE(1)                          ! t , bal_b, bal_b0
               OPEN(1,file='injection.txt',position='append')
               WRITE(1,FMT='(E13.6,E22.14,E22.14)') dump,eps,epm
               CLOSE(1)
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,FMT='(E13.6,E13.6,E13.6)') dump,rmp,rmq
               CLOSE(1)
            ENDIF
