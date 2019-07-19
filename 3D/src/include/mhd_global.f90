! Global quantities computed in MHD runs

            CALL mhdcheck(vx,vy,vz,ax,ay,az,hek,hok,hem,hom,dump,dt,1,1,0)
            CALL cross(vx,vy,vz,fx,fy,fz,eps,1) ! v.f
            CALL cross(ax,ay,az,mx,my,mz,epm,0) ! b.curl(fb)
            CALL maxabs(vx,vy,vz,rmp,0)
            CALL maxabs(ax,ay,az,rmq,1)
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

            IF (myrank.eq.0) THEN
               OPEN(1,file='uf.txt',position='append')
               WRITE(1,FMT='(E25.18,E25.18,E25.18)') dump,tmr,tms
               CLOSE(1)                          ! u^2, b^2 at kf
               OPEN(1,file='injection.txt',position='append')
               WRITE(1,FMT='(E13.6,E22.14,E22.14)') dump,eps,epm
               CLOSE(1)
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,FMT='(E13.6,E13.6,E13.6)') dump,rmp,rmq
               CLOSE(1)
            ENDIF
