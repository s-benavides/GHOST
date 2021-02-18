! Global quantities computed in HD runs in a rotating frame

            CALL qmhdcheck(vx,vy,vz,fx,fy,fz,hek,hok,NNx,NNz,dump,dt,1,0)
            CALL maxabs(vx,vy,vz,rmp,0)
!
! Calculates the energy perpendicular to B and Omega separately.
!
       
      tmp = 0.0D0
      tmq = 0.0D0                            !dpe,dpa
      CALL energy_arbdir(vx,vy,vz,NNx,0.0D0,NNz,tmp,tmq)
      tmr = 0.0D0
      tmt = 0.0D0
      CALL energy_arbdir(vx,vy,vz,omegax,omegay,omegaz,tmr,tmt)
        
      IF (myrank.eq.0) THEN
         OPEN(1,file='eperp_epara.txt',position='append')
           WRITE(1,FMT='(E25.18,E25.18,E25.18,E25.18,E25.18)') dump,tmp,tmq,tmr,tmt
         CLOSE(1)                 ! perp2B,para2B,perp2omega,para2omega
      ENDIF

!
! Component energies!
!

      CALL energy(vx,0.0D0*vy,0.0D0*vz,tmp,1) ! <|vx|^2>
      CALL energy(0.0D0*vx,vy,0.0D0*vz,tmq,1) ! <|vy|^2>
      CALL energy(0.0D0*vx,0.0D0*vy,vz,tmr,1) ! <|vz|^2>

      CALL energy_arbdir(vx,0.0D0*vy,0.0D0*vz,NNx,0.0D0,NNz,vxpe,vxpa) 
      CALL energy_arbdir(0.0D0*vx,vy,0.0D0*vz,NNx,0.0D0,NNz,vype,vypa) 
      CALL energy_arbdir(0.0D0*vx,0.0D0*vy,vz,NNx,0.0D0,NNz,vzpe,vzpa) 

      CALL energy_arbdir(vx,0.0D0*vy,0.0D0*vz,omegax,omegay,omegaz,vxpeo,vxpao)
      CALL energy_arbdir(0.0D0*vx,vy,0.0D0*vz,omegax,omegay,omegaz,vypeo,vypao)
      CALL energy_arbdir(0.0D0*vx,0.0D0*vy,vz,omegax,omegay,omegaz,vzpeo,vzpao)

!
      IF (myrank.eq.0) THEN
         OPEN(1,file='components.txt',position='append')            !vx  vy  vz
           WRITE(1,FMT='(E25.18,E25.18,E25.18,E25.18,E25.18,E25.18, &
                E25.18,E25.18,E25.18,E25.18,E25.18,E25.18,E25.18,E25.18 &
                ,E25.18,E25.18)') dump,tmp,tmq,tmr,&
                vxpe,vxpa,vype,vypa,vzpe,vzpa,vxpeo,vxpao,vypeo,vypao,vzpeo,vzpao
         CLOSE(1)                 
      ENDIF


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


            IF (myrank.eq.0) THEN
               OPEN(1,file='uf.txt',position='append')
               WRITE(1,FMT='(E25.18,E25.18)') dump,tmr
               CLOSE(1)                          ! u^2
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,FMT='(E13.6,E13.6)') dump,rmp
               CLOSE(1)
            ENDIF
