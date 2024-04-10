! Global quantities computed in HD runs
            CALL hdcheck(ps,fk,dump,nu,hnu,hek,hok,kdn,kup)
!           CALL maxabs(ps,rmp,ki,kj)
!           IF (myrank.eq.0) THEN
!              OPEN(1,file='maximum.txt',position='append')
!              WRITE(1,'(E13.6,E13.6)') (t-1)*dt,rmp
!              CLOSE(1)
!           ENDIF