! Global quantities computed in MHD runs

            CALL mhdcheck(ps,az,fk,mk,dump,dt)
            CALL maxabs(ps,rmp,ki,kj)
            CALL maxabs(az,rmq,ki,kj)
            IF (myrank.eq.0) THEN
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,'(E13.6,E13.6,E13.6)') dump,rmp,rmq
               CLOSE(1)
            ENDIF
