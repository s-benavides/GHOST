! Global quantities computed in HD runs

            CALL hdcheck(vx,vy,vz,fx,fy,fz,hek,hok,dump,dt,1,0)
            CALL maxabs(vx,vy,vz,rmp,0)
            IF (myrank.eq.0) THEN
               OPEN(1,file='maximum.txt',position='append')
               WRITE(1,FMT='(E13.6,E13.6)') dump,rmp
               CLOSE(1)
            ENDIF
