! Global quantities computed in MHD runs
#ifdef CFL_
        ta = time
#else
        ta = t
#endif
            CALL mhdcheck(vx,vy,vz,ax,ay,az,ta,dt,1,1,0)
            CALL cross(vx,vy,vz,fx,fy,fz,eps,1)
            CALL cross(ax,ay,az,mx,my,mz,epm,0)
            CALL maxabs(vx,vy,vz,rmp,0)
            CALL maxabs(ax,ay,az,rmq,1)
            IF (myrank.eq.0) THEN
               OPEN(1,file='injection.txt',position='append')
#ifdef CFL_
               WRITE(1,FMT='(E13.6,E22.14,E22.14)') time,eps,epm
#else
               WRITE(1,FMT='(E13.6,E22.14,E22.14)') (t-1)*dt,eps,epm
#endif
               CLOSE(1)
               OPEN(1,file='maximum.txt',position='append')
#ifdef CFL)
               WRITE(1,FMT='(E13.6,E13.6,E13.6)') time,rmp,rmq
#else
               WRITE(1,FMT='(E13.6,E13.6,E13.6)') (t-1)*dt,rmp,rmq

#endif
               CLOSE(1)
            ENDIF
