        CALL maxabs(vx,vy,vz,rmp,0)  ! max vort = k U
        CALL maxabs(ax,ay,az,rmq,1)  ! max j = k B
        tmp = sqrt((omegax)**2+(omegay)**2+(omegaz)**2)
        tmq = sqrt((bx0)**2+(by0)**2+(bz0)**2)

        dt = cfl/max(tmp,rmq,tmp,kcut*tmq,nu*kcut**2,mu*kcut**2)
      CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


