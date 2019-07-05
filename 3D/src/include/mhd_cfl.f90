        CALL maxabs(vx,vy,vz,rmp,0)  ! max vort = k U
        CALL maxabs(ax,ay,az,rmq,1)  ! max j = k B
        kcut = real(nx,kind=GP)/3.0_GP    !1/dx        

        dt = cfl/max(rmp,rmq,nu*kcut**(2*hek),&
                   mu*kcut**(2*hem),hnu,hmu)
        CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


