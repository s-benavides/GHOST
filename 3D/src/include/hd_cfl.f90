        CALL maxabs(vx,vy,vz,rmp,0)  ! max vort = k U
        kcut = real(nx,kind=GP)/3.0_GP    !1/dx        

        dt = cfl/max(rmp,nu*kcut**(2*hek),hnu)
        CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


