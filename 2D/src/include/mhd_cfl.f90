        CALL maxabs(ps,rmp)  ! max vorticity
        CALL maxabs(az,rmq)  ! max j = k B
        kcut = real(n,kind=GP)/3.0_GP !1/dx

        dt = cfl/max(rmp,rmq,nu*kcut**(2*hek),&
                   mu*kcut**(2*hem),hnu,hmu)
        CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


