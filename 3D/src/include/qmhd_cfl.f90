        CALL maxabs(vx,vy,vz,rmp,0)  ! max vort = k U
        tmp = sqrt((omegax)**2+(omegay)**2+(omegaz)**2)
        kcut = real(nx,kind=GP)/3.0_GP    !1/dx        
        tmq = kcut**2*NN ! Lorentz force term

        dt = cfl/max(rmp,tmp,nu*kcut**(2*hek),hnu,tmq)
        CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


