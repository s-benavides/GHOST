        CALL maxabs(vx,vy,vz,rmp,0)  ! max vort = k U
        tmp = sqrt((omegax)**2+(omegay)**2+(omegaz)**2)
        kcut = max(real(nx,kind=GP)/Lx,real(ny,kind=GP)/Ly,real(nz,kind=GP)/Lz)/3.0_GP    !1/dx        
        tmq = NNx**2+2*NNx*NNz+NNz**2 ! Lorentz force term

        dt = cfl/max(rmp,tmp,nu*kcut**(2*hek),hnu,tmq)
        CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


