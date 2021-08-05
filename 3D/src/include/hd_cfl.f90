        CALL maxabs(vx,vy,vz,rmp,0)  ! max vort = k U
        kcut = max(real(nx,kind=GP)/Lx,real(ny,kind=GP)/Ly,real(nz,kind=GP)/Lz)/3.0_GP !1/dx

        dt = cfl/max(rmp,nu*kcut**(2*hek),hnu)
        CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


