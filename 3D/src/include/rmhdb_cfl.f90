        CALL maxabs(vx,vy,vz,rmp,0)  ! max vort = k U
        CALL maxabs(ax,ay,az,rmq,1)  ! max j = k B
        tmp = sqrt((omegax)**2+(omegay)**2+(omegaz)**2)
        tmq = sqrt((bx0)**2+(by0)**2+(bz0)**2)
        kcut = max(real(nx,kind=GP)/Lx,real(ny,kind=GP)/Ly,real(nz,kind=GP)/Lz)/3.0_GP !1/dx

        dt = cfl/max(rmp,rmq,tmp,kcut*tmq,nu*kcut**(2*hek),&
                   mu*kcut**(2*hem),hnu,hmu)
        CALL MPI_BCAST(dt,1,GC_REAL,0,MPI_COMM_WORLD,ierr)


