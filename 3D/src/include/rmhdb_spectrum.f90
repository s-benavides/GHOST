! Spectra computed in MHD runs with B_0

            CALL spectrum(vx,vy,vz,ext,1,1,odir)
            CALL specpara(vx,vy,vz,ext,1,1,odir)
            CALL specperp(vx,vy,vz,ext,1,1,odir)
            CALL spectrum(ax,ay,az,ext,0,1,odir)
            CALL specpara(ax,ay,az,ext,0,1,odir)
            CALL specperp(ax,ay,az,ext,0,1,odir)
! Uncomment to write 2D axisymmetric spectra
           CALL spec2D(vx,vy,vz,ext,odir,1,0)
           CALL spec2D_yavg(vx,vy,vz,ext,odir,1)
! Calculating spectrum for magnetic energy, too.
           CALL rotor3(ay,az,c1,1) ! computing bx
           CALL rotor3(ax,az,c2,2)
           CALL rotor3(ax,ay,c3,3)
           CALL spec2D_yavg(c1,c2,c3,ext,odir,0)
!           CALL spec2D(ax,ay,az,ext,odir,0,0)
! Uncomment to save data for spatio-temporal spectra
!           CALL write_fourier(vx,'vx',ext,odir)
!           CALL write_fourier(vy,'vy',ext,odir)
!           CALL write_fourier(vz,'vz',ext,odir)
!           CALL rotor3(ay,az,C1,1)
!           CALL rotor3(ax,az,C2,2)
!           CALL rotor3(ax,ay,C3,3)
!           CALL write_fourier(C1,'bx',ext,odir)
!           CALL write_fourier(C2,'by',ext,odir)
!           CALL write_fourier(C3,'bz',ext,odir)
