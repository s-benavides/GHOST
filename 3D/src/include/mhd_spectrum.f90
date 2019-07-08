! Spectra computed in MHD runs

            CALL spectrum(vx,vy,vz,ext,1,1,odir)
            CALL spectrum(ax,ay,az,ext,0,1,odir)
            CALL crosspec(vx,vy,vz,ax,ay,az,ext,odir)
