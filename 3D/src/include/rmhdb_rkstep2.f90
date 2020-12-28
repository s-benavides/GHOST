! Step 2 of Runge-Kutta for the MHD equations with B_0 and rotation
! Computes the nonlinear terms and evolves the equations in dt/o

         CALL rotor3(ay,az,C7,1)
         CALL rotor3(ax,az,C8,2)
         CALL rotor3(ax,ay,C9,3)
         IF (myrank.eq.0) THEN          ! b = b + B_0
            C7(1,1,1) = bx0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            C8(1,1,1) = by0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
            C9(1,1,1) = bz0*real(nx,kind=GP)*real(ny,kind=GP)*real(nz,kind=GP)
         ENDIF
         CALL prodre3(vx,vy,vz,C10,C11,C12) ! calculates curl(v) x v
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend               ! Coriolis force
!$omp parallel do if (iend-ista.lt.nth) private (k)
            DO j = 1,ny
               DO k = 1,nz
                  C10(k,j,i) = C10(k,j,i)-2*omegaz*vy(k,j,i)
                  C11(k,j,i) = C11(k,j,i)+2*omegaz*vx(k,j,i)
               END DO
            END DO
         END DO
         CALL prodre3(C7,C8,C9,C13,C14,C15) ! calculates curl(b+B_0) x(b+ B_0)
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            CALL entrans(C1,C2,C3,C13,C14,C15,ext,2,odir) ! flux of v . (J x (b+B_0))
            CALL entpara(C1,C2,C3,C13,C14,C15,ext,2,odir) 
            CALL entperp(C1,C2,C3,C13,C14,C15,ext,2,odir)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C16,1)  ! computes-(v.grad)v+(b.grad)b-grad(p)
                                                      ! [or -curl(v)xv+curl(b)xb-grad(p)] in Fourier
                                                      ! space, with the pressure chosen to satisfy the
                                                      ! incompressibility condition.
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C17,2)
         CALL nonlin3(C10,C11,C12,C13,C14,C15,C10,3)
         CALL vector3(vx,vy,vz,C7,C8,C9,C11,C12,C13)  ! computes v x (b+B_0)
         CALL gauge3(C11,C12,C13,C7,1)  ! Computes the nonlinear terms in the induction
                                        ! equation for the vector potential, imposing a
                                        ! gauge that satisfies the condition div(A)=0.

         CALL gauge3(C11,C12,C13,C8,2)
         CALL gauge3(C11,C12,C13,C9,3)
         CALL diss(vx,vx,hek,hok,nu,hnu)
         CALL diss(vy,vy,hek,hok,nu,hnu)
         CALL diss(vz,vz,hek,hok,nu,hnu)
         CALL diss(ax,ax,hem,hom,mu,hmu)
         CALL diss(ay,ay,hem,hom,mu,hmu)
         CALL diss(az,az,hem,hom,mu,hmu)
         IF ((trans.eq.1).and.(times.eq.0).and.(bench.eq.0).and.(o.eq.ord)) &
            THEN
            CALL entrans(C1,C2,C3,C16,C17,C10,ext,1,odir) ! flux of v.(-v.grad(v) + Jx(b+B)0) -grad(p))
            CALL entpara(C1,C2,C3,C16,C17,C10,ext,1,odir)
            CALL entperp(C1,C2,C3,C16,C17,C10,ext,1,odir)
            CALL entrans(C4,C5,C6,C7,C8,C9,ext,0,odir)
            CALL entpara(C4,C5,C6,C7,C8,C9,ext,0,odir)
            CALL entperp(C4,C5,C6,C7,C8,C9,ext,0,odir)
         ENDIF

         rmp = 1./real(o,kind=GP)
!$omp parallel do if (iend-ista.ge.nth) private (j,k)
         DO i = ista,iend
!$omp parallel do if (iend-ista.ge.nth) private (k)
         DO j = 1,ny
         DO k = 1,nz

            IF ((kn2(k,j,i).le.kmax).and.(kn2(k,j,i).ge.tiny)) THEN
               vx(k,j,i) = C1(k,j,i)+dt*(vx(k,j,i)+C16(k,j,i) &
              +fx(k,j,i))*rmp
               vy(k,j,i) = C2(k,j,i)+dt*(vy(k,j,i)+C17(k,j,i) &
              +fy(k,j,i))*rmp
               vz(k,j,i) = C3(k,j,i)+dt*(vz(k,j,i)+C10(k,j,i) &
              +fz(k,j,i))*rmp
               ax(k,j,i) = C4(k,j,i)+dt*(ax(k,j,i)+C7(k,j,i)  &
              +mx(k,j,i))*rmp
               ay(k,j,i) = C5(k,j,i)+dt*(ay(k,j,i)+C8(k,j,i)  &
              +my(k,j,i))*rmp
               az(k,j,i) = C6(k,j,i)+dt*(az(k,j,i)+C9(k,j,i)  &
              +mz(k,j,i))*rmp
            ELSE
               vx(k,j,i) = 0.0_GP
               vy(k,j,i) = 0.0_GP
               vz(k,j,i) = 0.0_GP
               ax(k,j,i) = 0.0_GP
               ay(k,j,i) = 0.0_GP
               az(k,j,i) = 0.0_GP
            ENDIF

         END DO
         END DO
         END DO
