import numpy as np

N = input("What do you want N to be?")
N = float(N)
Theta = input("What do you want theta to be (in degrees)?")
Theta = float(Theta)*(np.pi/180)

kf = 9
inj = 1
u = ((4/5)*inj/kf)**(1/3.)

# Calculate magnitude of (Nx,Nz)
Nmag =  np.sqrt(N*u*kf)

Nx = Nmag*np.sin(Theta)
Nz = Nmag*np.cos(Theta)
print("(Nx,Nz) = (%.7f,%.7f)" % (Nx,Nz))
