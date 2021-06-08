import numpy as np

O = input("What do you want O to be?")
O = float(O)
Theta = input("What do you want theta to be (in degrees)?")
Theta = float(Theta)*(np.pi/180)

Oz = O*np.sin(Theta)
Ox = O*np.cos(Theta)

print("(omegax,omegaz) = (%.7f,%.7f)" % (Ox,Oz))
