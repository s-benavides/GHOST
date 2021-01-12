import numpy as np
import matplotlib.pyplot as plt
import string

def cdiff(x):
    return 0.5*np.array([x[i+1]-x[i-1] for i in range(1,len(x)-1)])

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Path to the data
runname = raw_input("Folder name: ")
path = '../'+runname+'/run/'

# Reads balance.txt
t,ene,denk,henk,inj,jenk= np.loadtxt(path+'balance.txt',unpack=True)
params = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=116,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

rand = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=142,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[5]
print('rand = ',rand)


nu =  float(params[4])
hnu =  float(params[5])
print("nu = ",nu,"hnu = ",hnu)

# Energy balance
disstot = nu*denk[1:-1]+hnu*henk[1:-1]+jenk[1:-1]
if rand==0:
	print("NONRANDOM FORCING")
	injtot = inj[1:-1]
elif rand==1:
	print("RANDOM FORCING")
	injtot = inj[1:-1]/2.  # BECAUSE RANDOM FORCING
dE = cdiff(ene)/2.
dt = cdiff(t)
deriv = dE/dt
slope = injtot- disstot

differ = deriv-slope

print("mean difference",np.mean(differ),"std",np.std(differ))

t = t[1:-1]

# Plot individually:
plt.figure(1)
plt.clf()
plt.plot(t,slope,'.g',label="Inj-Diss")
plt.plot(t,deriv,'.k',label="dE/dt")
plt.legend()
plt.xlabel('Time')

plt.figure(2)
plt.clf()
plt.semilogy(t,abs(differ),ls='',marker='.',label="dis")
plt.legend()
plt.title("Difference")
plt.xlabel('Time')

plt.show()
