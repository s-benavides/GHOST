import numpy as np
import matplotlib.pyplot as plt
import string

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

def cdiff(x):
    return 0.5*np.array([x[i+1]-x[i-1] for i in range(1,len(x)-1)])


# Path to the data
runname = raw_input("Folder name: ")
path = '../'+runname+'/run/'

# Reads balance.txt
t,enk,enm,denk,denm,henk,henm,jxb = np.loadtxt(path+'balance.txt',unpack=True)
t2,injk,injm = np.loadtxt(path+'injection.txt',unpack=True)
params_nu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=110,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
params_mu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=95,skip_header=61,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

rand = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=136,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[5]
print('rand = ',rand)


nu =  float(params_nu[4])
hnu =  float(params_nu[5])
mu =  float(params_mu[4])
hmu =  float(params_mu[5])
print("nu = ",nu,"hnu = ",hnu)
print("mu = ",mu,"hmu = ",hmu)

# Energy balance
disstot = nu*denk[1:-1]+hnu*henk[1:-1]+mu*denm[1:-1]+hmu*henm[1:-1]
injtot = injk[1:-1]+injm[1:-1] # RANDOM FORCING
if rand==0:
	print("NONRANDOM FORCING")
elif rand==1:
	injtot = injtot/2.
	print("RANDOM FORCING")
ene = enk+enm
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
