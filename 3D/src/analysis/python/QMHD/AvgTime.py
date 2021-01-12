import numpy as np
import matplotlib.pyplot as plt
import string
import bunch_err

# Plots energy as a function of time
# Assumes balance.txt is the output of an HD run
# Execute in ipython with '%run plot_energy.py'

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Path to the data
Otemp = raw_input("What O value? ")
Btemp = raw_input("What Nx value? ")
Ttemp = raw_input("What Nz value? ")

runname = 'O'+Otemp+'Nx'+Btemp+'Nz'+Ttemp

path = '../'+runname+'/run/'

rand = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=142,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[5]
	
sstep = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=142,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[3]

cstep = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=142,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[4]

# Reads balance.txt
t,enk,denk,henk,injk,jenk = np.loadtxt(path+'balance.txt',unpack=True)
t2,ufk = np.loadtxt(path+'uf.txt',unpack=True)
params_nu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=116,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

if len(t)!=len(t2):
        tmin = np.min([len(t),len(t2)])
else:
        tmin=len(t)

if rand==1:
       	injtot = injk/2. # RANDOM FORCING
elif rand==0:
	injtot = injk

injtot = np.mean(injtot)

# nu,kf
nu = float(params_nu[4])
hnu = float(params_nu[5])
hek = float(params_nu[6])
hok = float(params_nu[7])
kdn = float(params_nu[2])
kup = float(params_nu[3])
kf = (kdn+kup)/2.
	
print('rand',rand,'sstep',sstep,'cstep',cstep,'nu',nu,'hnu',hnu,'hek',hek,'kup',kup)

# Plots
plt.figure(1)
plt.title(runname)
plt.plot(enk,'-k',label='E')
	
plt.figure(2)
plt.title(runname)
plt.plot(hnu*henk[:tmin]/injtot,'-k',label = 'henk/<injtot>')
	
plt.figure(3)
plt.plot(jenk[:tmin]/injtot,'-b',label = 'jenk/<injtot>')

mufk = np.mean(ufk)
print('run: %s, mean ufk: %.3e' % (runname,mufk))		
minjk=np.mean(injtot)
print('run: %s, mean inj: %f4' % (runname,minjk))
Re_rms=np.sqrt(np.mean(ufk))/(nu*(kf)**(2*hek-1))
print('run: %s, Re_rms: %f4' % (runname,Re_rms))
	
plt.figure(1)
plt.xlabel("Output #")
plt.ylabel(r"$E_{tot}$")
plt.legend()
plt.tight_layout()

plt.figure(2)
plt.xlabel("Output #")
plt.ylabel("Hypovisc")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()

plt.figure(3)
plt.xlabel("Output #")
plt.ylabel("Joule Diss")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()

plt.show()

start = input('Enter x-axis value you want to start averaging: ')
start_fl = start/float(sstep/cstep)        #starting flux and spectra number

# Bunching algorithm error index:
bunch_err.bunch_err(etot[start:])
err_ind = input('Enter iteration to take error: ')

np.savetxt('rundat/AvgTime'+runname+'.txt',[start,start_fl,err_ind], delimiter ='   ')

