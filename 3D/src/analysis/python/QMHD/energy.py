import numpy as np
import matplotlib.pyplot as plt
import string
import os.path

# Plots energy as a function of time
# Assumes balance.txt is the output of an HD run
# Execute in ipython with '%run plot_energy.py'

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Path to the data
num = input("how many runs to compare?")
runnames = []

# For saving:
svfig = raw_input("Save figures? (Y or N): ")
if svfig in ['Y']:
    svname = raw_input("Name for save: ")


for i in range(num):
        folder = raw_input("Folder name: ")
        runnames.append(folder)

for jj,runname in enumerate(runnames):
    path = '../'+runname+'/run/'

    rand = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=141,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[5]
    print('rand = ',rand)

    # Reads balance.txt
    t,enk,denk,henk,injk,jenk = np.loadtxt(path+'balance.txt',unpack=True)
    t2,ufk = np.loadtxt(path+'uf.txt',unpack=True)
#    if os.path.exists(path+'ebperp_ebpara.txt'):
#    t5,ebperpb,ebparab,ebperpo,ebparao = np.loadtxt(path+'ebperp_ebpara.txt',unpack=True)
    
    params_nu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=116,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

    tmin = np.min([len(t),len(t2)])

    etot = enk
    
    if rand==1:
        injtot = injk/2. # RANDOM FORCING
    elif rand==0:
        injtot = injk

    # Averaging injection
    injtot = np.mean(injtot)

    # nu,kf
    nu = float(params_nu[4])
    hnu = float(params_nu[5])
    hek = float(params_nu[6])
    hok = float(params_nu[7])
    kdn = float(params_nu[2])
    kup = float(params_nu[3])
    kf = (kdn+kup)/2.

    # Plots
    if (len(runnames)==1):
        plt.figure(1)
        plt.title(runname)
        plt.plot(t,etot,'-k',label='E')
        
        plt.figure(2)
        plt.title(runname)
        plt.plot(t,hnu*henk/injtot,'-k',label = 'henk/<injtot>')
        
        plt.figure(3)
        plt.title(runname)
        plt.plot(t,jenk/injtot,'-k',label = 'henk/<injtot>')
        
        plt.figure(4)
        plt.title(runname)
        plt.plot(t,denk/injtot,'-k',label='denk/<injtot>')

#		plt.figure(5)
#		plt.title(runname)
#		plt.plot(t4,eperpb/enk,'-b',label = 'eperpb/enk')
#		plt.plot(t4,eperpo/enk,'--r',label = 'eperpo/enk')
#
#                plt.figure(6)
#                plt.title(runname)
#                plt.plot(t4,eparab/enk,'-b',label = 'eparab/enk')
#                plt.plot(t4,eparao/enk,'--r',label = 'eparao/enk')

#		if os.path.exists(path+'ebperp_ebpara.txt'):
#			plt.figure(7)
#                	plt.title(runname)
#                	plt.plot(t5,ebperpb/enm[-len(ebperpb):],'-b',label = 'ebperpb/enb')
#                	plt.plot(t5,ebperpo/enm[-len(ebperpb):],'--r',label = 'ebperpo/enb')
#	
#                	plt.figure(8)
#                	plt.title(runname)
#                	plt.plot(t5,ebparab/enm[-len(ebperpb):],'-b',label = 'ebparab/enb')
#                	plt.plot(t5,ebparao/enm[-len(ebperpb):],'--r',label = 'ebparao/enb')		


    else:
        plt.figure(1)
        plt.plot(t,etot,label=runname)
        
        plt.figure(2)
        plt.plot(t,hnu*henk/injtot,label = runname)
    
        plt.figure(3)
        plt.plot(t,jenk/injtot,label=runname)
    
        plt.figure(4)
        plt.plot(t,denk/injtot,label=runname)

    mufk = np.mean(ufk)
    print('run: %s, mean ufk: %.3e' % (runname,mufk))
            
    print('run: %s, mean injtot: %f4' % (runname,injtot))

    Re_rms=np.sqrt(np.mean(ufk))/(nu*(kf)**(2*hek-1))
    print('run: %s, Re_rms: %f4' % (runname,Re_rms))

plt.figure(1)
plt.xlabel("Time")
plt.ylabel(r"$E_{tot}$")
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/Etot_'+svname+'.png')


plt.figure(2)
plt.xlabel("Time")
plt.ylabel("Hypovisc")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/hyper_'+svname+'.png')

plt.figure(3)
plt.xlabel("Time")
plt.ylabel("Joule Diss")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/joule_'+svname+'.png')

plt.figure(4)
plt.xlabel("Time")
plt.ylabel("Viscous Diss")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/visc_'+svname+'.png')

plt.show()
