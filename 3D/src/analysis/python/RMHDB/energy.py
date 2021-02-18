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

	rand = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=136,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[5]
	print('rand = ',rand)

	# Reads balance.txt
	t,enk,enm,denk,denm,henk,henm,jxb = np.loadtxt(path+'balance.txt',unpack=True)
	t2,injk,injm = np.loadtxt(path+'injection.txt',unpack=True)
	t3,ufk,ufm = np.loadtxt(path+'uf.txt',unpack=True)
	t4,eperpb,eparab,eperpo,eparao = np.loadtxt(path+'eperp_epara.txt',unpack=True)
	if os.path.exists(path+'ebperp_ebpara.txt'):
		t5,ebperpb,ebparab,ebperpo,ebparao = np.loadtxt(path+'ebperp_ebpara.txt',unpack=True)

	params_nu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=110,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
	params_mu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=95,skip_header=61,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

	tmin = np.min([len(t),len(t2),len(t3),len(t4)])

	etot = enk+enm	
		
	if rand==1:
        	injtot = injk/2. + injm/2. # RANDOM FORCING
	elif rand==0:
		injtot = injk+injm	

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

        # mu,kf
        mu = float(params_mu[4])
        hmu = float(params_mu[5])
        hem = float(params_mu[6])
        hom = float(params_mu[7])
        mkdn = float(params_mu[2])
        mkup = float(params_mu[3])
        mkf = (mkdn+mkup)/2.

	# Plots
	if (len(runnames)==1):
		plt.figure(1)
		plt.title(runname)
		plt.plot(t,etot,'-k',label='E')
		plt.plot(t,enm,'--b',label='ME')
		plt.plot(t,enk,'--g',label='KE')
		
		plt.figure(4)
		plt.title(runname)
		plt.plot(t[:tmin],hnu*henk[:tmin]/injtot,'-k',label = 'henk/<injtot>')
		plt.plot(t[:tmin],hmu*henm[:tmin]/injtot,'-b',label = 'henm/<injtot>')


		plt.figure(5)
		plt.title(runname)
		plt.plot(t4,eperpb/enk,'-b',label = 'eperpb/enk')
		plt.plot(t4,eperpo/enk,'--r',label = 'eperpo/enk')

                plt.figure(6)
                plt.title(runname)
                plt.plot(t4,eparab/enk,'-b',label = 'eparab/enk')
                plt.plot(t4,eparao/enk,'--r',label = 'eparao/enk')

		if os.path.exists(path+'ebperp_ebpara.txt'):
			plt.figure(7)
                	plt.title(runname)
                	plt.plot(t5,ebperpb/enm[-len(ebperpb):],'-b',label = 'ebperpb/enb')
                	plt.plot(t5,ebperpo/enm[-len(ebperpb):],'--r',label = 'ebperpo/enb')
	
                	plt.figure(8)
                	plt.title(runname)
                	plt.plot(t5,ebparab/enm[-len(ebperpb):],'-b',label = 'ebparab/enb')
                	plt.plot(t5,ebparao/enm[-len(ebperpb):],'--r',label = 'ebparao/enb')		


	else:
		plt.figure(1)
		plt.plot(t,etot,label=runname)
	
		plt.figure(4)
		plt.plot(t[:tmin],hnu*henk[:tmin]/injtot,label = runname)
	
		plt.figure(2)
		plt.plot(t,enk/etot,label=runname)
	
		plt.figure(3)
		plt.plot(t,enm/etot,label=runname)
	
		plt.figure(4)
		plt.plot(t[:tmin],hmu*henm[:tmin]/injtot,'--',label = runname)
		

                plt.figure(5)
                plt.plot(t4[:tmin],eperpb[:tmin]/enk[:tmin],label = runname)

		plt.figure(7)
                plt.plot(t4[:tmin],eperpo[:tmin]/enk[:tmin],label = runname)

                plt.figure(6)
                plt.plot(t4[:tmin],eparab[:tmin]/enk[:tmin],label = runname)

		plt.figure(8)
                plt.plot(t4[:tmin],eparao[:tmin]/enk[:tmin],label = runname)

		if os.path.exists(path+'ebperp_ebpara.txt'):
			pltb=1
	                plt.figure(9)
	                plt.plot(t5,ebperpb/enm[-len(ebperpb):],label = runname)
	
	                plt.figure(10)
	                plt.plot(t5,ebperpo/enm[-len(ebperpb):],label = runname)
	
	                plt.figure(11)
	                plt.plot(t5,ebparab/enm[-len(ebperpb):],label = runname)
		
        	        plt.figure(12)
        	        plt.plot(t5,ebparao/enm[-len(ebperpb):],label = runname)


	mufm = np.mean(ufm)
	print('run: %s, mean ufb %.3e' % (runname,mufm))
	mufk = np.mean(ufk)
	print('run: %s, mean ufk: %.3e' % (runname,mufk))
			
	print('run: %s, mean injtot: %f4' % (runname,injtot))
	
	Re_rms=np.sqrt(np.mean(ufk))/(nu*(kf)**(2*hek-1))
	print('run: %s, Re_rms: %f4' % (runname,Re_rms))

	Rem_rms=np.sqrt(np.mean(ufk))/(mu*(kf)**(2*hem-1))
	print('run: %s, Rem_rms: %f4' % (runname,Rem_rms))
	
plt.figure(1)
plt.xlabel("Time")
plt.ylabel(r"$E_{tot}$")
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/Etot_'+svname+'.png')


if (len(runnames)>1):
        plt.figure(2)
        plt.xlabel("Time")
        plt.title(r'$KE/E_{tot}$')
        plt.legend()
	plt.tight_layout()
	if svfig in ['Y']: plt.savefig('./figs/KE_'+svname+'.png')

        plt.figure(3)
        plt.title(r"$ME/E_{tot}$")
        plt.xlabel("Time")
        plt.legend()
	plt.tight_layout()
	if svfig in ['Y']: plt.savefig('./figs/ME_'+svname+'.png')

	plt.figure(7)
        plt.xlabel("Time")
        plt.title(r'$KE^{\Omega}_{\perp}/KE$')
        plt.legend()
        plt.tight_layout()
	if svfig in ['Y']: plt.savefig('./figs/eperpo_'+svname+'.png')

        plt.figure(8)
        plt.xlabel("Time")
        plt.title(r'$KE^{\Omega}_{\parallel}/KE$')
        plt.legend()
        plt.tight_layout()
	if svfig in ['Y']: plt.savefig('./figs/eparapo_'+svname+'.png')
	if pltb==1:
		plt.figure(9)
        	plt.xlabel("Time")
        	plt.title(r'$ME^{B}_{\perp}/ME$')
        	plt.legend()
        	plt.tight_layout()
		if svfig in ['Y']: plt.savefig('./figs/ebperpb_'+svname+'.png')
	
        	plt.figure(10)
        	plt.xlabel("Time")
        	plt.title(r'$ME^{\Omega}_{\perp}/ME$')
        	plt.legend()
        	plt.tight_layout()
		if svfig in ['Y']: plt.savefig('./figs/ebperpo_'+svname+'.png')

                plt.figure(11)
                plt.xlabel("Time")
                plt.title(r'$ME^{B}_{\parallel}/ME$')
                plt.legend()
                plt.tight_layout()
		if svfig in ['Y']: plt.savefig('./figs/ebparapb_'+svname+'.png')

                plt.figure(12)
                plt.xlabel("Time")
                plt.title(r'$ME^{\Omega}_{\parallel}/ME$')
                plt.legend()
                plt.tight_layout()
		if svfig in ['Y']: plt.savefig('./figs/ebparapo_'+svname+'.png')

		

elif os.path.exists(path+'ebperp_ebpara.txt'):
        plt.figure(7)
        plt.xlabel("Time")
        plt.title(r'$ME_{\perp}/ME$')
        plt.legend()
        plt.tight_layout()
	if svfig in ['Y']: plt.savefig('./figs/ebperp_'+svname+'.png')

        plt.figure(8)
        plt.xlabel("Time")
        plt.title(r'$ME_{\parallel}/ME$')
        plt.legend()
        plt.tight_layout()
	if svfig in ['Y']: plt.savefig('./figs/ebpara_'+svname+'.png')
	

plt.figure(4)
plt.xlabel("Time")
plt.ylabel("Hypovisc")
#plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()
if svfig in ['Y']: plt.savefig('./figs/hyper_'+svname+'.png')

plt.figure(5)
plt.xlabel("Time")
plt.ylabel(r"$KE_\perp$")
plt.legend()
plt.tight_layout()
if (len(runnames)>1):
	plt.ylabel(r"$KE^{B}_{\perp}/KE$")
	if svfig in ['Y']: plt.savefig('./figs/eperpb_'+svname+'.png')
else:
	if svfig in ['Y']: plt.savefig('./figs/eperp_'+svname+'.png')	

plt.figure(6)
plt.xlabel("Time")
plt.ylabel(r"$KE_\parallel$")
plt.legend()
plt.tight_layout()
if (len(runnames)>1):
	plt.ylabel(r"$KE^{B}_{\parallel}/KE$")
	if svfig in ['Y']: plt.savefig('./figs/eparab_'+svname+'.png')
else:
	if svfig in ['Y']: plt.savefig('./figs/epara_'+svname+'.png')
plt.show()
