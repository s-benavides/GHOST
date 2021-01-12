import numpy as np
import matplotlib.pyplot as plt
import string

# Plots energy as a function of time
# Assumes balance.txt is the output of an HD run
# Execute in ipython with '%run plot_energy.py'

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Path to the data
num = input("how many runs to compare?")
runnames = []

pltname = raw_input('name for plot? ')

for i in range(num):
        folder = raw_input("Folder name: ")
        runnames.append(folder)

for jj,runname in enumerate(runnames):
  	path = '../'+runname+'/run/'

	# Reads balance.txt
	if (('MHD' in runname) and ('test2' not in runname)):
		t,enk,enm,denk,denm,henk,henm,jxb = np.loadtxt(path+'balance.txt',unpack=True)
	        t2,injk,injm = np.loadtxt(path+'injection.txt',unpack=True)
		t3,ufk,ufm = np.loadtxt(path+'uf.txt',unpack=True)
	        params_nu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=110,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
	        params_mu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=95,skip_header=61,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
	
		etot = enk+enm	

	        injtot = injk/2. + injm/2. # RANDOM FORCING
	
	        # nu,kf
	        nu = float(params_nu[4])
	        hek = float(params_nu[6])
	        kdn = float(params_nu[2])
	        kup = float(params_nu[3])
	        kf = (kdn+kup)/2.
#	        print('nu',nu)
#	        print('hek',hek)
#	        print('kdn',kdn)
#	        print('kup',kup)
	
	        # mu,kf
	        mu = float(params_mu[4])
	        hem = float(params_mu[6])
	        mkdn = float(params_mu[2])
	        mkup = float(params_mu[3])
	        mkf = (mkdn+mkup)/2.
#	        print('mu',mu)
#	        print('hem',hem)
#	        print('mkdn',mkdn)
#	        print('mkup',mkup)
	elif ('test2' in runname):
		t,etot,denk,denm = np.loadtxt(path+'balance.txt',unpack=True)
		t,enk,enm = np.loadtxt(path+'energy.txt',unpack=True)
                t2,injk,injm = np.loadtxt(path+'injection.txt',unpack=True)
                params_nu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=105,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
                params_mu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=89,skip_header=57,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
		print(params_nu)
		print(params_mu)
                injtot = injk/2. + injm/2. # RANDOM FORCING

                # nu,kf
                nu = float(params_nu[4])
                hek = 1.
                kdn = float(params_nu[2])
                kup = float(params_nu[3])
                kf = (kdn+kup)/2.
                print('nu',nu)
                print('hek',hek)
                print('kdn',kdn)
                print('kup',kup)

                # mu,kf
                mu = float(params_mu[4])
                hem = 1.
                mkdn = float(params_mu[2])
                mkup = float(params_mu[3])
                mkf = (mkdn+mkup)/2.
                print('mu',mu)
                print('hem',hem)
                print('mkdn',mkdn)
                print('mkup',mkup)
	else:
		t, enk,denk,henk,injtot = np.loadtxt(path+'balance.txt',unpack=True)
		t2, ufk = np.loadtxt(path+'uf.txt',unpack=True)
		params = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=110,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
		# nu,kf
        	nu = float(params[4])
		hek = float(params[6])
        	kdn = float(params[2])
        	kup = float(params[3])
		kf = (kdn+kup)/2.
#		print('nu',nu)
#		print('hek',hek)
#		print('kdn',kdn)
#		print('kup',kup)
#		print('kf',kf)

		injtot = injtot/2. #Random forcing

		etot = enk

	# Plots
	if (len(runnames)==1):
		plt.figure(1)
		plt.title(runname)
		plt.plot(t,etot,'-k',label='E')
		
		if ('MHD' in runname):
			plt.plot(t,enm,'--b',label='ME')
			plt.plot(t,enk,'--g',label='KE')

		
	else:
		plt.figure(1)
		plt.plot(t,etot,label=runname)
		
		if ('MHD' in runname):
			plt.figure(2)
			plt.plot(t,enk,label=runname)
			plt.figure(3)
			plt.plot(t,enm,label=runname)

	minjk=np.mean(injtot)
	print('run: %s, mean inj: %f4' % (runname,minjk))

	if ('test2' in runname):
		menk = np.mean(enk)
		Re_rms=np.sqrt(menk)/(nu*(kf)**(2*hek-1))
	else:
		Re_rms=np.sqrt(np.mean(ufk))/(nu*(kf)**(2*hek-1))
	print('run: %s, Re_rms: %f4' % (runname,Re_rms))
	
plt.figure(1)
plt.xlabel("Time")
plt.ylabel(r"$E_{tot}$")
plt.legend()
plt.tight_layout()
plt.savefig('./figs/E_'+pltname+'.png')

if ((len(runnames)>1) and ('MHD' in runname)):
        plt.figure(2)
        plt.xlabel("Time")
        plt.title('KE')
        plt.legend()
	plt.tight_layout()
	plt.savefig('./figs/KE_'+pltname+'.png')

        plt.figure(3)
        plt.title("ME")
        plt.xlabel("Time")
        plt.legend()
	plt.tight_layout()
	plt.savefig('./figs/ME_'+pltname+'.png')

plt.show()
# Saves plot to an EPS file
#plt.savefig('figure.eps', format='eps', dpi=1000)
