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

for i in range(num):
        folder = raw_input("Folder name: ")
        runnames.append(folder)

for jj,runname in enumerate(runnames):
  	path = '../'+runname+'/run/'

	rand = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=136,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[5]
	print('rand = ',rand)

	# Reads balance.txt
	if ('O0' in runname):
		t,enk,enm,denk,denm,henk,henm,jxb = np.loadtxt(path+'balance.txt',unpack=True)
	        t2,injk,injm = np.loadtxt(path+'injection.txt',unpack=True)
		t3,ufk,ufm = np.loadtxt(path+'uf.txt',unpack=True)
	        params_nu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=110,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
	        params_mu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=95,skip_header=61,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

		if len(t)!=len(t2):
			tmin = np.min([len(t),len(t2)])
		else:
			tmin=len(t)

		etot = enk+enm	
		
		if rand==1:
	        	injtot = injk/2. + injm/2. # RANDOM FORCING
		elif rand==0:
			injtot = injk+injm	

	        # nu,kf
	        nu = float(params_nu[4])
	        hnu = float(params_nu[5])
	        hek = float(params_nu[6])
	        hok = float(params_nu[7])
	        kdn = float(params_nu[2])
	        kup = float(params_nu[3])
	        kf = (kdn+kup)/2.
#	        print('nu',nu)
#	        print('hek',hek)
#	        print('kdn',kdn)
#	        print('kup',kup)
	
	        # mu,kf
	        mu = float(params_mu[4])
	        hmu = float(params_mu[5])
	        hem = float(params_mu[6])
	        hom = float(params_mu[7])
	        mkdn = float(params_mu[2])
	        mkup = float(params_mu[3])
	        mkf = (mkdn+mkup)/2.
#	        print('mu',mu)
#	        print('hem',hem)
#	        print('mkdn',mkdn)
#	        print('mkup',mkup)
	else:
		t, enk,denk,henk,injtot = np.loadtxt(path+'balance.txt',unpack=True)
		t2, ufk = np.loadtxt(path+'uf.txt',unpack=True)
		params = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=110,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
		
		if len(t)!=len(t2):
                        tmin = np.min([len(t),len(t2)])
		else:
			tmin = len(t)
		# nu,kf
        	nu = float(params[4])
        	hnu = float(params[5])
		hek = float(params[6])
		hok = float(params[7])
        	kdn = float(params[2])
        	kup = float(params[3])
		kf = (kdn+kup)/2.
#		print('nu',nu)
#		print('hek',hek)
#		print('kdn',kdn)
#		print('kup',kup)
#		print('kf',kf)

                if rand==1:
                        injtot = injtot/2. # RANDOM FORCING

		etot = enk

	# Plots
	if (len(runnames)==1):
		plt.figure(1)
		plt.title(runname)
		plt.plot(t,etot,'-k',label='E')
		
		plt.figure(4)
		plt.title(runname)
		plt.plot(t[:tmin],hnu*henk[:tmin]/injtot[:tmin],'-k',label = 'henk')

		if ('O0' in runname):
			plt.figure(1)
			plt.plot(t,enm,'--b',label='ME')
			plt.plot(t,enk,'--g',label='KE')
		
			plt.figure(4)
			plt.plot(t[:tmin],hmu*henm[:tmin]/injtot[:tmin],'-b',label = 'henm')
		
	else:
		plt.figure(1)
		plt.plot(t,etot,label=runname)
	
		plt.figure(4)
		plt.plot(t[:tmin],hnu*henk[:tmin]/injtot[:tmin],label = runname)
	
		if ('O0' in runname):
			plt.figure(2)
			plt.plot(t,enk,label=runname)
			plt.figure(3)
			plt.plot(t,enm,label=runname)
			plt.figure(4)
			plt.plot(t[:tmin],hmu*henm[:tmin]/injtot[:tmin],'--',label = runname)
			mufm = np.mean(ufm)
			print('run: %s, mean ufb %.3e' % (runname,mufm))
	mufk = np.mean(ufk)
	print('run: %s, mean ufk: %.3e' % (runname,mufk))
			
	minjk=np.mean(injtot)
	print('run: %s, mean inj: %f4' % (runname,minjk))
	
	Re_rms=np.sqrt(np.mean(ufk))/(nu*(kf)**(2*hek-1))
	print('run: %s, Re_rms: %f4' % (runname,Re_rms))
	
plt.figure(1)
plt.xlabel("Time")
plt.ylabel(r"$E_{tot}$")
plt.legend()
plt.tight_layout()

if ((len(runnames)>1) and ('O0' in runname)):
        plt.figure(2)
        plt.xlabel("Time")
        plt.title('KE')
        plt.legend()
	plt.tight_layout()

        plt.figure(3)
        plt.title("ME")
        plt.xlabel("Time")
        plt.legend()
	plt.tight_layout()

plt.figure(4)
plt.xlabel("Time")
plt.ylabel("Hypovisc")
plt.ylim(0,1.0)
plt.legend()
plt.tight_layout()

plt.show()
# Saves plot to an EPS file
#plt.savefig('figure.eps', format='eps', dpi=1000)
