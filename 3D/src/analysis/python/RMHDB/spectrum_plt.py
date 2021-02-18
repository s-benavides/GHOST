import numpy as np
import glob as glob
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import string

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

cmap = cm.get_cmap('viridis')

# Plots energy spectra as a function of time
# Execute in ipython with '%run plot_spectrum.py'

# Path to the data
numspec = input("how many runs to compare?")

pltname= raw_input("name of plot outputs? ")

runnames = []
for i in range(numspec):
	folder = raw_input("Folder name: ")
	runnames.append(folder)

for jj,runname in enumerate(runnames):
  path = '../'+runname+'/outs/'
  params_nu = np.genfromtxt('../'+runname+'/run/parameter.inp',comments='!',skip_footer=110,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
  # nu,kf
  nu = float(params_nu[4])
  hnu = float(params_nu[5])

  # Reads and plots all spectra in the directory.
  # We only plot one every five spectra, starting
  # from the second.
  kspec = sorted(glob.glob(path+'kspectrum.*.txt'))
  if ('RHD' in runname):
	kspecperp = sorted(glob.glob(path+'kspecperp.*.txt'))
  	kspecpara = sorted(glob.glob(path+'kspecpara.*.txt'))
  if (('MHDB' in runname) or ('RMHDB' in runname)):
        mspec = sorted(glob.glob(path+'mspectrum.*.txt'))
        kspecperp = sorted(glob.glob(path+'kspecperp.*.txt'))
        kspecpara = sorted(glob.glob(path+'kspecpara.*.txt'))
        mspecperp = sorted(glob.glob(path+'mspecperp.*.txt'))
        mspecpara = sorted(glob.glob(path+'mspecpara.*.txt'))
  nfiles = np.size(kspec)
  if len(runnames)==1:
    for i in range(1, nfiles, max(nfiles/6,1)): 
        plt.figure(1)
    	plt.title(runname)
        ene = np.loadtxt(kspec[i])
        plt.loglog(ene[:,0],ene[:,1],color=cmap(i/float(nfiles)),label=str(i),lw=2)
#        plt.loglog(ene[:50,0],10*ene[:50,0]**(-5./3.),'--k',lw=1)

	if ('RHD' in runname):
        	plt.figure(2)
    		plt.title(runname)
        	enperp = np.loadtxt(kspecperp[i])
        	plt.loglog(enperp[:,0],enperp[:,1],color=cmap(i/float(nfiles)),label=str(i),ls='-',lw=2)
#        	plt.loglog(enperp[:10,0],10*enperp[10:80,0]**(-5./3.),'--k',lw=1)

	if (('MHDB' in runname) or ('RMHDB' in runname)):
		plt.figure(3)
    		plt.title(runname)
        	ene = np.loadtxt(mspec[i])
        	plt.loglog(ene[:,0],ene[:,1],color=cmap(i/float(nfiles)),label=str(i),lw=2)
#        	plt.loglog(ene[10:80,0],10*ene[10:80,0]**(-5./3.),'--k',lw=1)
		plt.figure(4)
    		plt.title(runname)
                ene = np.loadtxt(mspecperp[i])
                plt.loglog(ene[:,0],ene[:,1],color=cmap(i/float(nfiles)),label=str(i),lw=2)
#                plt.loglog(ene[10:80,0],10*ene[10:80,0]**(-5./3.),'--k',lw=1)

  else:
    i=-1
    plt.figure(1)
    ene = np.loadtxt(kspec[i])
    if (runname in ['RHD_hyper8']):    
	plt.loglog(ene[:,0],ene[:,1],color = (jj/float(len(runnames)-1),0,0,1),label='Omegaz = 10.0, hnu = %.3e' % hnu,lw=2)
    elif (runname in ['RHD_hyper3','RHD_hyper4','RHD_hyper5']):
    	plt.loglog(ene[:,0],ene[:,1],color = (jj/float(len(runnames)-1),0,0,1),label='Omegaz = 100.0, hnu = %.3e' % hnu,lw=2)
    else:
    	plt.loglog(ene[:,0],ene[:,1],color = (jj/float(len(runnames)-1),0,0,1),label=runname,lw=2)
    #plt.loglog(ene[:10,0],10*ene[:10,0]**(-5./3.),'--k',lw=1)

    if ('RHD' in runname):
	    plt.figure(2)
	    enperp = np.loadtxt(kspecperp[i])
	    enpara = np.loadtxt(kspecpara[i])
	    plt.loglog(enperp[:,0],enperp[:,1],label=runname+r"$\perp$",lw=2,ls='-')
	    plt.loglog(enpara[:,0],enpara[:,1],label=runname+r"$\parallel$",lw=2,ls='--')
#	    plt.loglog(enperp[:10,0],10*enperp[:10,0]**(-5./3.),'--k',lw=1)
	
    if (('MHDB' in runname) or ('RMHDB' in runname)):
	    plt.figure(3)	
            ene = np.loadtxt(mspec[i])
            plt.loglog(ene[:,0],ene[:,1],label=runname,lw=2)
#	    plt.loglog(ene[10:60,0],10*ene[10:60,0]**(-5./3.),'--k',lw=1)
            
	    plt.figure(4)
            ene = np.loadtxt(mspecperp[i])
	    plt.loglog(ene[:,0],ene[:,1],label=runname,lw=2)
#	    plt.loglog(ene[10:60,0],10*ene[10:60,0]**(-5./3.),'--k',lw=1)

plt.figure(1)
plt.xlabel('k')
plt.ylabel(r'$E_{u}(k)$')
plt.ylim([1e-7,10])
plt.legend()
plt.tight_layout()
plt.savefig('./figs/KE_spec_'+pltname+'.png')

if ('RHD' in runname):
	plt.figure(2)
	#plt.xlim(1,nmax/3)
	plt.xlabel('k')
	plt.ylabel(r'$E_\perp(k)$')
	plt.ylim([1e-7,10])
	plt.legend()
	plt.tight_layout()
	plt.savefig('./figs/KEperp_spec_'+pltname+'.png')

if (('MHDB' in runname) or ('RMHDB' in runname)):
        plt.figure(3)
	plt.xlabel('k')
	plt.ylabel(r'$E_b(k)$')
	plt.ylim([1e-7,10])
	plt.legend()
	plt.tight_layout()
	plt.savefig('./figs/ME_spec_'+pltname+'.png')

	plt.figure(4)
	plt.xlabel('k')
	plt.ylabel(r"$E_{b,\perp}(k)$")
	plt.ylim([1e-7,10])
	plt.legend()
	plt.tight_layout()
	plt.savefig('./figs/MEperp_spec_'+pltname+'.png')
plt.show()	
