import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
import pickle
import glob as glob
import string
from datetime import date

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Get all the AvgTimeO*B*.txt files
avglist = sorted(glob.glob('rundat/*T90*.txt'))

runnames = []
for file in avglist:
  run = file.split('rundat/AvgTime')[1]
  run = run.split('.txt')[0]
  runnames.append(run)

Data = dict([])
for i,run in enumerate(runnames):
        Data_E = dict([])
        print("Working on run %s " % run)
        path = '../'+run+'/outs/'
	
	 # Average start indices
        [start,start_fl] = np.loadtxt('rundat/AvgTime'+run+'.txt')

        start = int(start)
        start_fl = int(start_fl)

	inds = []
	files = sorted(glob.glob(path+'kspectrum.*.txt'))
	for file in files:
		num = file.split(path+'kspectrum.')[1]
  		num = num.split('.txt')[0]
		if (int(num)>=start_fl):
 			inds.append(num)

	numfiles = int(len(inds))
	print('%s files to average' % numfiles)

        # Reads flux files
        if ('B0' in run):  # In case ROTH or HD was run instead of RMHDB
		# Averaging
                for ind in inds:
			kspec = sorted(glob.glob(path+'kspectrum.'+str(ind)+'.txt'))[0]
			mspec = np.nan
	
	        	kspecperp = sorted(glob.glob(path+'kspecperp.'+str(ind)+'.txt'))[0]
			mspecperp = np.nan
	
	        	kspecpara = sorted(glob.glob(path+'kspecpara.'+str(ind)+'.txt'))[0]
			mspecpara = np.nan


                        # Load data:
                        kspec_avg += np.loadtxt(kspec)[:,1]/float(numfiles)

                        kspecperp_avg += np.loadtxt(kspecperp)[:,1]/float(numfiles)

                        kspecpara_avg += np.loadtxt(kspecpara)[:,1]/float(numfiles)


		mspec_avg = np.empty(np.shape(kspec_avg))
		mspecperp_avg = np.empty(np.shape(kspec_avg))
		mspecpara_avg = np.empty(np.shape(kspec_avg))
		mspec_avg[:] = np.nan
		mspecperp_avg[:] = np.nan
		mspecpara_avg[:] = np.nan
	else:
		# Averaging
                for ind in inds:
                        kspec = sorted(glob.glob(path+'kspectrum.'+str(ind)+'.txt'))[0]
                        mspec = sorted(glob.glob(path+'mspectrum.'+str(ind)+'.txt'))[0]

                        kspecperp = sorted(glob.glob(path+'kspecperp.'+str(ind)+'.txt'))[0]
                        mspecperp = sorted(glob.glob(path+'mspecperp.'+str(ind)+'.txt'))[0]

                        kspecpara = sorted(glob.glob(path+'kspecpara.'+str(ind)+'.txt'))[0]
                        mspecpara = sorted(glob.glob(path+'mspecpara.'+str(ind)+'.txt'))[0]

                        # Load data:
                        kspec_avg =+ np.loadtxt(kspec)[:,1]/float(numfiles)
                        mspec_avg =+ np.loadtxt(mspec)[:,1]/float(numfiles)

                        kspecperp_avg =+ np.loadtxt(kspecperp)[:,1]/float(numfiles)
                        mspecperp_avg =+ np.loadtxt(mspecperp)[:,1]/float(numfiles)

                        kspecpara_avg =+ np.loadtxt(kspecpara)[:,1]/float(numfiles)
                        mspecpara_avg =+ np.loadtxt(mspecpara)[:,1]/float(numfiles)

	Data_E['ks'] = np.loadtxt(kspec)[:,0]
	Data_E['kspec'] = kspec_avg
	Data_E['mspec'] = mspec_avg

	Data_E['kspecperp'] = kspecperp_avg
	Data_E['mspecperp'] = mspecperp_avg

	Data_E['kspecpara'] = kspecpara_avg
	Data_E['mspecpara'] = mspecpara_avg
        Data[run] = Data_E

# Saving data:
print "Saving Data"
name = 'rundat/SpectraAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))
