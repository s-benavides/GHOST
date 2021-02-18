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
avglist = sorted(glob.glob('rundat/*T90.txt'))
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
        [start,start_fl,err_ind] = np.loadtxt('rundat/AvgTime'+run+'.txt')

        start = int(start)
        start_fl = int(start_fl)

	rinfo  = np.genfromtxt('../'+run+'/run/parameter.inp',comments='!',skip_footer=136,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

        rand = rinfo[5]	

	inds = []
        files = sorted(glob.glob(path+'ktransfer.*.txt'))
        for file in files:
                num = file.split(path+'ktransfer.')[1]
                num = num.split('.txt')[0]
                if (int(num)>=start_fl):
                        inds.append(num)

        numfiles = int(len(inds))
        print('%s files to average' % numfiles)


        # Reads flux files
	# Load data
	t,injk,injm = np.loadtxt('../'+run+'/run/injection.txt',unpack=True)

	injtot = injk+injm

	injtot = np.mean(injtot[start:])

        if rand!=0:
        	injtot = 0.5*injtot

	fields = ['ktransfer','mtransfer','jtransfer','ktranpara','mtranpara','jtranpara','ktranperp','mtranperp','jtranperp']

	# Averaging
	for field in fields:
#		print(field)
        	for ii,ind in enumerate(inds):
	        	# Load file names
                	flux = sorted(glob.glob(path+field+'.'+str(ind)+'.txt'))[0]
			# Load and Average File Names:
			if ii==0:
                		flux_avg = np.loadtxt(flux)[:,1]/float(numfiles)/injtot
			else:
				try: # Sometimes the saved file is corrupt if computation stopped in the middle of a save
                			flux_avg += np.loadtxt(flux)[:,1]/float(numfiles)/injtot
				except (KeyboardInterrupt, SystemExit):
					raise
				except: print("Error with file %s" % flux)
		Data_E[field] = flux_avg
	
	Data_E['ks'] = np.loadtxt(flux)[:,0]
        Data[run] = Data_E

# Saving data:
print "Saving Data"
name = 'rundat/FluxesAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))
