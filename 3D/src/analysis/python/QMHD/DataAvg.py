import numpy as np
import matplotlib.pyplot as plt
#from matplotlib import rc
import pickle
import glob as glob
import string
from datetime import date
import bunch_err
import field_calc

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Get all the AvgTimeO*B*.txt files
avglist = sorted(glob.glob('rundat/*T90*.txt'))

runnames = []
for file in avglist:
  run = file.split('rundat/AvgTime')[1]
  run = run.split('.txt')[0]
  if (('O50B10' in run)or('_rnd' in run)or('Pr_0d001' in run)or('_c' in run)or('movie' in run)):
        pass
  else:
        runnames.append(run)

# Drops runs that don't have perp/para outputs
no_ebdata = []
for run in runnames:
  path = '../'+run+'/run/'
  if not glob.glob(path+'ebperp_ebpara.txt'):
    print("RUN %s doesn't have ebperp_ebpara data!" % run)
    no_ebdata.append(run)

Data = dict([])
for i,run in enumerate(runnames):
	Data_E = dict([])
	print("Working on run %s " % run)
	path = '../'+run+'/run/'

	# Some run info
	rinfo  = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=136,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

	rand = rinfo[5]

	sstep = rinfo[3]

	cstep = rinfo[4]

	# Average start indices
	[start,start_fl,err_ind] = np.loadtxt('rundat/AvgTime'+run+'.txt')

	start = int(start)
	start_fl = int(start)
	err_ind = int(err_ind)

	# Reading rotation and B-field
	b0x,b0y,b0z = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=65,skip_header=110,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

	omegax,omegay,omegaz = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=56,skip_header=123,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

	# Reads balance.txt
        t,enk,enm,denk,denm,henk,henm,jxb = np.loadtxt(path+'balance.txt',unpack=True)
        t2,injk,injm = np.loadtxt(path+'injection.txt',unpack=True)
        t3,ufk,ufm = np.loadtxt(path+'uf.txt',unpack=True)
	t4,eperp_B,epara_B,eperp_Omega,epara_Omega = np.loadtxt(path+'eperp_epara.txt',unpack=True)
	if run not in no_ebdata:
		t5,ebperp_B,ebpara_B,ebperp_Omega,ebpara_Omega = np.loadtxt(path+'ebperp_ebpara.txt',unpack=True)
	else:
		t5,ebperp_B,ebpara_B,ebperp_Omega,ebpara_Omega = np.loadtxt(path+'eperp_epara.txt',unpack=True)*np.nan
	t6,khel,mhel = np.loadtxt(path+'helicity.txt',unpack=True)
        params_nu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=110,skip_header=41,converters={2:  lambda val: float(val.translate(rule))},usecols=2)
        params_mu = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=95,skip_header=61,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

        if len(t)!=len(t2):
                tmin = np.min([len(t),len(t2)])
        else:
                tmin=len(t)

        if rand!=0:
		print("randomi forcing")
                injk = injk/2.
                injm = injm/2. # RANDOM FORCING

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

	# list of observables to average:
	olist = {'enk':enk,'enm':enm,'denk':denk,'denm':denm,'henk':henk,'henm':henm,'jxb':jxb,'injk':injk,'injm':injm,'ufk':ufk,'ufm':ufm,'eperp_B':eperp_B,'epara_B':epara_B,'eperp_Omega':eperp_Omega,'epara_Omega':epara_Omega,'ebperp_B':ebperp_B,'ebpara_B':ebpara_B,'ebperp_Omega':ebperp_Omega,'ebpara_Omega':ebpara_Omega,
		'khel':khel,'mhel':mhel}
	
        # AVERAGING
	Data_E = dict([])
	for ii,obs in enumerate(olist):
		avg = np.nanmean(olist[obs][start:])
		err = bunch_err.bunch_err(olist[obs][start:],err_ind=err_ind)
		Data_E[obs]=[avg,err]

        # Some calculations

	# vx^2, vy^2, vz^2
	reso = 256
	# Spatial resolution
	NX = reso
	NY = reso
	NZ = reso
	shape = (NX,NY,NZ)

	path = '../'+run+'/outs/'
        tf = np.loadtxt('../'+run+'/run/time_field.txt')
        outnum = str(int(tf[-1][0])) #raw_input("out num? ") #sys.argv[1]
        outnum ="{:0>4s}".format(outnum)

        # Reads binary files
        #psi = np.fromfile(path+'ps.'+outnum+'.out',dtype=np.float32).reshape(shape,order='F')
	field_list = ['vx','vy','vz']#,'ax','ay','az','bx','by','bz','jx','jy','jz']
	for field in field_list:
        	filelist = sorted(glob.glob(path+field+'.'+outnum+'.out'))
        	if len(filelist)==0:
       	 	        print("Need to calculate output for %s. Calculating..." % field)
	                field_calc.field_calc(run,field,outnum,reso=reso)

		mean = 0.0
		out = np.fromfile(path+field+'.'+outnum+'.out').reshape(shape)
        	mean = np.sum(np.abs(out)**2)/float(reso)**3
		Data_E[field]=mean		

#        Re_rms=np.sqrt(avg_ufk)/(nu*(kf)**(2*hek-1))
        Re_rms=np.sqrt(Data_E['ufk'][0])/(nu*(kf)**(2*hek-1))
        Rem_rms=np.sqrt(Data_E['ufm'][0])/(nu*(kf)**(2*hek-1))
        Pr = nu/mu

	Data_E['Re_rms'] = Re_rms
	Data_E['Pr'] = Pr
	Data_E['nu'] = nu
	Data_E['hnu'] = hnu
	Data_E['mu'] = mu
	Data_E['hmu'] = hmu
	Data_E['hek'] = hek
	Data_E['hok'] = hok
	Data_E['hem'] = hem
	Data_E['hom'] = hom
	Data_E['kdn'] = kdn
	Data_E['kup'] = kup
	Data_E['kf'] = kf
	Data_E['b0x'] = b0x
	Data_E['b0y'] = b0y
	Data_E['b0z'] = b0z
	Data_E['omegax'] = omegax
	Data_E['omegay'] = omegay
	Data_E['omegaz'] = omegaz
	Data_E['rand'] = rand
		
	Data[run] = Data_E
	
# Saving data:
print "Saving Data"
name = 'rundat/DataAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))	
