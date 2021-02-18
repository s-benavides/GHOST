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
avglist = sorted(glob.glob('rundat/AvgTimeO*.txt'))

runnames = []
for file in avglist:
    run = file.split('rundat/AvgTime')[1]
    run = run.split('.txt')[0]
    if (('_2DF' in run)or('_old' in run)or('_nu' in run)or('O0N' in run)):
        pass
    else:
        runnames.append(run)

Data = dict([])
for i,run in enumerate(runnames):
    Data_E = dict([])
    print(" ---------------- Working on run %s ---------------- " % run)
    path = '../'+run+'/run/'

    # Some run info
    rand = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=142,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[5]

    sstep = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=142,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[3]

    cstep = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=142,skip_header=15,converters={2:  lambda val: float(val.translate(rule))},usecols=2)[4]

    omegaz = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=60,skip_header=125,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

    Nx, Nz = np.genfromtxt(path+'parameter.inp',comments='!',skip_footer=56,skip_header=130,converters={2:  lambda val: float(val.translate(rule))},usecols=2)

    print("omegaz = %s, Nx = %s, Nz = %s" % (omegaz,Nx,Nz))

    # Average start indices
    [start,start_fl,err_ind] = np.loadtxt('rundat/AvgTime'+run+'.txt')

    start = int(start)
    start_fl = int(start)
    err_ind = int(err_ind)

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

    # nu,kf
    nu = float(params_nu[4])
    hnu = float(params_nu[5])
    hek = float(params_nu[6])
    hok = float(params_nu[7])
    kdn = float(params_nu[2])
    kup = float(params_nu[3])
    kf = (kdn+kup)/2.

    #print('rand',rand,'sstep',sstep,'cstep',cstep,'nu',nu,'hnu',hnu,'hek',hek,'kup',kup)

    # list of observables to average:
    olist = {'enk':enk,'denk':denk,'henk':henk,'injk':injtot,'jenk':jenk,'ufk':ufk}
                 
    # AVERAGING
    Data_E = dict([])
    for ii,obs in enumerate(olist):
        avg = np.nanmean(olist[obs][start:])
        err = bunch_err.bunch_err(olist[obs][start:],err_ind=err_ind)
        Data_E[obs]=[avg,err]

    # Some calculations

    ##################################### WORKING HERE
    
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
    field_list = ['vx','vy','vz']#,'wx','wz']
    for field in field_list:
        filelist = sorted(glob.glob(path+field+'.'+outnum+'.out'))
        if len(filelist)==0:
            print("Need to calculate output for %s. Calculating..." % field)
            field_calc.field_calc(run,field,outnum,reso=reso)

        mean = 0.0
        out = np.fromfile(path+field+'.'+outnum+'.out').reshape(shape)
        mean = np.mean(np.abs(out)**2)
        Data_E[field]=mean

    mufk = np.sqrt(np.mean(ufk))
    Re_kf=np.sqrt(np.mean(ufk))/(nu*(kf)**(2*hek-1))
    print('run: %s, Re_rms: %f4' % (run,Re_kf))
    Ro_kf = (mufk*kf)/(2*omegaz)
    print("Ro(u(kf)) = %s" % Ro_kf)
    N_kf = (Nx**2+2*Nx*Nz+Nz**2)/(kf*mufk)
    print("N(u(kf)) = %s" % N_kf)

    u = ((4/5.)*np.mean(injtot)/kf)**(1/3.)
    Re_inj=u/(nu*(kf)**(2*hek-1))
    print('run: %s, Re_rms: %f4' % (run,Re_inj))
    Ro_inj = (u*kf)/(2*omegaz)
    print("Ro(inj) = %s" % Ro_inj)
    N_inj = (Nx**2+2*Nx*Nz+Nz**2)/(kf*u)
    print("N(inj) = %s" % N_inj)

    Data_E['nu'] = nu
    Data_E['hnu'] = hnu
    Data_E['hek'] = hek
    Data_E['hok'] = hok
    Data_E['kdn'] = kdn
    Data_E['kup'] = kup
    Data_E['kf'] = kf
    Data_E['Nx'] = Nx
    Data_E['Nz'] = Nz
    Data_E['omegaz'] = omegaz
    Data_E['rand'] = rand
    Data_E['Re_kf'] = Re_kf
    Data_E['Ro_kf'] = Ro_kf
    Data_E['N_kf'] = N_kf
    Data_E['Re_inj'] = Re_inj
    Data_E['Ro_inj'] = Ro_inj
    Data_E['N_inj'] = N_inj
    
    Data[run] = Data_E

# Saving data:
print "Saving Data"
name = 'rundat/DataAvg_'+str(date.today())+'.p'
pickle.dump(Data, open(str(name), 'wb'))	
