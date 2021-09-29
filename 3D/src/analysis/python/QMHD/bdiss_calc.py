import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys
import glob as glob
import field_calc
import h5py

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

runs1 = list(glob.glob('../*_bdiss/'))
runs = []
for run in runs1:
    run_t = str(run).split('/')[1]
    runs.append(run_t)
print(runs)
for runname in runs:
    #if '_10Re' not in runname:
    print("Working on %s" % runname)
    # Path to the binary data
    #runname = raw_input("Folder name: ")
    path = '../'+runname+'/outs/'

    if ('_10Re') in str(runname):
        print('10Re run!')
        Nx=Ny=Nz=512
        Lx=Ly=Lz=2*np.pi
    elif ('_4Lx') in str(runname):
        Ny=Nz = 256
        Nx = 1024
        Ly=Lz = 2*np.pi
        Lx = 8*np.pi
    elif ('_2Lx') in str(runname):
        Ny=Nz = 256
        Nx = 512
        Ly=Lz = 2*np.pi
        Lx = 4*np.pi
    else:
        Nx=Ny=Nz=256
        Lx=Ly=Lz=2*np.pi

    # Spatial resolution
    shape = (Nx,Ny,Nz)

    tf = np.loadtxt('../'+runname+'/run/time_field.txt')
    #print("Last output: %s" % int(tf[-1][0]))
    #outnum = raw_input("out num? ") #sys.argv[1]
    outnum = str(int(tf[-1][0]))
    outnum ="{:0>4s}".format(outnum)

    # If any calculations need to be made:
    otypes_calc = ['bdiss'] # also possible: bal_B0
    for otype in otypes_calc:
            filelist = sorted(glob.glob(path+otype+'.'+outnum+'.out'))
            if len(filelist)==0:
                    print("Need to calculate output for %s. Calculating..." % otype)
                    field_calc.field_calc(runname,otype,outnum,Nx=Nx,Ny=Ny,Nz=Nz,Lx=Lx,Ly=Ly,Lz=Lz)

