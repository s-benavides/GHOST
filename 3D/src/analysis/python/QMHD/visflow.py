import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys
import glob as glob
import field_calc
import h5py

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
runname = raw_input("Folder name: ")
path = '../'+runname+'/outs/'

if ('_3Re' or '_10Re') in runname:
    Nx=Ny=Nz=512
    Lx=Ly=Lz=2*np.pi
elif ('_4Lx') in runname:
    Ny=Nz = 256
    Nx = 1024
    Ly=Lz = 2*np.pi
    Lx = 8*np.pi
elif ('_2Lx') in runname:
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
print("Last output: %s" % int(tf[-1][0]))

outnum = raw_input("out num? ") #sys.argv[1]
outnum ="{:0>4s}".format(outnum)

# If any calculations need to be made:
otypes_calc = ['wz','wx','bdiss'] # also possible: bal_B0
for otype in otypes_calc:
        filelist = sorted(glob.glob(path+otype+'.'+outnum+'.out'))
        if len(filelist)==0:
                print("Need to calculate output for %s. Calculating..." % otype)
                field_calc.field_calc(runname,otype,outnum,Nx=Nx,Ny=Ny,Nz=Nz,Lx=Lx,Ly=Ly,Lz=Lz)

# otypes that we want to view or save
#otypes = ['wz']#raw_input("out type? ") #str(sys.argv[2])
otypes = ['wz','wx','bdiss']#raw_input("out type? ") #str(sys.argv[2])

level = int(Nx/2) #128 #input("level? ")#int(sys.argv[3])

# Reads binary files
outs  = dict([])
datbars_yz = dict([])
datbars_xy = dict([])
for otype in otypes:
    out = np.fromfile(path+otype+'.'+outnum+'.out').reshape(shape)#,order='F')
    outs[otype] = out
    omax_xy = np.max(out[level,:,:])
    omax_yz = np.max(out[:,:,level])
    omin_xy = np.min(out[level,:,:])
    omin_yz = np.min(out[:,:,level])
    datbars_xy[otype] = np.max([abs(omax_xy),abs(omin_xy)])
    datbars_yz[otype] = np.max([abs(omax_yz),abs(omin_yz)])
	#datbar = 100

# Show a horizontal cut of the field in the middle of the box
for otype in otypes:
	plt.figure()
	plt.imshow(outs[otype][level,:,:],vmin = -datbars_xy[otype],vmax = datbars_xy[otype],cmap='bwr')
	plt.title(runname+' '+otype+' out = '+str(outnum))
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.tight_layout()
	plt.savefig('./figs/'+runname+'_'+otype+'_'+str(outnum)+'_yx.png',bbox_inches='tight',dpi=300)
	
	plt.figure()
	plt.imshow(outs[otype][:,:,level],vmin = -datbars_yz[otype],vmax = datbars_yz[otype],cmap='bwr')
	plt.title(runname+' '+otype+' out = '+str(outnum))
	plt.colorbar()
	plt.xlabel('y')
	plt.ylabel('z')
	plt.tight_layout()
	plt.savefig('./figs/'+runname+'_'+otype+'_'+str(outnum)+'_zy.png',bbox_inches='tight',dpi=300)
plt.show()
