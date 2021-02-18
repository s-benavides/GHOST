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

reso = 256 #input("Resolution? :")

# Spatial resolution
NX = reso
NY = reso
NZ = reso
shape = (NX,NY,NZ)

tf = np.loadtxt('../'+runname+'/run/time_field.txt')
print("Last output: %s" % int(tf[-1][0]))

outnum = raw_input("out num? ") #sys.argv[1]
outnum ="{:0>4s}".format(outnum)

# If any calculations need to be made:
otypes_calc = ['bx','by','bz','ME','bal_B0']
for otype in otypes_calc:
        filelist = sorted(glob.glob(path+otype+'.'+outnum+'.out'))
        if len(filelist)==0:
                print("Need to calculate output for %s. Calculating..." % otype)
                field_calc.field_calc(runname,otype,outnum,reso=reso)

# otypes that we want to view or save
otypes = ['ME','bal_B0']#raw_input("out type? ") #str(sys.argv[2])

level = 128 #input("level? ")#int(sys.argv[3])

# Reads binary files
outs  = dict([])
datbars = dict([])
for otype in otypes:
	out = np.fromfile(path+otype+'.'+outnum+'.out').reshape(shape)#,order='F')
	outs[otype] = out

	omax = np.max(out[:,:,level])
	omin = np.min(out[:,:,level])
	datbars[otype] = np.max([abs(omax),abs(omin)])
	#datbar = 100

# Show a horizontal cut of the field in the middle of the box
for otype in otypes:
	datbar = datbars[otype]
	plt.figure()
	plt.imshow(outs[otype][level,:,:],vmin = -datbar,vmax = datbar,cmap='bwr')
	plt.title(runname+' '+otype+' out = '+str(outnum))
	plt.colorbar()
	plt.xlabel('x')
	plt.ylabel('y')
	plt.tight_layout()
	plt.savefig('./figs/'+runname+'_'+otype+'_'+str(outnum)+'_yx.png',bbox_inches='tight',dpi=300)
	
	plt.figure()
	plt.imshow(outs[otype][:,:,level],vmin = -datbar,vmax = datbar,cmap='bwr')
	plt.title(runname+' '+otype+' out = '+str(outnum))
	plt.colorbar()
	plt.xlabel('y')
	plt.ylabel('z')
	plt.tight_layout()
	plt.savefig('./figs/'+runname+'_'+otype+'_'+str(outnum)+'_zy.png',bbox_inches='tight',dpi=300)
plt.show()
