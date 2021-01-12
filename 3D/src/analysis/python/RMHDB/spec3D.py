import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys

# Reads spec2D files

# Path to the binary data
runname = raw_input("Folder name: ")
path = '../'+runname+'/outs/'

reso = 256 #input("Resolution? :")

# Spatial resolution
shape = (reso,reso,reso/2+1)

tf = np.loadtxt('../'+runname+'/run/time_field.txt')
print("Last output: %s" % int(tf[-1][0]))

outnum = raw_input("out num? ") #sys.argv[1]
outnum ="{:0>4s}".format(outnum)

otypes = ['mspec3D','kspec3D']
#otype = raw_input("out type? ") #str(sys.argv[2])

for otype in otypes:
	# Reads binary files
	#Opens file to read (and also accept binary output I think)
	f = open(path+otype+'.'+outnum+'.out','rb')         
	#Fortran files begin with integer telling you how many bytes are recorded in the write. We skip this integer and go straight to the data.
	f.seek(4)                            
	#gets the data in a 1x(nx*ny/nfile) vector, later to be transformed into 2D array, once stitched together
	out = np.fromfile(f,dtype = 'd',count=shape[0]*shape[1]*shape[2]).reshape(shape)        
	#out[out>0] = np.log10(out[out>0])
	#out[out==0] = out[out==0]+np.min(out)
	f.close()
	
	#Saves new output
	name = "./rundat/"+otype+"_"+runname+".out"
	print("Saving %s" % name)
	ofile = open(name,"wb")
	output = bytearray(out)
	ofile.write(output)
	
#
#omax = np.max(out)
#omin = np.min(out)
#datbar = np.max([abs(omax),abs(omin)])
#datbar = 100

# Show a horizontal cut of the field in the middle of the box
#plt.figure(1)
#plt.imshow(out[128,:,:],vmin = omin,vmax = omax,cmap='viridis')
#plt.colorbar()
#plt.ylim([-reso/3,reso/3])
#plt.xlim([0,reso/3])
#plt.xlabel(r'$k_x$')
#plt.ylabel(r'$k_z$')
#plt.tight_layout()
#plt.savefig(runname+'_'+otype+'.png')
#plt.show()
