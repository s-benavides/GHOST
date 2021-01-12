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
shape = (reso,reso/2+1)

tf = np.loadtxt('../'+runname+'/run/time_spec.txt')
print("Last output: %s" % int(tf[-1][0]))

outnum = raw_input("out num? ") #sys.argv[1]
outnum ="{:0>4s}".format(outnum)

otype = 'kspec2D_yavg' #raw_input("out type? ") #str(sys.argv[2])

# Reads binary files
f = open(path+otype+'.'+outnum+'.out','rb')         #Opens file to read (and also accept binary output I think)
f.seek(4)                            #Fortran files begin with integer telling you how many bytes are recorded in the write. We skip this integer and go straight to the data.
out = np.fromfile(f,dtype = 'd',count=shape[0]*shape[1]).reshape(shape)        #gets the data in a 1x(nx*ny/nfile) vector, later to be transformed into 2D array, once stitched together
#print(np.min(out),np.max(out))
#out[out>0] = np.log10(out[out>0])
#out[out==0] = out[out==0]+np.min(out)
f.close()

out = np.flip(out,0)

#Saves new output
#ofile = open("test.out","wb")
#output = bytes(out)
#ofile.write(output)

omax = np.max(out)
omin = np.min(out)
datbar = np.max([abs(omax),abs(omin)])
#datbar = 100

# Show a horizontal cut of the field in the middle of the box
plt.figure(1)
plt.imshow(out,interpolation='none',vmin = omin,vmax = omax,cmap='viridis',extent=[0,reso/2,-reso/2,reso/2])
plt.colorbar()
plt.ylim([-15,15])
plt.xlim([0,15])
plt.xlabel(r'$k_x$')
plt.ylabel(r'$k_z$')
plt.tight_layout()
plt.savefig('./figs/'+runname+'_'+otype+'.png')
plt.show()
