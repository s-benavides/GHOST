import numpy as np
import matplotlib.pyplot as plt
import scipy
import sys
import glob as glob
import string

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Making rule for reading nui,nun,mu
rule = string.maketrans('d', '0')

# Get all the AvgTimeO*B*.txt files
avglist = sorted(glob.glob('rundat/*B20T90.txt'))

runnames = []
for file in avglist:
  run = file.split('rundat/AvgTime')[1]
  run = run.split('.txt')[0]
  runnames.append(run)


reso = 256 #input("Resolution? :")

# Spatial resolution
NX = reso
NY = reso
NZ = reso
shape = (NX,NY,NZ)

means = []
omegas =[]
for runname in runnames:
	print("Working on %s" % runname)
	path = '../'+runname+'/outs/'
	tf = np.loadtxt('../'+runname+'/run/time_field.txt')
	print("Last output: %s" % int(tf[-1][0]))
	
	outnum = int(tf[-1][0]) #raw_input("out num? ") #sys.argv[1]
	outnum ="{:0>4s}".format(outnum)
	omega = #input("omega? ")	
	omegas.append(omega)

	# Reads binary files
	#psi = np.fromfile(path+'ps.'+outnum+'.out',dtype=np.float32).reshape(shape,order='F')
	out = np.fromfile(path+'vx.'+outnum+'.out').reshape(shape)#,order='F')
	mean = 0.0
	nums = 0
	for i in range(reso):
		for j in range(reso):
			for k in range(reso):
				mean+=(bs(out[i,j,k])**2)2.
				nums+=1
	print(nums,256**3)
	means.append(mean/nums)	
	

plt.figure()
plt.plot(omegas,means,'ok')
plt.show()
