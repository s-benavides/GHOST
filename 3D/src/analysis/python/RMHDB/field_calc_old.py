import numpy as np
import glob as glob
import string as string
import matplotlib.pyplot as plt

# Reads binary files in a directory and does
# postprocessing (e.g., to visualize later with
# VAPOR). Note that for big runs, GHOST can do some
# automatic postprocessing (e.g., compute vorticity)
# at run time.
# Execute in ipython with '%run postprocess.py'

# Path to the binary data
runname = raw_input("Folder name: ")
path = '../'+runname+'/outs/'

tf = np.loadtxt('../'+runname+'/run/time_field.txt')
print("Last output: %s" % int(tf[-1][0]))

start = input('Starting at? ')
end = input('Ending at? ')

# Box size
Lx = 2*np.pi
Ly = 2*np.pi
Lz = 2*np.pi

reso = 256 #input("Resolution? :")

# Spatial resolution
NX = reso
NY = reso
NZ = reso
dx = Lx/NX
dy = Ly/NY
dz = Lz/NZ
shape = (NX,NY,NZ)

#otype = 'wz'
#otype = 'wx'
#otype = 'bz'
otype = raw_input("What output? wz,wx,bz? ")

# Reads binary files, computes vertical vorticity
# using centered finite differences, and saves in 
# a new binary file named 'wz.NNNN.out'
if 'wz' in otype:
	filelist = sorted(glob.glob(path+'vx.*.out'))
	for ii,file in enumerate(filelist):
	  ind = file.split(path+'vx.')[1]
	  ind = ind.split('.out')[0]
	  if start<=int(ind)<=end: 
		  print("working on file %s" % ind)
		  vx = np.fromfile(file).reshape(shape)#,order='F')
		  str = 'vy.'+ind+'.out'
		  vy = np.fromfile(path+str).reshape(shape)#,order='F')
		  adv = np.roll(vy,-1,axis=2)
		  ret = np.roll(vy,1,axis=2)
		  vy = (adv-ret)/(2*dx) # dv_y/dx
		  adv = np.roll(vx,-1,axis=1)
		  ret = np.roll(vx,1,axis=1)
		  vx = (adv-ret)/(2*dy) # dv_x/dy 
		  wz = vy-vx
		  wz.tofile(path+'wz.'+ind+'.out')
elif 'wx' in otype:
        filelist = sorted(glob.glob(path+'vz.*.out'))
        for ii,file in enumerate(filelist):
          ind = file.split(path+'vz.')[1]
          ind = ind.split('.out')[0]
          if start<=int(ind)<=end:
                  print("working on file %s" % ind)
                  vz = np.fromfile(file).reshape(shape)#,order='F')
                  str = 'vy.'+ind+'.out'
                  vy = np.fromfile(path+str).reshape(shape)#,order='F')
                  adv = np.roll(vy,-1,axis=0)
                  ret = np.roll(vy,1,axis=0)
                  vy = (adv-ret)/(2*dz) # dv_y/dz
                  adv = np.roll(vz,-1,axis=1)
                  ret = np.roll(vz,1,axis=1)
                  vz = (adv-ret)/(2*dy) # dv_z/dy
                  wx = vy-vz
                  wx.tofile(path+'wx.'+ind+'.out')
elif 'bz' in otype:
	filelist = sorted(glob.glob(path+'ax.*.out'))
        for ii,file in enumerate(filelist):
          ind = file.split(path+'ax.')[1]
          ind = ind.split('.out')[0]
	  if start<=int(ind)<=end: 
	          print("working on file %s" % ind)
	          vx = np.fromfile(file).reshape(shape,order='F')
	          str = 'ay.'+ind+'.out'
	          vy = np.fromfile(path+str).reshape(shape,order='F')
	          adv = np.roll(vy,-1,axis=0)
	          ret = np.roll(vy,1,axis=0)
	          vy = (adv-ret)/(2*dx) # dv_y/dx
	          adv = np.roll(vx,-1,axis=1)
	          ret = np.roll(vx,1,axis=1)
	          vx = (adv-ret)/(2*dy) # dv_x/dy
	          wz = vy-vx
	          wz.tofile(path+'bz.'+ind+'.out')
