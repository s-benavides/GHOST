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

def field_calc(run,otype,outnum,reso=256):
	# Path to the binary data
	path = '../'+run+'/outs/'
	
	# Box size
	Lx = 2*np.pi
	Ly = 2*np.pi
	Lz = 2*np.pi
	
	# Spatial resolution
	NX = reso
	NY = reso
	NZ = reso
	dx = Lx/NX
	dy = Ly/NY
	dz = Lz/NZ
	shape = (NX,NY,NZ)
	
	# Reads binary files, computes vertical vorticity
	# using centered finite differences, and saves in 
	# a new binary file named 'wz.NNNN.out'
	if 'wz' in otype:
		vx = np.fromfile(path+'vx.'+outnum+'.out').reshape(shape)	
		vy = np.fromfile(path+'vy.'+outnum+'.out').reshape(shape)	
		adv = np.roll(vy,-1,axis=2)
		ret = np.roll(vy,1,axis=2)
		vy = (adv-ret)/(2*dx) # dv_y/dx
		adv = np.roll(vx,-1,axis=1)
		ret = np.roll(vx,1,axis=1)
		vx = (adv-ret)/(2*dy) # dv_x/dy 
		wz = vy-vx
		wz.tofile(path+'wz.'+outnum+'.out')
	elif 'wx' in otype:
                vz = np.fromfile(path+'vz.'+outnum+'.out').reshape(shape)
                vy = np.fromfile(path+'vy.'+outnum+'.out').reshape(shape)
	        adv = np.roll(vy,-1,axis=0)
	        ret = np.roll(vy,1,axis=0)
	        vy = (adv-ret)/(2*dz) # dv_y/dz
	        adv = np.roll(vz,-1,axis=1)
	        ret = np.roll(vz,1,axis=1)
	        vz = (adv-ret)/(2*dy) # dv_z/dy
	        wx = vz-vy
	        wx.tofile(path+'wx.'+outnum+'.out')
        elif 'wy' in otype:
                vx = np.fromfile(path+'vx.'+outnum+'.out').reshape(shape)
                vz = np.fromfile(path+'vz.'+outnum+'.out').reshape(shape)
                adv = np.roll(vx,-1,axis=0)
                ret = np.roll(vx,1,axis=0)
                vx = (adv-ret)/(2*dz) # dv_x/dz
                adv = np.roll(vz,-1,axis=2)
                ret = np.roll(vz,1,axis=2)
                vz = (adv-ret)/(2*dx) # dv_z/dx
                wy = vx-vz
                wy.tofile(path+'wy.'+outnum+'.out') # dv_x/dz - dv_z/dx
	elif 'bz' in otype:
                vx = np.fromfile(path+'ax.'+outnum+'.out').reshape(shape)
                vy = np.fromfile(path+'ay.'+outnum+'.out').reshape(shape)
                adv = np.roll(vy,-1,axis=2)
                ret = np.roll(vy,1,axis=2)
                vy = (adv-ret)/(2*dx) # dv_y/dx
                adv = np.roll(vx,-1,axis=1)
                ret = np.roll(vx,1,axis=1)
                vx = (adv-ret)/(2*dy) # dv_x/dy
                wz = vy-vx
                wz.tofile(path+'bz.'+outnum+'.out')
        elif 'bx' in otype:
                vz = np.fromfile(path+'az.'+outnum+'.out').reshape(shape)
                vy = np.fromfile(path+'ay.'+outnum+'.out').reshape(shape)
                adv = np.roll(vy,-1,axis=0)
                ret = np.roll(vy,1,axis=0)
                vy = (adv-ret)/(2*dz) # dv_y/dz
                adv = np.roll(vz,-1,axis=1)
                ret = np.roll(vz,1,axis=1)
                vz = (adv-ret)/(2*dy) # dv_z/dy
                wx = vz-vy
                wx.tofile(path+'bx.'+outnum+'.out')
        elif 'by' in otype:
                vx = np.fromfile(path+'ax.'+outnum+'.out').reshape(shape)
                vz = np.fromfile(path+'az.'+outnum+'.out').reshape(shape)
                adv = np.roll(vx,-1,axis=0)
                ret = np.roll(vx,1,axis=0)
                vx = (adv-ret)/(2*dz) # dv_x/dz
                adv = np.roll(vz,-1,axis=2)
                ret = np.roll(vz,1,axis=2)
                vz = (adv-ret)/(2*dx) # dv_z/dx
                wy = vx-vz
                wy.tofile(path+'by.'+outnum+'.out') # dv_x/dz - dv_z/dx
        elif 'jz' in otype:
                vx = np.fromfile(path+'bx.'+outnum+'.out').reshape(shape)
                vy = np.fromfile(path+'by.'+outnum+'.out').reshape(shape)
                adv = np.roll(vy,-1,axis=2)
                ret = np.roll(vy,1,axis=2)
                vy = (adv-ret)/(2*dx) # dv_y/dx
                adv = np.roll(vx,-1,axis=1)
                ret = np.roll(vx,1,axis=1)
                vx = (adv-ret)/(2*dy) # dv_x/dy
                wz = vy-vx
                wz.tofile(path+'jz.'+outnum+'.out')
        elif 'jx' in otype:
                vz = np.fromfile(path+'bz.'+outnum+'.out').reshape(shape)
                vy = np.fromfile(path+'by.'+outnum+'.out').reshape(shape)
                adv = np.roll(vy,-1,axis=0)
                ret = np.roll(vy,1,axis=0)
                vy = (adv-ret)/(2*dz) # dv_y/dz
                adv = np.roll(vz,-1,axis=1)
                ret = np.roll(vz,1,axis=1)
                vz = (adv-ret)/(2*dy) # dv_z/dy
                wx = vz-vy
                wx.tofile(path+'jx.'+outnum+'.out')
        elif 'jy' in otype:
                vx = np.fromfile(path+'bx.'+outnum+'.out').reshape(shape)
                vz = np.fromfile(path+'bz.'+outnum+'.out').reshape(shape)
                adv = np.roll(vx,-1,axis=0)
                ret = np.roll(vx,1,axis=0)
                vx = (adv-ret)/(2*dz) # dv_x/dz
                adv = np.roll(vz,-1,axis=2)
                ret = np.roll(vz,1,axis=2)
                vz = (adv-ret)/(2*dx) # dv_z/dx
                wy = vx-vz
                wy.tofile(path+'jy.'+outnum+'.out') # dv_x/dz - dv_z/dx
	elif 'KE' in otype:
                vx = np.fromfile(path+'vx.'+outnum+'.out').reshape(shape)
                vy = np.fromfile(path+'vy.'+outnum+'.out').reshape(shape)
                vz = np.fromfile(path+'vz.'+outnum+'.out').reshape(shape)
		KE = 0.5*(vx**2+vy**2+vz**2)
                KE.tofile(path+'KE.'+outnum+'.out')
	elif 'ME' in otype:
                vx = np.fromfile(path+'bx.'+outnum+'.out').reshape(shape)
                vy = np.fromfile(path+'by.'+outnum+'.out').reshape(shape)
                vz = np.fromfile(path+'bz.'+outnum+'.out').reshape(shape)
                KE = 0.5*(vx**2+vy**2+vz**2)
                KE.tofile(path+'ME.'+outnum+'.out')
	elif 'bal_B0' in otype:		
                vx = np.fromfile(path+'vx.'+outnum+'.out').reshape(shape)
                vy = np.fromfile(path+'vy.'+outnum+'.out').reshape(shape)
                vz = np.fromfile(path+'vz.'+outnum+'.out').reshape(shape)
                bx = np.fromfile(path+'bx.'+outnum+'.out').reshape(shape)
                by = np.fromfile(path+'by.'+outnum+'.out').reshape(shape)
                bz = np.fromfile(path+'bz.'+outnum+'.out').reshape(shape)
		adv = np.roll(vx,-1,axis=2)
                ret = np.roll(vx,1,axis=2)
                vx = (adv-ret)/(2*dx) # dv_x/dx
                adv = np.roll(vy,-1,axis=2)
                ret = np.roll(vy,1,axis=2)
                vy = (adv-ret)/(2*dx) # dv_y/dx
                adv = np.roll(vz,-1,axis=2)
                ret = np.roll(vz,1,axis=2)
                vz = (adv-ret)/(2*dx) # dv_z/dx
		bx = bx*vx+by*vy+bz*vz
		bx.tofile(path+'bal_B0.'+outnum+'.out')
	return
