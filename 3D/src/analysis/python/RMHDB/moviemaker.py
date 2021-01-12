import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import glob as glob
import sys
import imageio
#imageio.plugins.ffmpeg.download()

# Reads a binary file and plots a cut in the x-y plane.
# Execute in ipython with '%run plot_bindata.py'

# Path to the binary data
runname = raw_input("Folder name: ")
path = '../'+runname+'/outs/'

# Out type (e.g. 'wz','ax','spec2D_yavg',...)
#otype = 'kspec2D_yavg'
otype = 'vx'#'wz'
print("Making movie of %s for run %s" % (otype,runname))

# Spatial resolution
N = 256

if otype in ['kspec2D_yavg']:
	shape = (int(N/2+1),int(N/2+1))
else:
	shape = (int(N),int(N),int(N))
	plane = 'xy'
	#plane = 'xz'
	level = int(N/2.6)

# File path
filelist = sorted(glob.glob(path+otype+'.*.out'))
nfiles = np.size(filelist)
print("nfiles = %s" % nfiles)
 

if otype in ['kspec2D_yavg']:
	frames_per_second = int(nfiles/float(5))  # so it lasts 5 seconds
	writer = imageio.get_writer('./movies/'+runname+'_'+otype+'.wmv', codec='msmpeg4',quality=10.0, fps = frames_per_second)
else:
	frames_per_second = 5
	writer = imageio.get_writer('./movies/'+runname+'_'+otype+'_'+plane+'_plane.wmv', codec='msmpeg4',quality=10.0, fps = frames_per_second)


# Reads binary files
for ii in range(nfiles):
	print("Working on snapshot %s" % (ii+1))
	# Reads binary files
	if otype in ['kspec2D_yavg']:
		#Opens file to read (and also accept binary output I think)
		f = open(filelist[ii],'rb')
		#Fortran files begin with integer telling you how many bytes are recorded in the write. We skip this integer and go straight to the data.
		f.seek(4)
		#gets the data in a 1x(nx*ny/nfile) vector, later to be transformed into 2D array, once stitched together
		out = np.fromfile(f,dtype = 'd',count=shape[0]**2).reshape(shape)  
		out[out>0] = np.log10(out[out>0])
		out[out==0] = out[out==0]+np.min(out)
		f.close()
	else:
		out = np.fromfile(filelist[ii]).reshape(shape,order='F')
#	if ii==0:
#                if otype not in ['kspec2D_yavg']:
#                        outd = np.fromfile(filelist[-1]).reshape(shape,order='F')
#                        datmin = np.amin(outd)
#                        datmax = np.amax(outd)
#                        datbar = max(abs(datmin),abs(datmax))
#                        print(datmin,datmax)
#                else:
	datmin = np.amin(out)
	datmax = np.amax(out)
	datbar = max(abs(datmin),abs(datmax))
	fig = plt.figure(1)
	if otype in ['kspec2D_yavg']:
		im = plt.imshow(out,cmap=cm.viridis,vmin=datmin,vmax=datmax)
	else:
		if 'xy' in plane:
			im = plt.imshow(out[level,:,:],cmap=cm.bwr,vmin=-datbar,vmax=datbar)
		elif 'xz' in plane:
			im = plt.imshow(out[:,level,:],cmap=cm.bwr,vmin=-datbar,vmax=datbar)
	cbar = plt.colorbar(im)
	if otype in ['kspec2D_yavg']:
		plt.xlabel(r"$k_\perp$")
		plt.ylabel(r"$k_\parallel$")
		plt.ylim([0,N/3])
		plt.xlim([0,N/3])
		plt.title("KE, run %s, out # %s" % (runname,ii+1))
	else:
		if 'xy' in plane:
			plt.xlabel('x')
			plt.ylabel('y')
			plt.title("%s, run %s, z = %s, out # %s" % (otype,runname,level,ii+1))
		elif 'xz' in plane:
			plt.xlabel('x')
			plt.ylabel('z')
			plt.title("%s, run %s, y = %s, out # %s" % (otype,runname,level,ii+1))
	plt.tight_layout()
	fig.canvas.draw()
	img = np.fromstring(fig.canvas.tostring_rgb(), dtype=np.uint8, sep='')
	img = img.reshape(fig.canvas.get_width_height()[::-1] + (3,))
	writer.append_data(img)
	plt.close()
if otype in ['kspec2D_yavg']:
	print("done making %s for run %s" % (otype,runname))
else:
	print("done making %s in plane %s for run %s" % (otype,plane,runname))
