import os,glob,pathlib
import numpy as np

dirs='./'

skip = False

tf_paths = list(pathlib.Path().glob(dirs+'/O*/run/time_field.txt'))

for path in tf_paths:
#	if (('O50B40T90' in str(path)) or ('O10B40T90' in str(path))):
#		print("Skipping  %s" % os.path.split(path)[0][:-4])
    cont=True
    if skip:
        print("Skipping")
        pass
    else:
        run = os.path.split(path)[0][:-4]
        print("Working on %s" % os.path.split(path)[0][:-4])
        # Find last output
        if np.size(np.loadtxt(str(path))[-1])==1:
            print("Run just started, skipping.")
        else:
            last_out = int(np.loadtxt(str(path))[-1][0])
            last_out = "{0:0>4s}".format(str(last_out))
            print("Last out: %s" % last_out)
            outsdir = os.path.split(path)[0][:-3]+"outs"
            insdir = os.path.split(path)[0][:-3]+"ins"

        # Fields to remove:
        outs = ['vx','vy','vz','wx','wy','wz']
        #outs = ['wx','wy','wz','bx','by','bz','jx','jy','jz']

        #for deldir in [outsdir,insdir]:
    ############# OUTPUT folder
        deldir = outsdir
        print(" --- Deleting in %s --- " % deldir)
        # Fields:
        count = 0
        for out in outs:
            fdel = list(set(glob.glob(deldir+'/'+out+'.*.out'))-set(glob.glob(deldir+'/'+out+'.'+last_out+'.out')))
            lastfiles = list(set(glob.glob(deldir+'/'+out+'.'+last_out+'.out')))
            # Check if any of the last files are not the right size
            for ffile in lastfiles:
                osize = pathlib.Path(ffile).stat().st_size
                #print(run,osize)
                if ((osize!=134217728)&(osize!=268435456)&(osize!=1073741824)&(osize!=536870912)):
                    print(" ~%~%~%~%~%~ Not complete file! Not deleting and stopping.  ~%~%~%~%~%~ ")
                    cont = False
            # If it's all good, then deletes.
            for ffile in fdel:
                if not cont:
                    print('not deleting %s' % ffile)
                else:
                    os.remove(ffile)
                    count+=1
        if count>0:
            print("Done! Removed %s out files" % count)

############# INPUT folder
#	deldir = insdir
#	print(" --- Deleting in %s --- " % deldir)
#	# Fields:
#	count=0
#	fdel = glob.glob(deldir+'/*.out')
#	for ffile in fdel:
#		os.remove(ffile)
#		count+=1
#	if count>0:
#		print("Done! Removed %s out files" % count)

