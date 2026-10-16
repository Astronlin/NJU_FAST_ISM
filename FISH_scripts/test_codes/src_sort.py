## First version: Lingrui Lin; Nov 27, 2021
## Sort sources of each beam by srcname

import numpy as np
import os
import sys
import glob

filelist = sys.argv[1]
fitslist = np.loadtxt(filelist,dtype='str')
srclist  = [fitsfile.split('/')[-1][:8] for fitsfile in fitslist]
srclist  = np.unique(srclist)
path     = sys.argv[2]

for i in range(len(srclist)):
    srcname  = srclist[i]
    dirname  = '{}/src_{}'.format(path,srcname)
    os.system(r'mkdir '+dirname)
    srcfiles = ' '.join(glob.glob(path+'/{}*'.format(srcname)))
    os.system('mv {} {}'.format(srcfiles,dirname))

#for i in range(1,20):
#    os.system('python3 src_sort.py M{}.LIST /Leo_disk/lingrui/my_projects/FAST_explore/Tcal_files/all_Tcal/zzy_all/M{}'.format(str(i).rjust(2,'0'),str(i).rjust(2,'0'))) 
#
#for i in range(1,20):
#    os.system('ls $PWD/M{}/*/* > M{}.LIST'.format(str(i).rjust(2,'0'),str(i).rjust(2,'0')))  
