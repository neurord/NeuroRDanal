#sig2-fig.py
#Generate signature traces, plot them, and create files for plotting in othe software
#ARGS="subdir/fileroot,LTPmol,LTDmol,tstart tend" where tstart and tend used to calculate basal value
#ARGS="Model_SPNspineAChm4R_Gshydr5_GapD,Aphos CKpCam Epac1 Pkc,2ag,10 20"
from __future__ import print_function
from __future__ import division

import os
import numpy as np
from matplotlib import pyplot
import sys  
import glob
from NeuroRDanal import h5utils
from NeuroRDanal import plot_h5 as pu5
from collections import OrderedDict
import plot_h5 as pu5

#specify either 'mean' or column numbers, or region names
#These should probably be additional parameters/arguments
#coltype=[1,2,3,4,5,6]
coltype=['nonspine','sa1[0]']#,'sa1[1]','sa1[2]','sa1[3]','sa1[4]','sa1[5]','sa1[6]','sa1[7]','sa1[8]','sa1[9]']
trials=3 #make trials=1 if coltype=='mean'
normYN=1
textsize=8
plotYN=1

try:
    args = ARGS.split(",")
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

#Identify set of filenames and parameter values
if len(args[1]):
    pattern=args[0]+'-'+args[1]
else:
    pattern=args[0]
fnames = glob.glob(pattern+"*plas.txt")
print(pattern+"*plas.txt, ", len(fnames), 'files found')
fname_roots=sorted(set([fnm[0:fnm.rfind('_')] for fnm in fnames]))
fname_endings=set([fnm[fnm.rfind('_'):] for fnm in fnames])
parval=[fnm[fnm.find('-')+1:].replace('-',' ') for fnm in fname_roots]

#assign other inputs (args) to parameters
ltp_molecules=args[2].split()
ltd_molecules=args[3].split()
tstart,tend=args[4].split()

all_molecules=ltp_molecules+ltd_molecules
num_mols=len(all_molecules)

#Read in all files and assemble relevant columns into large data matrix
all_peaks={}
all_sig_array={}
for molnum,mol in enumerate(all_molecules):
    for fnum,fname in enumerate(fname_roots):
        temp=[]
        fnm=fname+'_'+mol+'plas.txt'
        f = open(fnm, 'r+')
        header=f.readline()
        f.close()
        head_names=header.split()
        if coltype=='mean':
            col_num=[head_names.index(x) for x in head_names if coltype in x][0:-1]
        elif type(coltype) is list and type(coltype[0]) is str:
            for col in coltype:
                temp.append([head_names.index(x) for x in head_names if '_t' in x and col in x])
            col_num=[val for sublist in temp for val in sublist]
        else:
            col_num=coltype
        alldata=np.loadtxt(fnm,skiprows=1)
        time=alldata[:,0]
        data_cols=alldata[:,col_num]
        #very specific kluge because one of the files started much later. 
        if fname.split('-')[-1]=='blockPKA' or fname.split('-')[-1]=='blockPKA2':
            extra=int(30/time[1])
        else:
            extra=0
        pstrt=int(int(tstart)/time[1])
        pend=int(int(tend)/time[1])
        if fnum==0:
            sig_array=np.zeros((len(fname_roots),len(time),len(col_num)))
            peaks=np.zeros((len(fname_roots),len(col_num)))
        #
        ########################### normalize and create arrays of molecules ##########################
        #
        if normYN:
            basal=np.mean(data_cols[pstrt+extra:pend+extra],axis=0)
        else:
            basal=np.zeros(len(col_num))
        if np.shape(data_cols[extra:,:])[0]<np.shape(sig_array)[1]:
            print (fname,'zero padding end by', extra*time[1])
            sig_array[fnum]=np.pad(data_cols[extra:,:],((0,extra),(0,0)),mode='constant')-basal
        else:
            sig_array[fnum]=data_cols[extra:,:]-basal
        peaks[fnum]=np.max(sig_array[fnum],axis=0)
    all_sig_array[mol]=sig_array
    all_peaks[mol]=peaks

########################### Signature part ##########################
#place into function  (or two)
sig_ltp=np.zeros((len(fname_roots),len(time),len(col_num)))
sig_ltd=np.zeros((len(fname_roots),len(time),len(col_num)))

#for each file, calculate two signatures:
#   (A): discrim1Const+discrim11*mol1trace+discrim12*mol2trace etc.
#        apply to each trial and repeat for discrim2
#        1st normalize by baseline subtraction (and baseline division for non-zero molecules)


#   (B): sum normalized molecules separately for dend and spine, separately for LTP and LTD
#        apply to each trial, 1st normalize by baseline subtraction (and peak factor)
for molnum,mol in enumerate(all_molecules):
    region_peak=np.max(all_peaks[mol],axis=0)
    for reg in range(num_regions):
        colset=[tr+reg*trials for tr in range(trials)]
        region_peak[colset]=np.max(region_peak[colset])
    print (mol,region_peak)
    for fnum,fname in enumerate(fname_roots):
        if mol in ltp_molecules:
            sig_ltp[fnum]=sig_ltp[fnum]+all_sig_array[mol][fnum]/region_peak
        else:
            sig_ltd[fnum]=sig_ltd[fnum]+all_sig_array[mol][fnum]/region_peak
#
#show figures (and create output files for publication figures)
figtitle=fnm[0:fnm.find('-')]
label=[[] for p in parval]
domain=[head_names[x].split(mol) for x in col_num]
for par in range(len(parval)):
    label[par]=[parval[par]+' '+dom[0][0:6]+dom[1] for dom in domain]
sign_title=args[2]+' vs '+args[3]

if plotYN:
  if len(ltd_molecules):
    pu5.plot_signature(label,sig_ltp,time,figtitle,sign_title,textsize,thresh,sig_ltd)
  else:
    pu5.plot_signature(label,sig_ltp,time,figtitle,sign_title,textsize,thresh)

        
    
