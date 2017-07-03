#sig2.py
#evaluate various features to use for signature
#ARGS="subdir/fileroot,LTPmol,LTDmol,tstart tend"
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

coltype='mean'
normYN=0
textsize=8

try:
    args = ARGS.split(",")
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

fnames = glob.glob(args[0]+"*plas.txt")
print(args[0]+"*plas.txt, ", len(fnames), 'files found')
fname_roots=sorted(set([fnm[0:fnm.rfind('_')] for fnm in fnames]))
fname_endings=set([fnm[fnm.rfind('_'):] for fnm in fnames])
parval=[fnm[fnm.find('-')+1:].replace('-',' ') for fnm in fname_roots]

ltp_molecules=args[1].split()
ltd_molecules=args[2].split()
num_ltpmols=len(ltp_molecules)
num_ltdmols=len(ltd_molecules)
all_molecules=ltp_molecules+ltd_molecules
num_mols=len(all_molecules)

tstart,tend=args[3].split()

if len(args[4]):
    thresh=args[4].split()
else:
    thresh=['0', '0', '0', '0']

all_peaks={}
all_sig_array={}
for molnum,mol in enumerate(all_molecules):
    for fnum,fname in enumerate(fname_roots):
        fnm=fname+'_'+mol+'plas.txt'
        f = open(fnm, 'r+')
        header=f.readline()
        f.close()
        head_names=header.split()
        col_num=[head_names.index(x) for x in head_names if coltype in x]
        #print(fname.split('/')[-1],mol,[head_names[col] for col in col_num])
        alldata=np.loadtxt(fnm,skiprows=1)
        time=alldata[:,0]
        data_cols=alldata[:,col_num[0:-1]]
        #very specific kluge because one of the files started much later. 
        if fname.split('-')[-1]=='blockPKA':
            extra=int(30/time[1])
        else:
            extra=0
        pstrt=int(int(tstart)/time[1])
        pend=int(int(tend)/time[1])
        if fnum==0:
            sig_array=np.zeros((len(fname_roots),len(time),len(col_num[0:-1])))
            peaks=np.zeros((len(fname_roots),len(col_num[0:-1])))
        #
        ########################### Signature part ##########################
        #
        if normYN:
            basal=np.mean(data_cols[pstrt+extra:pend+extra],axis=0)
        else:
            basal=np.zeros(len(col_num[0:-1]))
        sig_array[fnum]=data_cols[extra:,:]-basal
        peaks[fnum]=np.max(sig_array[fnum],axis=0)
    all_sig_array[mol]=sig_array
    all_peaks[mol]=peaks
#
sig_ltp=np.zeros((len(fname_roots),len(time),len(col_num[0:-1])))
sig_ltd=np.zeros((len(fname_roots),len(time),len(col_num[0:-1])))
auc=OrderedDict()
#Other possible signature values
#1. time above threshold = number of points in above_thresh array
#2. contiguous time above threshold - how to calculcate this?  A. low pass filter values, B. low pass filter the above thresh array
#3. auc of contiguous time above threshold = np.sum(sig[points extraced in item 2])
#4. 4 time samples of LTP molecules, e.g. 10 s mean (100 points) surrounding 60, 90, 120, 150 s after stim - use in discriminant analysis
#    This would be done prior to peak normalization and summing
time_above=OrderedDict()

#This variant sums molecules after normalizing to peak across paradigms 
for molnum,mol in enumerate(all_molecules):
    peak_norm=np.max(all_peaks[mol],axis=0)
    for fnum,fname in enumerate(fname_roots):
        if mol in ltp_molecules:
            sig_ltp[fnum]=sig_ltp[fnum]+all_sig_array[mol][fnum]/peak_norm
        else:
            sig_ltd[fnum]=sig_ltd[fnum]+all_sig_array[mol][fnum]/peak_norm
#AUC
for fnum,fname in enumerate(fname_roots):
    auc_set=np.zeros((2,2))
    for region in range(len(col_num[0:-1])):
        T_LTP=(float(thresh[1]),float(thresh[0]))[region==0]
        T_LTD=(float(thresh[3]),float(thresh[2]))[region==0]
        above_thresh=[x for x in range(len(sig_ltp[fnum,:,region])) if sig_ltp[fnum,x,region]>T_LTP and x>pend]
        auc_set[0][region]=np.sum(sig_ltp[fnum,above_thresh,region])*time[1]
        above_thresh=[x for x in range(len(sig_ltd[fnum,:,region])) if sig_ltd[fnum,x,region]>T_LTD and x>pend]
        auc_set[1][region]=np.sum(sig_ltp[fnum,above_thresh,region])*time[1]
    auc[fname[fname.find('-'):]]=auc_set
print('######### auc:')
for key,value in auc.items():
    #print('{0:30}  {1:.2f}  {2:.2f}  {3:.2f}  {4:.2f}'.format(key,value.flatten()[0],value.flatten()[1],value.flatten()[2], value.flatten()[3]))
    print('{0:30}  {1:30}'.format(key,np.round(value.flatten(),2)))

figtitle=fnm[0:fnm.find('-')]
auc_label=[[] for p in parval]
domain=[head_names[x].split(coltype)[0] for x in col_num[0:-1]]
for par in range(len(parval)):
    auc_label[par]=[parval[par]+' '+dom[0:6] for dom in domain]
sign_title=args[1]+' vs '+args[2]

#if len(ltd_molecules):
#    pu5.plot_signature(auc_label,sig_ltp,time,figtitle,sign_title,textsize,thresh,sig_ltd)
#else:
#    pu5.plot_signature(auc_label,sig_ltp,time,figtitle,sign_title,textsize,thresh)
