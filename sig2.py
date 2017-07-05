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
import pandas as pd

coltype='mean'
normYN=0
textsize=8
time_thresh=10 #units are sec.  Make this parameter later

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
auc_contig=OrderedDict()
time_above=OrderedDict()
contig_time_above=OrderedDict()
# 4 time samples of LTP molecules, e.g. 10 s mean (100 points) surrounding 60, 90, 120, 150 s after stim - use in discriminant analysis
#    This would be done prior to peak normalization and summing

#This variant sums molecules after normalizing to peak across paradigms 
for molnum,mol in enumerate(all_molecules):
    peak_norm=np.max(all_peaks[mol],axis=0)
    for fnum,fname in enumerate(fname_roots):
        if mol in ltp_molecules:
            sig_ltp[fnum]=sig_ltp[fnum]+all_sig_array[mol][fnum]/peak_norm
        else:
            sig_ltd[fnum]=sig_ltd[fnum]+all_sig_array[mol][fnum]/peak_norm
#various measures
for fnum,fname in enumerate(fname_roots):
    auc_set=np.zeros((2,2))
    auc_contig_set=np.zeros((2,2))
    time_above_set=np.zeros((2,2))
    contig_time_above_set=np.zeros((2,2))
    for tnum,sig in enumerate([sig_ltp,sig_ltd]):
      for region in range(len(col_num[0:-1])):
        reg_thresh=(float(thresh[2*tnum+1]),float(thresh[2*tnum]))[region==0]
        dur_thresh=int(time_thresh/dt[0])
        #1. auc above threshold, #2. time above threshold
        above_thresh=[x for x in range(len(sig[fnum,:,region])) if sig[fnum,x,region]>reg_thresh and x>pend]
        time_above_set[tnum][region]=len(above_thresh)
        auc_set[tnum][region]=np.sum(sig_ltp[fnum,above_thresh,region])*time[1]
        #3. contiguous time above threshold, #4, auc for contiguous time above threshold
        df_sig = pd.DataFrame(sig[fnum, :, region])
        above_thresh=(df_sig[pend:] > reg_thresh).rolling(dur_thresh).sum() >= dur_thresh
        contig_time_above_set[tnum][region]=above_thresh.any()
        above=[x+pend for x in above_thresh if x == True]
        auc_contig_set[tnum][region]=sum(sig[fnum,above,region])
    auc[fname[fname.find('-'):]]=auc_set
    time_above[fname[fname.find('-'):]]=time_above_set
    contig_time_above[fname[fname.find('-'):]]=contig_time_above_set
    auc_contig[fname[fname.find('-'):]]=auc_contig_set
print('######### auc:')
for measure,measure_name in zip([auc,time_above,contig_time_above,auc_contig],['auc','time_above','contig_time_above','auc_contig']):
    print ("########## measure type", measure_name,"###############")
    for key,value in measure.items():
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
