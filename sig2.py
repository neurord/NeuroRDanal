#sig2.py
#evaluate various features to use for signature
#ARGS="subdir/fileroot,LTPmol,LTDmol,tstart tend, t_ltp_dend t_ltp_sp t_ltd_dend t_ltd_sp time_thresh"
#ARGS="Model_SPNspineAChm4R_Gshydr5_GapD,Aphos CKpCam Epac1 Pkc,2ag,10 20,1.8 1.5 0.4 0.2 11"
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

coltype=[1,2,3,4,5,6]#'mean'
trials='3' #make trials=1 if coltype=='mean'
normYN=0
textsize=8
plotYN=0
print_measures=0

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
    if len(thresh)>4:
        time_thresh=int(thresh[4])
    else:
        time_thresh=10 #units are sec.
else:
    thresh=['0', '0', '0', '0']
    time_thresh=10 #units are sec.

all_peaks={}
all_sig_array={}
for molnum,mol in enumerate(all_molecules):
    for fnum,fname in enumerate(fname_roots):
        fnm=fname+'_'+mol+'plas.txt'
        f = open(fnm, 'r+')
        header=f.readline()
        f.close()
        head_names=header.split()
        if type(coltype) is str:
            col_num=[head_names.index(x) for x in head_names if coltype in x][0:-1]
        else:
            col_num=coltype
        alldata=np.loadtxt(fnm,skiprows=1)
        time=alldata[:,0]
        data_cols=alldata[:,col_num]
        #very specific kluge because one of the files started much later. 
        if fname.split('-')[-1]=='blockPKA':
            extra=int(30/time[1])
        else:
            extra=0
        pstrt=int(int(tstart)/time[1])
        pend=int(int(tend)/time[1])
        if fnum==0:
            sig_array=np.zeros((len(fname_roots),len(time),len(col_num)))
            peaks=np.zeros((len(fname_roots),len(col_num)))
        #
        ########################### Signature part ##########################
        #
        if normYN:
            basal=np.mean(data_cols[pstrt+extra:pend+extra],axis=0)
        else:
            basal=np.zeros(len(col_num))
        sig_array[fnum]=data_cols[extra:,:]-basal
        peaks[fnum]=np.max(sig_array[fnum],axis=0)
    all_sig_array[mol]=sig_array
    all_peaks[mol]=peaks
#
# 4 time samples of LTP molecules: mean over a window (e.g. 10s or 100 points) surrounding 60, 90, 120, 150 s after stim - use in discriminant analysis
# Include baseline to be able to use amount above or ratio above baseline. 
win=int(time_thresh/dt[0]/2)
t1=int(60/dt[0])+pstrt-win
t2=int(90/dt[0])+pstrt-win
t3=int(120/dt[0])+pstrt-win
t4=int(150/dt[0])+pstrt-win
sampletimes=[(pstrt,pend),(t1,t1+2*win), (t2,t2+2*win), (t3,t3+2*win), (t4,t4+2*win)]
header='filename '
for fnum,fname in enumerate(fname_roots):
    #loop over trials. create new output line for each trial.  figure out what columns to concatenate for header and to average over for time samples
    samples=[]
    for mol in all_sig_array.keys():
        for t,timepoint in enumerate(sampletimes):
            if fnum==0:
                for col in range(len(col_num)):
                    header=header+str(mol)+'_'+str(int((timepoint[0]+win)*dt[0]))+'c'+str(col)+' '
            samples.append(list(np.mean(all_sig_array[mol][fnum][timepoint[0]:timepoint[1]],axis=0)))
    if fnum==0:
        outfname=fname_roots[0].split('-')[0]+'time_samples.txt'
        f=open(outfname, 'w')
        f.write(header+'\n')
    outputrow=[np.round(val,2) for sublist in samples for val in sublist]
    out=" ".join(str(e) for e in outputrow)
    f.write(fname+' '+out+'\n')
    print(fname,out)
f.close()

sig_ltp=np.zeros((len(fname_roots),len(time),len(col_num)))
sig_ltd=np.zeros((len(fname_roots),len(time),len(col_num)))
auc=OrderedDict()
auc_contig=OrderedDict()
time_above=OrderedDict()
contig_time_above=OrderedDict()
#This variant sums molecules after normalizing to peak across paradigms 
for molnum,mol in enumerate(all_molecules):
    peak_norm=np.max(all_peaks[mol],axis=0)
    for fnum,fname in enumerate(fname_roots):
        if mol in ltp_molecules:
            sig_ltp[fnum]=sig_ltp[fnum]+all_sig_array[mol][fnum]/peak_norm
        else:
            sig_ltd[fnum]=sig_ltd[fnum]+all_sig_array[mol][fnum]/peak_norm

newsig=np.zeros((2,len(fname_roots),len(time),len(col_num)))
#various measures
dur_thresh=int(time_thresh/dt[0])
for fnum,fname in enumerate(fname_roots):
    auc_set=np.zeros((2,2))
    auc_contig_set=np.zeros((2,2))
    time_above_set=np.zeros((2,2))
    contig_time_above_set=np.zeros((2,2))
    for tnum,sig in enumerate([sig_ltp,sig_ltd]):
      for region in range(len(col_num)):
        reg_thresh=(float(thresh[2*tnum+1]),float(thresh[2*tnum]))[region==0]
        #1. auc above threshold, #2. time above threshold
        above_thresh=[x for x in range(len(sig[fnum,:,region])) if sig[fnum,x,region]>reg_thresh and x>pend]
        time_above_set[tnum][region]=len(above_thresh)
        auc_set[tnum][region]=np.sum(sig[fnum,above_thresh,region])*time[1]
        #3. contiguous time above threshold, #4, auc for contiguous time above threshold
        contig_above=h5utils.rolling(above_thresh,dur_thresh)
        contig_time_above_set[tnum][region]=len(contig_above)
        auc_contig_set[tnum][region]=np.sum(sig[fnum,contig_above,region])*time[1]
        if len(contig_above):
            newsig[tnum,fnum,contig_above,region]=sig[fnum,contig_above,region]
    auc[fname[fname.find('-'):]]=auc_set
    time_above[fname[fname.find('-'):]]=time_above_set
    contig_time_above[fname[fname.find('-'):]]=contig_time_above_set
    auc_contig[fname[fname.find('-'):]]=auc_contig_set

#print and display measures
figtitle=fnm[0:fnm.find('-')]
auc_label=[[] for p in parval]
domain=[head_names[x].split(mol)[0] for x in col_num]
for par in range(len(parval)):
    auc_label[par]=[parval[par]+' '+dom[0:6] for dom in domain]
sign_title=args[1]+' vs '+args[2]

if plotYN:
  if len(ltd_molecules):
    pu5.plot_signature(auc_label,sig_ltp,time,figtitle,sign_title,textsize,thresh,sig_ltd)
    pu5.plot_signature(auc_label,newsig[0],time,figtitle,sign_title,textsize,thresh,newsig[1])
  else:
    pu5.plot_signature(auc_label,sig_ltp,time,figtitle,sign_title,textsize,thresh)
    pu5.plot_signature(auc_label,newsig[0],time,figtitle,sign_title,textsize,thresh)

if print_measures:
  print('######### auc using :', thresh, "duration threshold", time_thresh)
  for measure,measure_name in zip([auc,time_above,contig_time_above,auc_contig],['auc','time_above','contig_time_above','auc_contig']):
    print ("########## measure type", measure_name,"###############")
    for key,value in measure.items():
        #print('{0:30}  {1:.2f}  {2:.2f}  {3:.2f}  {4:.2f}'.format(key,value.flatten()[0],value.flatten()[1],value.flatten()[2], value.flatten()[3]))
        print('{0:30}  {1:30} {2:.2f} {3:.2f}'.format(key,np.round(value.flatten(),2),value[0,0]-value[1,0], value[0,1]-value[1,1]))

        
    
