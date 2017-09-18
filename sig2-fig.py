#sig2-fig.py
#Generate signature traces and plot them, or create files of time samples for analysis in other software
#ARGS="subdir/fileroot,parameter,LTPmol,LTDmol,tstart tend,4 thresholds,sample times" where tstart and tend used to calculate basal value
#ARGS="Model_SPNspineAChm4R_Gshydr5_GapD,,Aphos CKpCam Epac1 Pkc,2ag,10 20"
#for file output:
#ARGS="subdir/fileroot,parameter,LTPmol,LTDmol,tstart tend,sample durs,sample times,file_suffix"
#ARGS="Model_SPNspineAChm4R_Gshydr5_GapD,,Aphos CKpCam Epac1 Pkc,2ag,10 20,30 60,60 120,ACEP2"
#last 3 parameters and 2nd parameter is optional
#############TODO: add in the other signatures, such as AUC and contiguous time above from sig2.py?  OR, add that as option?  Then replace sig2.py


from __future__ import print_function
from __future__ import division

import os
import numpy as np
from matplotlib import pyplot
import sys  
import glob
from NeuroRDanal import h5utils
from NeuroRDanal import plot_h5 as pu5

#specify either 'mean' or column numbers, or region names
#These should probably be additional parameters/arguments
#coltype=[1,2,3,4,5,6]
coltype=['nonspine','sa1[0]']#,'sa1[1]','sa1[2]','sa1[3]','sa1[4]','sa1[5]','sa1[6]','sa1[7]','sa1[8]','sa1[9]']
trials=3 #make trials=1 if coltype=='mean'
normYN=1
textsize=8
print_measures=1

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

#if provide file suffix, then files will be written, and plots will not be generated
if len(args[7]):
    fname_suffix=args[7]
    plotYN=0
    normYN=0
else:
    fname_suffix=''
    plotYN=1

if len(args[5]):
    if plotYN:
        thresh=args[5]
    else:
        time_thresh=args[5].split()
else:
    if plotYN:
        thresh=[0,0,0,0]
    else:
        time_thresh=[60,30] #units are sec.
if len(args[6]):
    samp_times=args[6].split()
else:
    samp_times=[60,120]

all_molecules=sorted(ltp_molecules+ltd_molecules)
num_mols=len(all_molecules)

#Read in all files and assemble relevant columns into large data matrix
all_peaks={}
all_sig_array={}
for molnum,mol in enumerate(all_molecules):
    for fnum,fname in enumerate(fname_roots):
        fnm=fname+'_'+mol+'plas.txt'
        f = open(fnm, 'r+')
        header=f.readline()
        f.close()
        head_names=header.split()
        temp=[]
        if coltype=='mean':
            col_num=[head_names.index(x) for x in head_names if coltype in x][0:-1]
        elif type(coltype) is list and type(coltype[0]) is str:
            for col in coltype:
                temp.append([head_names.index(x) for x in head_names if '_t' in x and col in x])
            col_num=[val for sublist in temp for val in sublist]
        else:
            col_num=coltype
        if molnum==0:
            print ('fname:', fname, 'columns:', col_num)
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
        ########################### normalize (optional) and create arrays of molecules ##########################
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

if plotYN:
    #calculate signatures
    sig_ltp,sig_ltd=h5utils.calc_sig(all_sig_array,all_peaks,all_molecules,ltp_molecules,trials,fname_roots,time,col_num)
    #show figures (and create output files for publication figures?)
    figtitle=fnm[0:fnm.find('-')]
    label=[[] for p in parval]
    domain=[head_names[x].split(mol) for x in col_num]
    for par in range(len(parval)):
        label[par]=[parval[par]+' '+dom[0][0:6]+dom[1] for dom in domain]
    sign_title=args[2]+' vs '+args[3]
    #
    pyplot.ion()
    if len(ltd_molecules):
        pu5.plot_signature(label,sig_ltp,time,figtitle,sign_title,textsize,thresh,sig_ltd)
    else:
        pu5.plot_signature(label,sig_ltp,time,figtitle,sign_title,textsize,thresh)
else:
    #or create output files for analysis of signatures
    h5utils.sig_outfile_multi(time_thresh,time,samp_times,pstrt,pend,fname_roots,fname_suffix,col_num,trials,all_sig_array)
    #
    #h5utils.sig_outfile_single(time_thresh,time,samp_times,pstart,pend,fname_roots,fname_suffix,col_num,trials,all_sig_array)