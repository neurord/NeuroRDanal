#sig2-fig.py
#Generate signature traces and plot them, or create files of time samples for analysis in other software
#ARGS="subdir/fileroot,parameter,LTPmol,LTDmol,tstart tend,4 thresholds (optional)" where tstart and tend used to calculate basal value
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

#specify either 'mean' or column numbers, or region names - must match sig.py output
#These should probably be additional parameters/arguments
#coltype=[1,2,3,4,5,6]
#coltype=[str(x) for x in np.linspace(1.0,19.0,10)] #for spatial samples in long dendrite
coltype=['sa1[0]','sa1[1]','sa1[2]','sa1[3]','sa1[4]','sa1[5]','sa1[6]','sa1[7]','sa1[8]','sa1[9]']
coltype=['nonspine','sa1[0]']
suffix='plas'  #either dend to use spatial average files or plas to use spine outputs from sig.py  
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
    pattern=args[0]+'-'+'*-'.join(args[1].split())
else:
    pattern=args[0]
globpattern=pattern+"*"+suffix+".txt"
fnames = glob.glob(globpattern)
print(globpattern, len(fnames), 'files found')
fname_roots=sorted(set([fnm[0:fnm.rfind('_')] for fnm in fnames]))
fname_endings=set([fnm[fnm.rfind('_'):] for fnm in fnames])
parval=[fnm[fnm.find('-')+1:].replace('-',' ') for fnm in fname_roots]

#assign other inputs (args) to parameters
ltp_molecules=args[2].split()
ltd_molecules=args[3].split()
tstart,tend=args[4].split()

#if provide file suffix, then sample times file will be written, and plots will not be generated
if len(args)>7:
    fname_suffix=args[7]
    plotYN=0
    normYN=0
    if len(args[5]):
        mean_dur=args[5].split()
    else:
        mean_dur=[60,30] #units are sec
    if len(args[6]):
        samp_times=args[6].split()
    else:
        samp_times=[60,120]
else: #if no file suffix, optionally provide threshold for calculating AUC.  Normalize only if specified above
    fname_suffix=''
    plotYN=1
    if len(args)>5:  #optional specification of thresholds
        thresh=args[5].split()
        if len(thresh)==1:
            thresh=[thresh for i in range(4)]
        elif len(thresh)==2:
            thresh=thresh+thresh
    else:
        thresh=[0,0,0,0]
    print("thresholds", thresh)

all_molecules=sorted(ltp_molecules+ltd_molecules)
num_mols=len(all_molecules)

#Read in all files and assemble relevant columns into large data matrix
all_peaks={}
all_sig_array={}
for molnum,mol in enumerate(all_molecules):
    for fnum,fname in enumerate(fname_roots):
        fnm=fname+'_'+mol+suffix+'.txt'
        f = open(fnm, 'r+')
        header=f.readline()
        f.close()
        head_names=header.split()
        temp=[]
        if coltype=='mean':
            col_num=[head_names.index(x) for x in head_names if coltype in x][0:-1]
        elif type(coltype) is list and type(coltype[0]) is str:
            for col in coltype:
                temp.append([head_names.index(x) for x in head_names if '_t' in x and x.startswith(col)])
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
            time=time[:-extra]
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

domain=[head_names[x].split(mol) for x in col_num]
if plotYN:
    #calculate signatures
    sig_ltp,sig_ltd=h5utils.calc_sig(all_sig_array,all_peaks,all_molecules,ltp_molecules,trials,fname_roots,time,col_num)
    #show figures (and create output files for publication figures?)
    figtitle=fnm[0:fnm.find('-')]
    label=[[] for p in parval]
    for par in range(len(parval)):
        label[par]=[parval[par]+' '+dom[0][0:6]+dom[1] for dom in domain]
    sign_title=args[2]+' vs '+args[3]
    #
    pyplot.ion()
    if len(ltd_molecules):
        pu5.plot_signature(label,sig_ltp,time[1],figtitle,sign_title,textsize,thresh,sig_ltd)
    else:
        pu5.plot_signature(label,sig_ltp,time[1],figtitle,sign_title,textsize,thresh)
else:
    #or create output files for analysis of signatures. either 1 file with multiple points
    h5utils.sig_outfile_multi(mean_dur,time[1],samp_times,pstrt,pend,fname_roots,fname_suffix,col_num,trials,all_sig_array)
    #or mulitple files, each with basal and 1 other point
    #h5utils.sig_outfile_single(mean_dur,time,samp_times,pstart,pend,fname_roots,fname_suffix,col_num,trials,all_sig_array)

#To write the signature file traces:
if any([col.startswith('sa1[') for col in coltype]):
    newdomain=[[dom[0].replace('sa1[','sp').replace(']',''),dom[1]] for dom in domain]
    meanlabl=np.unique([dom[0] for dom in newdomain])
    reg=[dom.split('sp')[1] for dom in meanlabl]
    xvalues=sorted([float(x.replace('ine','-1')) for x in reg])
else:
    newdomain=domain
    meanlabl=np.unique([dom[0] for dom in domain])
    xvalues=sorted([float(dom) for dom in meanlabl])
if np.min(xvalues)<0:
    xvalues=[x-np.min(xvalues) for x in xvalues]

if plotYN:
    #set up arrays for storing mean over trials
    num_regions=int(len(col_num)/trials)
    mean_sigs_ltp=np.zeros((len(fname_roots),len(time),num_regions))
    mean_sigs_ltd=np.zeros((len(fname_roots),len(time),num_regions))
    #abbreviated column labels
    newlabel=[''.join(newdomain[x]) for x in range(len(newdomain))]
    sig=sig_ltp-sig_ltd
    if suffix=='plas':
        new_suffix='sp'
    elif suffix=='dend':
        new_suffix='dnd'
    else:
        new_suffix='unk'
    for ii,file_cond in enumerate(fname_roots):
        for array,plastype,molecules in zip([sig_ltp,sig_ltd,sig],['ltp','ltd','sig'],[ltp_molecules,ltd_molecules,ltp_molecules+ltd_molecules]):
            #Automatic file naming and column headers, Write the file
            fname=file_cond+'_'+plastype+''.join(molecules)+new_suffix+'.txt'
            param_name='_'.join(file_cond.split('-')[1:]).replace('stim','')
            header='time '+' '.join([plastype+'_'+param_name+labl for labl in newlabel])
            f=open(fname, 'w')
            f.write(' '.join(molecules)+'\n')
            f.write(header+'\n')
            np.savetxt(f,np.column_stack((time,array[ii])), fmt='%.4f', delimiter=' ')
            f.close()
            #average over trials, write to file
            mean_sig=np.zeros(np.shape(mean_sigs_ltp)[1:])
            for reg in range(num_regions):
                colset=[tr+reg*trials for tr in range(trials)]
                mean_sig[:,reg]=np.mean(array[ii,:,colset],axis=0)
            if plastype=='ltp':
                mean_sigs_ltp[ii,:,:]=mean_sig
            elif plastype=='ltd':
                mean_sigs_ltd[ii,:,:]=mean_sig
            fname=file_cond+'_'+plastype+''.join(molecules)+new_suffix+'avg.txt'
            header='time '+' '.join([plastype+'_'+param_name+labl for labl in meanlabl])
            f=open(fname, 'w')
            f.write(' '.join(molecules)+'\n')
            f.write(header+'\n')
            np.savetxt(f,np.column_stack((time,mean_sig)), fmt='%.4f', delimiter=' ')
            f.close()

    for array,mean_array,mols in zip([sig_ltp,sig_ltd],[mean_sigs_ltp,mean_sigs_ltd],[ltp_molecules,ltd_molecules]):
        #3D plot for long morphology
        pu5.plot3D(mean_array,parval,figtitle,mols,xvalues,time)
        #pu5.plot3D(array,parval,figtitle,mols,xvalues,time)
    pu5.plot3D(mean_sigs_ltp-mean_sigs_ltd,parval,figtitle,ltp_molecules+ltd_molecules,xvalues,time)
