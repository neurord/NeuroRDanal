#sig.py
#calculate LTP/LTD signature from two sets of molecules, separately for spines and dendrites
#USAGE: in python, type ARGS="subdir/fileroot,par1 par2,LTPmol1 LTPmol2,LTDmol1 LTdmol2,basal_start basal_end, T_LTPd T_LTPsp T_LTDd T_LTDsp,sum_name1 sumname2"
#then execfile('sig.py')
#DO NOT PUT ANY SPACES NEXT TO THE COMMAS, DO NOT USE TABS, rows is optional
#LTPmol1 LTPmol2, etc are the names of molecles which produce LTP is sufficiently high (and hinder LTD)
#LTDmol1 LTDmol2, etc are the names of molecles which produce LTD is sufficiently high (and hinder LTP)
#sum_name1 and sum_name2 are prefixes for filenames with output molecule sum traces
#par1 and optionally par2 are specifications of parameter variations, as follows:
#The filenames to read in are constructed as "subdir/fileroot"+"-"+par1+"*"-"+par2+"*"
#DO NOT use hyphens in filenames except for preceding parameter name
#if no parameters specified, then fileroot needs to be full filename (excluding the extension)
#from outside python, type python sig.py "subdir/fileroot [par1 par2] [mol1 mol2]"
#Does not work with only a single file

from __future__ import print_function
from __future__ import division

import os
import numpy as np
from matplotlib import pyplot
import sys  
import glob
from NeuroRDanal import h5utils
from NeuroRDanal import plot_h5 as pu5
import h5py as h5

#######################################################
spinehead="head" #name of spine head from morphology file
dend="dend"
textsize=14 #make bigger for presentations
outputavg=1  #if 1 and args[6] is given will calculate molecule sum, if 2 and args[6] will create output file with molecule sum
window_size=1  #number of seconds on either side of peak value to average for maximum
norm=0 #set to 0 to eliminate baseline subtraction (good for output signatures), set to 1 for baseline subtraction (good for AUC calculation)
trialstats=1
spatialaverage=0
bins=10
trial_auc_ratio=0 #calculate ratio of auc to mean auc for dhpg=0.  Only relevant for Uchi simulations
endpt=2000 #1200 for AUC for Uchi sims, 2000 to find proper peak with bathDaCa; make this another parameter
LTPbas =665.53#567.15 #use to calculate ratio of peak versus basal, e.g. for Dp34 or Dp75
LTDbas=10347.56#10377.95
#######################################################
Avogadro=6.023e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15
msec_per_sec=1000

try:
    args = ARGS.split(",")
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

ftuples,parlist,params=h5utils.argparse(args)
figtitle=args[0].split('/')[-1]

ltp_molecules=args[2].split()
ltd_molecules=args[3].split()
if len(args[5]):
    thresh=args[5].split()
else:
    thresh=['0', '0', '0', '0']

try:
    data.close()
except Exception:
    pass

###################################################
parval=[]
numfiles=len(ftuples)
signature_array=[]
overall_baseline=[]
LTD_auc_all={}
for fnum,ftuple in enumerate(sorted(ftuples, key=lambda x:x[1])):
    data,maxvols,TotVol,trials,seeds,arraysize,p=h5utils.initialize(ftuple,numfiles,parval)
    if len(p):
        params=p[0]
        parval=p[1]
        parlist=[2]
    sig_data=[]
    #
    ##########################################################
    #   Extract region and structure voxels and volumes
    ##########################################################
    if maxvols>1 and fnum==0:
        molecules=data['model']['species'][:]
        structType=data['model']['grid'][:]['type']
        region_list,region_dict,region_struct_dict=h5utils.subvol_list(structType,data['model'])
        #
        try:
            head_index=list(region_dict.keys()).index(spinehead)
        except ValueError:
            head_index=-1
        if head_index>0:
            #create "group" dictionary-spinelist & spinevox, with one more for nonspine voxels
            spinelist,spinevox=h5utils.multi_spines(data['model'])
        else:
            spinelist=''
        if spatialaverage:
            spatial_dict=h5utils.spatial_average(bins,dend,data['model']['grid'])
            vox=[x['vox'] for x in spatial_dict.values()]
            if any(v==0 for v in vox):
                print ("**********Too many bins for spatial average****************")
    #
    ##########################################################
    #   Initialize some arrays and get molecule-region output set information
    ##########################################################
    if fnum==0:
    #
        #Get list of molecules for LTP and list for LTD.  identify which output sets and voxels they are in
        num_ltpmols=len(ltp_molecules)
        num_ltdmols=len(ltd_molecules)
        all_molecules=ltp_molecules+ltd_molecules
        num_mols=len(all_molecules)
        if numfiles>1:
            signature_array=[[] for mol in range(num_mols)]
        if trial_auc_ratio:
            for item in parlist[0]:
                LTD_auc_all[item]={}
    ######################################
    #Calculate region averages, such as indivivdual spines and non-spines
    ######################################
    #
    # do this for each file since they may have different number of samples, or different locations
    out_location,dt,rows=h5utils.get_mol_info(data,all_molecules,maxvols)
    #
    #Which "rows" should be used for baseline value, specifed in args[4]
    sstart,ssend=h5utils.sstart_end(all_molecules, args, 4, out_location,dt,rows)
    molecule_name_issue=0
    if maxvols>1:
        if spatialaverage:
            num_regions=len(spatial_dict.keys())
        elif len(spinelist):
            num_regions=len(spinelist)
        else:
            num_regions=1
        if outputavg:
            LTP_sum=np.zeros((len(trials),np.max(rows),num_regions))
            LTP_sumTot=np.zeros((len(trials),np.max(rows)))
            LTD_sum=np.zeros((len(trials),np.max(rows),num_regions))
            LTD_sumTot=np.zeros((len(trials),np.max(rows)))
        for imol,molecule in enumerate(all_molecules):
          if out_location[molecule]!=-1:
            molecule_pop,time=h5utils.get_mol_pop(data,out_location[molecule],maxvols,trials)
            OverallMean=np.sum(molecule_pop[:,:,:],axis=2)/(TotVol*mol_per_nM_u3)
            #calculate non-spine mean and individual spine means
            spineheader,spinemeans,spineMeanStd=h5utils.region_means_dict(molecule_pop,spinevox,time,molecule,trials)
            if spatialaverage:
                spacehead,spaceMeans,spaceMeanStd=h5utils.region_means_dict(molecule_pop,spatial_dict,time,molecule,trials)
                summeans=spaceMeans
                outheader=spacehead
            else:
                summeans=spinemeans
                outheader=spineheader
            if outputavg:
                if molecule in ltp_molecules:
                    LTP_sum=LTP_sum+summeans
                    LTP_sumTot=LTP_sumTot+OverallMean
                if molecule in ltd_molecules:
                    LTD_sum=LTD_sum+summeans
                    LTD_sumTot=LTD_sumTot+OverallMean
            if numfiles>1:
                #dimensions will be number of molecules x sample times x (1+ num spines)
                sig_data.append(np.mean(summeans,axis=0))
            else:
                #dimensions will be number of molecules x number of trials x sample times x 1+ num spines
                sig_data.append(summeans)
          else:
              if fnum==0 and molecule_name_issue==0:
                  print("Choose molecules from:", data['model']['species'][:])
                  molecule_name_issue=1
              sig_data.append(np.zeros(len(time)))
    ######################################
    #minimal processing needed if only a single voxel.
    ######################################
    else:
        voxel=0
        if outputavg:
            LTP_sum=np.zeros((len(trials),np.max(rows)))
            LTD_sum=np.zeros((len(trials),np.max(rows)))
            outheader="all"
        for mol in all_molecules:
          if out_location[mol]!=-1:
            outset = out_location[mol]['location'].keys()[0]
            imol=out_location[mol]['location'][outset]['mol_index']
            tempConc=np.zeros((len(trials),out_location[mol]['samples']))
            time=data[trials[0]]['output'][outset]['times'][:]/msec_per_sec
            #generate output files for these cases
            for trialnum,trial in enumerate(trials):
                tempConc[trialnum]=data[trial]['output'][outset]['population'][:,voxel,imol]/TotVol/mol_per_nM_u3
            if outputavg:
                if molecule in ltp_molecules:
                    LTP_sum=LTP_sum+tempConc
                if molecule in ltd_molecules:
                    LTD_sum=LTD_sum+tempConc
            if numfiles>1:
                 sig_data.append(np.mean(tempConc,axis=0))
                 #sig_data dimensions=number of molecules x sample times
            else:
                #sig_data dimensions=number of molecules (x number of trials) x sample times
                sig_data.append(tempConc)
          else:
              if fnum==0 and molecule_name_issue==0:
                  print("Choose molecules from:", all_molecules)
                  molecule_name_issue=1
              sig_data.append(np.zeros(len(time)))
              print("molecule", molecule, "not found in output data!!!!!!!!!!!")
    ######################################
    #Whether 1 voxel or multi-voxel, create array of means for all molecules, all files, all trials
    ##########################################
    if numfiles>1:
        for mol in range(num_mols):
            signature_array[mol].append(sig_data[mol])
    else:
        signature_array=sig_data
    #############################################
    # output of sum of molecules - useful for showing signature
    #############################################################
    if outputavg:
        if spatialaverage:
            pu5.plot3D(LTP_sum,trials,ftuple[0][0:-3],ltp_molecules,spatial_dict.keys(),time)
            if num_ltdmols:
                pu5.plot3D(LTD_sum,trials,ftuple[0][0:-3],ltd_molecules,spatial_dict.keys(),time)
        sum_array=[LTP_sum,LTD_sum]
        Tot_array=[LTP_sumTot,LTD_sumTot]
        if len(args)>6:
            sum_name=args[6].split()
        if len(args)<=6 or len(sum_name)==0:
            sum_name=[ltp_molecules[0],ltd_molecules[0]]
            outputavg=1
        if len(sum_name)<2:
            sum_name.append('ltdmol')
        if trialstats:
            basalstrt=sstart[0] 
            basalend=ssend[0] 
            if fnum==0:
                print('STATISTICS', sum_name[0],'trials, stderr  ',sum_name[1],'trials, stderr')
            LTP_basal=np.mean(LTP_sumTot[:,basalstrt:basalend],axis=1)
            LTD_basal=np.mean(LTD_sumTot[:,basalstrt:basalend],axis=1)
            overall_baseline.append(LTD_basal)
            print('basal',np.round(LTP_basal,1), np.round(np.std(LTP_basal)/np.sqrt(len(trials)),1),
                  np.round(LTD_basal,0), np.round(np.std(LTD_basal)/np.sqrt(len(trials)),1))
            LTP_peak=[LTP_sumTot[i,ssend[0]:endpt].max() for i in range(len(trials))]
            LTD_peak=[LTD_sumTot[i,ssend[0]:endpt].max() for i in range(len(trials))]
            LTP_min=[LTP_sumTot[i,ssend[0]:endpt].min() for i in range(len(trials))]
            LTD_min=[LTD_sumTot[i,ssend[0]:endpt].min() for i in range(len(trials))]
            normD=np.zeros((len(trials),np.shape(LTD_sumTot)[1]))
            normP=np.zeros((len(trials),np.shape(LTD_sumTot)[1]))
            for i in range(len(trials)):
                normD[i,:]=LTD_sumTot[i,:]-LTD_basal[i]
                normP[i,:]=LTP_sumTot[i,:]-LTP_basal[i]
            LTD_auc=[np.sum(normD[i][ssend[0]:endpt])*dt[0]/msec_per_sec for i in range(len(trials))]
            LTP_auc=[np.sum(normP[i][ssend[0]:endpt])*dt[0]/msec_per_sec for i in range(len(trials))]
            if trial_auc_ratio:
                LTD_auc_all[ftuple[1][0]][ftuple[1][1]]=LTD_auc
            print('peak',np.round(LTP_peak,1), np.round(LTD_peak,1), 'min',np.round(LTP_min,1), np.round(LTD_min,1))
            #print('peak ratio',np.round(np.array(LTP_peak)/LTPbas,3), np.round(np.array(LTD_peak)/LTDbas,3), 'min',np.round(np.array(LTP_min)/LTPbas,3), np.round(np.array(LTD_min)/LTDbas,3))
            print('auc', np.round(LTP_auc,2),np.round(np.mean(LTP_auc),2),np.round(LTD_auc,2),np.round(np.mean(LTD_auc),2))
        for num,nm in enumerate(sum_name):
            sum_header='time  '
            for item in outheader.split():
                newitem=[item.split('_')[-1]+nm+'_t'+str(t)+' ' for t in range(len(trials))]
                sum_header=sum_header+''.join(newitem)
            mean_head=[item.split('_')[-1]+nm+'mean ' for item in outheader.split()]
            stdev_head=[item.split('_')[-1]+nm+'stdev ' for item in outheader.split()]
            sum_header=sum_header+"".join(mean_head)+"".join(stdev_head)
            Overall_header=[nm+'_t'+str(t)+' ' for t in range(len(trials))]+[nm+"mean ",nm+"stdev"]
            sum_header=sum_header+"".join(Overall_header)
            num_trials=len(trials)
            cols=num_trials*num_regions
            sum_rows=np.shape(LTP_sum)[1]
            outdata=np.zeros((sum_rows,cols))
            outmean=np.zeros((sum_rows,num_regions))
            outstd=np.zeros((sum_rows,num_regions))
            outOverall=np.zeros((sum_rows,num_trials+2))
            for p in range(num_regions):
                outdata[:,p*num_trials:(p+1)*num_trials]=sum_array[num][:,:,p].T
                #outdata[:,p*num_trials:(p+1)*num_trials]=LTD_sum[:,:,p].T
                outmean[:,p]=np.mean(outdata[:,p*num_trials:(p+1)*num_trials],axis=1)
                outstd[:,p]=np.std(outdata[:,p*num_trials:(p+1)*num_trials],axis=1)
            outOverall[:,range(num_trials)]=Tot_array[num].T
            outOverall[:,num_trials]=np.mean(outOverall[:,range(num_trials)],axis=1)
            outOverall[:,num_trials+1]=np.std(outOverall[:,range(num_trials)],axis=1)
            outputdata=np.column_stack((time,outdata,outmean,outstd,outOverall))
            if outputavg>1:
                if spatialaverage:
                    suffix='dend'
                else:
                    suffix='plas'
                outfname=ftuple[0][0:-3]+'_'+nm+suffix+'.txt'
                f=open(outfname, 'w')
                f.write(sum_header+'\n')
                np.savetxt(f, outputdata, fmt='%.4f', delimiter=' ')
                f.close()
            ########### Calculate baseline, peak, min and ratios
            window=int(window_size/dt[0])
            base_sum=outOverall[sstart[0]:ssend[0],num_trials].mean()
            peakpt_sum=outOverall[ssend[0]:,num_trials].argmax()+ssend[0]
            #peaktime[pnum,imol]=peakpt*dt[imol]
            peak_sum=outOverall[peakpt_sum-window:peakpt_sum+window,num_trials].mean()
            minpt_sum=outOverall[ssend[0]:,num_trials].argmin()+ssend[0]
            min_sum=outOverall[minpt_sum-window:minpt_sum+window,num_trials].mean()
            if fnum==0 and num==0:
                print ('fname                     base   peak    inc     min     dec')
            print ('{0:25}  {1:.1f}  {2:.1f}  {3:.3f}  {4:.1f}  {5:.3f}'.format(''.join(parval[fnum])+nm,base_sum,peak_sum,peak_sum/base_sum,min_sum,min_sum/base_sum))
if trial_auc_ratio:
    print('overall baseline',np.round(np.mean(overall_baseline),3),'use on line 251 for more consistent baseline subtraction' )
    print ('auc calculated between',ssend[0]*dt,'and', endpt*dt, ', ratio with the 0dhpg case')
    print('dur  dhpg  trials                         mean stdev stderr')
    for dur,durdict in sorted(LTD_auc_all.items(), key=lambda x:x[1]):
        for dhpg in sorted(durdict.keys(), key=lambda x:int(x)):
            auc_ratio=durdict[dhpg]/np.mean(durdict['0'])
            print(dur,dhpg,auc_ratio, np.round(np.mean(auc_ratio),3), np.round(np.std(auc_ratio),3), np.round(np.std(auc_ratio)/np.sqrt(len(auc_ratio)),3))
#####################################################################
#Calculate signature
#####################################################################
auc_label=[]
sign_title=''
for mol in ltp_molecules:
    sign_title=sign_title+'+'+mol
for mol in ltd_molecules:
    sign_title=sign_title+'-'+mol
if maxvols==1:
    auc=np.zeros(len(parval))
    num_regions=1
else:
    auc=np.zeros((len(parval),num_regions))
ltp_above_thresh=np.zeros((len(parval),num_regions))
ltd_above_thresh=np.zeros((len(parval),num_regions))
lengths=[np.shape(signature_array[0][x])[0] for x in range(numfiles)]
sig_ltp=np.zeros((len(parval),np.max(lengths),num_regions))
sig_ltd=np.zeros((len(parval),np.max(lengths),num_regions))
#############################
#customize this part.  E.g.
#add values of LTP molecules, subtract LTD molecules; or accumulate each signature separately
def sig_subtract(sig_array,strt,send,num_regions,ltp_samples,normYN):
    if normYN:
        basal=np.mean(sig_array[strt:send],axis=0)
    else:
        basal=0
    sig_subtracted=sig_array-basal
    extra=ltp_samples-np.shape(sig_subtracted)[0]
    if extra:
        extra_zeros=np.zeros((extra,num_regions))
        sig_subtracted=np.vstack((sig_subtracted,extra_zeros))
    return sig_subtracted

for f in range(len(parval)):
    for each_mol in ltp_molecules:
        col=all_molecules.index(each_mol)
        sig_ltp[f]=sig_ltp[f]+sig_subtract(signature_array[col][f],sstart[col],ssend[col],num_regions,np.shape(sig_ltp[f])[0],norm)
    for each_mol in ltd_molecules:
        col=all_molecules.index(each_mol)
        sig_ltd[f]=sig_ltd[f]+sig_subtract(signature_array[col][f],sstart[col],ssend[col],num_regions,np.shape(sig_ltd[f])[0],norm)
#signature dimensions=num files/trials x sample times x (1+numspines)
#End customization
#############################
#area between signature and basal
if maxvols==1:
    for par in range(len(parval)):
        label=h5utils.join_params(parval[par],params)
        auc[par]=np.sum(sig_ltp[par,:])*dt[0]/msec_per_sec
        if len(ltd_molecules):
            auc[par]=auc[par]-np.sum(sig_ltd[par,:])*dt[0]/msec_per_sec
        auc_label.append(label+" auc="+str(np.round(auc[par],2)))
else:
    auc_label=[[] for sp in range(len(parval))]
    if spatialaverage:
        regionnames=spatial_dict.keys()
    else:
        regionnames=spinelist
    for par in range(len(parval)):
        label=h5utils.join_params(parval[par],params)
        #label=parval[par][0]
        for sp in range(num_regions):
            T_LTP=(float(thresh[1]),float(thresh[0]))[sp==0]
            T_LTD=(float(thresh[3]),float(thresh[2]))[sp==0]
            auc[par,sp]=np.sum(sig_ltp[par,:,sp])*dt[0]/msec_per_sec
            ltp_above_thresh[par,sp]=np.sum(sig_ltp[par,sig_ltp[par,:,sp]>T_LTP,sp]-T_LTP)*dt[0]/msec_per_sec
            if len(ltd_molecules):
                auc[par,sp]=auc[par,sp]-np.sum(sig_ltd[par,:,sp])*dt[0]/msec_per_sec
                ltd_above_thresh[par,sp]=np.sum(sig_ltd[par,sig_ltd[par,:,sp]>T_LTD,sp]-T_LTD)*dt[0]/msec_per_sec
            #auc_label[par].append(label+" auc="+str(np.round(auc[par,sp],1))+" "+spinelist[sp])
            auc_label[par].append(label+' '+str(np.round(auc[par,sp],1))+" "+regionnames[sp])
pyplot.ion()
if len(ltd_molecules):
    pu5.plot_signature(auc_label,sig_ltp,dt[0],figtitle,sign_title,textsize,thresh,sig_ltd)
    numcol=2
else:
    numcol=1
    pu5.plot_signature(auc_label,sig_ltp,dt[0],figtitle,sign_title,textsize,thresh)
if spatialaverage:
    pu5.plot3D(sig_ltp,parval,figtitle,ltp_molecules,spatial_dict.keys(),time)
    if len(ltd_molecules):
        pu5.plot3D(sig_ltd,parval,figtitle,ltd_molecules,spatial_dict.keys(),time)
        
print("area above threshold for LTP and LTD using", thresh)
for par in range(len(parval)):
    print('{0:20} {1:8} {2:8}'.format(''.join(parval[par]), np.round(ltp_above_thresh[par],3), np.round(ltd_above_thresh[par],3)))

    
