# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 10:58:05 2020

@author: kblackw1
in python window, type ARGS="subdir/fileroot,par1 par2,mol1 mol2,sstart ssend,rows" then execfile('path/to/file/nrdh5_anal.py')
#from outside python, type python nrdh5_analv2,py "subdir/fileroot [par1 par2] [mol1 mol2] [sstart ssend]"
#DO NOT PUT ANY SPACES NEXT TO THE COMMAS, DO NOT USE TABS, rows is optional
#mol1 mol2, etc are the names of molecles to process
#par1 and optionally par2 are specifications of parameter variations, as follows:
#The filenames to read in are constructed as "subdir/fileroot"+"-"+par1+"*"-"+par2+"*"
#DO NOT use hyphens in filenames except for preceding parameter name
#if no parameters specified, then fileroot needs to be full filename (excluding the extension)
#e.g. ARGS="../Repo/plc/Model_PLCassay,Ca GaqGTP,Ca GaqGTP Ip3,15 20" time units are sec
#e.g. ARGS="plc/Model_PLCassay_Ca1,Ca Gaq,GTP IP3"
#if mol ommitted, then all molecules processed.  if sstart ssend are ommitted, then calculates basal from 7.5-10% of runtime
additional parameters lines 27-45
"""
import numpy as np
import sys

import plot_h5V2 as pu5
from nrd_output import nrdh5_output
from nrd_group import nrdh5_group

#probably add some of these to args with defaults once real arg parser is written
submembname='sub'
dendname="dend"
spinehead="head"
stimspine=['sa1[0]'] #list of stimulated spines
spatial_bins=1  #number of spatial bins to subdivide dendrite to look at spatial gradients
window_size=0.1  #number of msec on either side of peak value to average for maximum
num_LTP_stim=1 #number of 100Hz trains - used to determine when stimulation is over and to search for molecule decay
#These control what output is printed or written
show_inject=0
write_output=0 #one file per molecules per input file
output_auc=0 #one file per molecule per set of input files
showplot=1 #0 for none, 1 for overall average, 2 for spine concentration
show_mol_totals=0
print_head_stats=0
textsize=14
feature_list=[]#'auc','amplitude']
#these molecules MUST be specified as plot_molecules
mol_pairs=[]#[['CKpCamCa4','ppERK']]#,['ppERK','pSynGap']]
pairs_timeframe=[]#[200,2000] #units are sec
basestart_time=0#2200 #make this value 0 to calculate AUC using baseline calculated from initial time period
write_trials=0

#for totaling subspecies, the default outputset is _main_, 
#can specify outset= something (line 84) to use different outputset

sub_species={'ras':['rasGap','RasGTPGap'], 
             'rap':['rap1Gap', 'Rap1GTPGap'],
             'Ras': ['pShcGrb2SosRas', 'CamCa4GRFRas', 'Raf1Ras', 'dRaf1Ras','dRaf1RasMEK', 'dRaf1RaspMEK','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK', 'bRafRas', 'bRafRasMEK','bRafRaspMEK', 'RasGTP', 'RasGDP', 'RasSynGap', 'RasGTPGap', 'RaspSynGap'],
             'Rap1GTP':['bRafRap1MEK', 'bRafRap1pMEK', 'bRafRap1', 'Raf1Rap1', 'Rap1GTP','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK'],
             'PKA':['PKA', 'PKAcAMP2', 'PKAcAMP4', 'PKAr'], 
             'erk':['ppERK','pERK'], 
             'RasGTP':['Raf1Ras', 'dRaf1Ras', 'dRaf1RasMEK', 'dRaf1RaspMEK', 'bRafRas', 'bRafRasMEK', 'bRafRaspMEK', 'RasGTP','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK'], 
             'RasSyn':['RasSynGap', 'RaspSynGap'], 
             'Rap1Syn':['Rap1SynGap', 'Rap1pSynGap'],
             'cAMP': ['EpacAMP', 'cAMP','PDE4cAMP','PDE2cAMP', 'PDE2cAMP2', 'PKAcAMP2', 'PKAcAMP4'], 
             'ERK':['pERK', 'ppERK', 'pERKMKP1', 'ppERKMKP1', 'ppMEKERK', 'ppMEKpERK', 'ppERKpShcGrb2Sos'], 
             'free_Syn':['SynGap', 'pSynGap'],
             'CamCa':['CamCa4','CamCa2C','CamCa2N','PP2BCamCa2C','PP2BCamCa2N', 'PP2BCamCa4']}

#subspecies of tot_species do NOT need to be specified as plot_molecules 
sub_species={'syn_act':['RasSynGap', 'RaspSynGap','Rap1SynGap', 'Rap1pSynGap','SynGap','pSynGap'],'syn_inact':['pSynden','CKpCamCa4SynGap'],'Rasgrf':['CamCa4GRFRas','CamCa4GRF'],'grf':['CamCa4GRF','CamCa4GRFRas'],'gef':['RasGTP','Rap1GTP']}
#weights are multiplies for calculating weighted sum, i.e., signature
#weight={'RasSynGap':-1, 'RaspSynGap':-1.4,'Rap1SynGap':-1, 'Rap1pSynGap':-1.9,'SynGap':-1,'pSynGap':-1.65,'pSynden':1,'CKpCamCa4SynGap':1}
weight={}  #setting weight to empty dictionary will make all weights = 1
tot_species=[]#['syn_act','syn_inact','grf','gef']   
############## END OF PARAMETERS #################
try:
    args = ARGS.split(",")
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

if len(args[2].split()):
    plot_molecules=args[2].split()
else:
    plot_molecules=None
    
figtitle=args[0].split('/')[-1]+args[1]

og=nrdh5_group(args,tot_species)
for fnum,ftuple in enumerate(og.ftuples):
    data=nrdh5_output(ftuple)
    data.rows_and_columns(plot_molecules,args)
    data.molecule_population()
    print(data.data['model']['grid'][:])
    if data.maxvols>1:
        data.region_structures(dendname,submembname,spinehead,stimspine) #stimspine is optional
        if spatial_bins>0:
            data.spatial_structures(spatial_bins,dendname)
        data.average_over_voxels()
    # need to add another total array for different regions (to use for signature)
    #Default outputset is _main_, can specify outset= something
    data.total_subspecies(tot_species,sub_species,args,weights=weight)
    og.conc_arrays(data)
    #Now, print or write some optional outputs
    if write_output:
        data.write_average()
    if 'event_statistics' in data.data['trial0']['output'].keys() and show_inject:
        print ("seeds", data.seeds," injection stats:")
        print('molecule             '+'    '.join(data.trials))
        for imol,inject_sp in enumerate(data.data['model']['event_statistics'][:]):
            inject_num=np.column_stack([data.data[t]['output']['event_statistics'][0] for t in data.trials])
            print (inject_sp.split()[-1].rjust(20),inject_num[imol])
    if print_head_stats:
        data.print_head_stats()
#extract some features from the group of data files
#EXTRACT FEATURES OF total array to add in sig.py functionality
#Default numstim = 1, so that parameter not needed for single pulse
#other parameter defaults:  lo_thresh_factor=0.2,hi_thresh_factor=0.8, std_factor=1
#another parameter default: end_baseline_start=0 (uses initial baseline to calculate auc).
#Specify specific sim time near end of sim if initialization not sufficient for a good baseline for auc calculation
og.trace_features(data.trials,window_size,std_factor=1,numstim=num_LTP_stim,end_baseline_start=basestart_time)

if len(feature_list):
    og.write_features(feature_list,args[0])

#################
#print all the features in nice format.
features=[k[0:7] for k in og.feature_dict.keys()]
print("file-params       molecule " +' '.join(features)+' ratio')
for fnum,ftuple in enumerate(og.ftuples):
    for imol,mol in enumerate(og.molecules+og.tot_species):
        outputvals=[str(np.round(odict[imol,fnum],3)) for feat,odict in og.mean_feature.items()]
        if og.mean_feature['baseline'][imol,fnum]>1e-5:
            outputvals.append(str(np.round(og.mean_feature['amplitude'][imol,fnum]/og.mean_feature['baseline'][imol,fnum],2)).rjust(8))
        else:
            outputvals.append('inf')
        print(ftuple[1],mol.rjust(16),'  ','  '.join(outputvals))

######################### Plots
if showplot:
    fig,col_inc,scale=pu5.plot_setup(data.molecules,og,len(stimspine),showplot)
    if showplot==2 and len(stimspine):
        figtitle=figtitle+' '+stimspine
    #else:
    #    figtitle=figtitle+' nonspine'
    fig.canvas.set_window_title(figtitle)
    pu5.plottrace(data.molecules,og,fig,col_inc,scale,stimspine,showplot,textsize=textsize)
    #also plot the totaled molecule forms
    if len(tot_species):
        pu5.plot_signature(tot_species,og,figtitle,col_inc,textsize=textsize)    #plot some feature values
    for feat in feature_list:
        pu5.plot_features(og,feat,figtitle)
    if spatial_bins and data.maxvols>1:
        pu5.spatial_plot(data,og)
    if len(mol_pairs):
        pu5.pairs(og,mol_pairs,pairs_timeframe)

'''
5. Move Nadia's pairs plots into plot_h5V2
7. Possibly bring in signature code from sig.py or sig2.py and eliminate one or both of those.
    Create separate class?  Or add weight into total_species.  Add in baseline subtract as option?
    perhaps separate function to compare to thresholds
    pu5.plot3D is used for signatures in sig.py.  How does this differ from spatial_plot?
    lines 341-360 calculates the signature
    lines 370-387 compares to thresholds
    #EXTRACT FEATURES OF total array to add in sig.py functionality
8. possible calculate some feature value relative to a control group (e.g. auc_ratio)
Once working
1. create real arg parser - to avoid position dependence - especially add numstim (or read from h5file)
2. possibly try to figure out how to extract number of stimuli
inj_keyword='<onset>'
#lookup - how to decode b'<onset>46000.0</onset>'
#injtime=
#duration=b'<duration>1000.0</duration>'
#stop = onset+duration
#sort by injtime
#check whether injtime<=stop - if so, separate pulse?  if not continuous

'''
