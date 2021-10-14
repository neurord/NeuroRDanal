# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 10:58:05 2020

@author: kblackw1
#from outside python, type python nrdh5_analv2,py subdir/fileroot -par par1 par2 -mol mol1 mol2 -start 100 200 -tot tot_species_file
from within python, type 
     ARGS="subdir/fileroot -par par1 par2 -mol mol1 mol2 -start 100 200 -tot tot_species_file"
     execfile('path/to/file/nrdh5_anal.py')

#-par:
#  par1 and optionally par2 are specifications of parameter variations, as follows:
#  The filenames to read in are constructed as "subdir/fileroot"+"-"+par1+"*"-"+par2+"*"
#  DO NOT use hyphens in filenames except for preceding parameter name
#-mol:
#  mol1 mol2, etc are the names of molecules to process
#  if mol ommitted, then all molecules processed.  
#-start:
#   start and end time for calculating basal concentration, in sec
#   if sstart ssend are ommitted, then calculates basal from 7.5-10% of runtime
#-tot:
#  tot_species_files containes a list of molecule forms to total, e.g. pPDE10 and pPDE10cAMP to calculate total pPDE10
#if no parameters specified, then fileroot needs to be full filename (excluding the extension)

#e.g. ARGS="../Repo/plc/Model_PLCassay,Ca GaqGTP,Ca GaqGTP Ip3,15 20" time units are sec
#e.g. ARGS="plc/Model_PLCassay_Ca1,Ca Gaq,GTP IP3"

additional parameters lines 27-48
"""
import numpy as np
import sys

from NeuroRDanal import plot_h5V2 as pu5
from NeuroRDanal.nrd_output import nrdh5_output
from NeuroRDanal.nrd_group import nrdh5_group
from NeuroRDanal.h5utilsV2 import parse_args

#probably should add most of these to args with defaults 
submembname='sub'
dendname="other"
spinehead="head"
stimspine=[] #list of stimulated spines
spatial_bins=0  #number of spatial bins to subdivide dendrite to look at spatial gradients
window_size=0.1  #number of msec on either side of peak value to average for maximum
#These control what output is printed or written
show_inject=0
write_output=0 #one file per molecules per input file
output_auc=0#one file per molecule per set of input files
showplot=1 #0 for none, 1 for overall average, 2 for spine concentration, 3 for spine and nonspine on seperate graphs 
show_mol_totals=0
print_head_stats=0
textsize=10
feature_list=[]#,'amplitude']
#these molecules MUST be specified as plot_molecules
mol_pairs=[]#[['CKpCamCa4','ppERK']]#,['ppERK','pSynGap']]
pairs_timeframe=[]#[200,2000] #units are sec
basestart_time=0#2200 #make this value 0 to calculate AUC using baseline calculated from initial time period
aucend=None#600#end time for calculating auc, or None to calculate end time automatically

############## END OF PARAMETERS #################
try:
    args = ARGS.split()
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

params=parse_args(args,do_exit)

if params.mol:
    plot_molecules=params.mol
else:
    plot_molecules=None
    
figtitle=params.fileroot.split('/')[-1]
if params.par:
    figtitle+=' '.join(params.par)

if params.tot:
    import importlib
    tot_name=params.tot
    nm=importlib.import_module(tot_name)
    tot_species=nm.tot_species
    weight=nm.weight
    sub_species=nm.sub_species
else:
    weight={}
    sub_species={}
    tot_species=[]

num_LTP_stim=params.num_stim

og=nrdh5_group(params.fileroot,params.par,tot_species)
for fnum,ftuple in enumerate(og.ftuples):
    data=nrdh5_output(ftuple)
    data.rows_and_columns(plot_molecules,params.start,params.end)
    data.molecule_population()
    #print(data.data['model']['grid'][:])
    if data.maxvols>1:
        data.region_structures(dendname,submembname,spinehead,stimspine) #stimspine is optional
        if spatial_bins>0:
            data.spatial_structures(spatial_bins,dendname)
        data.average_over_voxels()
    # need to add another total array for different regions (to use for signature)
    #Default outputset is _main_, can specify outset= something
    data.total_subspecies(tot_species,sub_species,params.start,weights=weight)
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
og.trace_features(data.trials,window_size,std_factor=1,numstim=num_LTP_stim,end_baseline_start=basestart_time,filt_length=31,aucend=aucend)

if len(feature_list):
    og.write_features(feature_list,params.fileroot,params.write_trials)
#################
#print all the features in nice format.
features=[k[0:7] for k in og.feature_dict.keys()]
print("file-params       molecule " +' '.join(features)+' ratio   CV')
for fnum,ftuple in enumerate(og.ftuples):
    for imol,mol in enumerate(list(og.molecules)+og.tot_species):
        outputvals=[str(np.round(odict[imol,fnum],3)) for feat,odict in og.mean_feature.items()]
        if og.mean_feature['baseline'][imol,fnum]>1e-5:
            outputvals.append(str(np.round(og.mean_feature['amplitude'][imol,fnum]/og.mean_feature['baseline'][imol,fnum],2)).rjust(8))
        else:
            outputvals.append('inf')
        print(ftuple[1],mol.rjust(16),'  ','  '.join(outputvals),np.std(og.feature_dict['auc'][imol,fnum])/og.mean_feature['auc'][imol,fnum])


######################### Plots
if showplot:
    fig,col_inc,scale=pu5.plot_setup(data.molecules,og,len(data.spinelist),showplot)
    if showplot==2 and len(stimspine):
        figtitle=figtitle+' '+' '.join(stimspine)
    if showplot==3:
        for spnum,sp in enumerate(data.spinelist):
            fig[spnum].suptitle(figtitle+' '+sp)
    else:
        fig.canvas.manager.set_window_title(figtitle)
    pu5.plottrace(data.molecules,og,fig,col_inc,scale,data.spinelist,showplot,textsize=textsize)
    #also plot the totaled molecule forms
    if len(tot_species):
        pu5.plot_signature(tot_species,og,figtitle,col_inc,textsize=textsize)    #plot some feature values
    for feat in feature_list:
        pu5.plot_features(og,feat,figtitle)
    if spatial_bins and data.maxvols>1:
        pu5.spatial_plot(data,og)
    if len(mol_pairs):
        pu5.pairs(og,mol_pairs,pairs_timeframe)
    if data.maxvols>1:
        fig,col_inc,scale=pu5.plot_setup(data.molecules,og,len(data.region_dict),3)
        for regnum,reg in enumerate(data.region_dict):
            fig[regnum].suptitle(figtitle+' '+reg)
        pu5.plotregions(data.molecules,og,fig,col_inc,scale,data.region_dict,textsize=textsize)

'''
1. Test total_traces for spatial model
2. Possibly bring in signature code from sig.py or sig2.py and eliminate one or both of those.
    Create separate class?  Or add weight into total_species.  Add in baseline subtract as option?
    perhaps separate function to compare to thresholds
    pu5.plot3D is used for signatures in sig.py.  How does this differ from spatial_plot?
    lines 341-360 calculates the signature
    lines 370-387 compares to thresholds
    #EXTRACT FEATURES OF total array to add in sig.py functionality

3. possibly calculate some feature value relative to a control group (e.g. auc_ratio)

4. possibly try to figure out how to extract number of stimuli
inj_keyword='<onset>'
#lookup - how to decode b'<onset>46000.0</onset>'
#injtime=
#duration=b'<duration>1000.0</duration>'
#stop = onset+duration
#sort by injtime
#check whether injtime<=stop - if so, separate pulse?  if not continuous

'''
