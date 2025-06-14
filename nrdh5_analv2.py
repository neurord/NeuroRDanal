# -*- coding: utf-8 -*-
"""
Created on Thu Apr 23 10:58:05 2020

@author: kblackw1
#from outside python, type python nrdh5_analv2,py subdir/fileroot -par par1 par2 -mol mol1 mol2 -start 100 200 -tot tot_species_file
from within python, type /fileroot -par par1 par2 -mol mol1 mol2 -start 100 200 -tot tot_species_file"
     exec(open(('path/to/file/nrdh5_anal.py').read())

#-par:
#  par1 and optionally par2 are specifications of parameter variations, as follows:
#  The filenames to read in are constructed as "subdir/fileroot"+"-"+par1+"*"-"+parpython  /local/vol00/Users/nminingouzobon/NeuroRDanal/nrdh5_analv2.py /local/vol00/Users/nminingouzobon/cofilin/trials_set/Model_Cof-Abasal*test -par Abasal*test -mol Cof -start 0 3002+"*"
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

#e.g. ARGS='Model_Cof -par HSJCF4trains -mol Cof pCof Cofactin -tot tot_species'
#e.g. ARGS='Model_Cof-HSJCF4trainsspacedcrtl -mol Cof pCof Cofactin -tot tot_species'

additional parameters lines 41-63

"""
#ARGS='/local/vol00/Users/klblackwell/sigpath/nadia_cofilin/Model_Cof -par HSJCF4train*crtl -savedir /local/vol00/Users/klblackwell/sigpath/nadia_cofilin/tmp_out -mol Ca Cof Cofactin RacPAK -start 0 300 -write_trials 1 -tot /local/vol00/Users/klblackwell/sigpath/nadia_cofilin/tot_species_minmax'
import numpy as np
import sys

import plot_h5V2 as pu5
from nrd_output import nrdh5_output
from nrd_group import nrdh5_group
from h5utilsV2 import parse_args,get_tot

#probably should add most of these to args with defaults 
submembname='sub'
dendname="dend"
spinehead="head"
stimspine=['sa1[0]'] #list of stimulated spines
spatial_bins=0  #number of spatial bins to subdivide dendrite to look at spatial gradients
window_size=0.1  #number of msec on either side of peak value to average for maximum
#These control what output is printed or written
show_inject=0
write_output=1#one file per molecules per input file
output_features=1#one file per molecule per set of input files
showplot=3 #0 for none, 1 for overall average, 2 for spine concentration, 3 for spine and nonspine on seperate graph, or for a region plot when there are no spines
show_mol_totals=0
print_head_stats=0
textsize=8
feature_list=[]#['auc','duration']#['duration','auc']
#these molecules MUST be specified as plot_molecules
mol_pairs=[] #[['pCof','actCof'],['CKpCamCa4','PKAphos']]#[['pCof','RacPAK']]#[['CKpCamCa4','ppERK']]#,['ppERK','pSynGap']]
pairs_timeframe=[100,600]#[200,2000] #units are sec
interest_region=['dend','sa1[0]'] #plot mol pairs and write features in which region?
basestart_time=0#2200 #make this value 0 to calculate AUC using baseline calculated from initial time period
aucend=None#600#end time for calculating auc, or None to calculate end time automatically
save_fig=True
zoom = [[0,300],[300,500],[300,1000]]

############## END OF PARAMETERS #################s
try:
    args = ARGS.split()
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False 
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    import os
    mydir=os.getcwd()
    sys.path.append(mydir)
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


tot_species,weight,sub_species,signature,thresh,min_max=get_tot(params)

num_LTP_stim=params.num_stim
iti=params.iti

og=nrdh5_group(params.fileroot,params.par,tot_species,params.savedir)
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
        data.write_average(og.savedir)
    if 'event_statistics' in data.data['trial0']['output'].keys() and show_inject:
        print ("seeds", data.seeds," injection stats:")
        print('molecule             '+'    '.join(data.trials))
        for imol,inject_sp in enumerate(data.data['model']['event_statistics'][:]):
            inject_num=np.column_stack([data.data[t]['output']['event_statistics'][0] for t in data.trials])
            print (inject_sp.split()[-1].rjust(20),inject_num[imol])
    if print_head_stats:
        data.print_head_stats()
if params.write_trials and len(interest_region):
    og.write_trace_trials(interest_region,params.fileroot)
        
#extract some features from the group of data files
#Default numstim = 1, so that parameter not needed for single pulse
#other parameter defaults:  lo_thresh_factor=0.2,hi_thresh_factor=0.8, std_factor=2
#another parameter default: end_baseline_start=0 (uses initial baseline to calculate auc).
#Specify specific sim time near end of sim if initialization not sufficient for a good baseline for auc calculation
og.trace_features(window_size,std_factor=2,numstim=num_LTP_stim,end_baseline_start=basestart_time,filt_length=31,aucend=aucend,iti=iti)
if len(feature_list) and output_features:
    og.write_features(feature_list,params.fileroot,interest_region,params.write_trials)
#################
#print all the features in nice format - Overall values only
regnum=0 #change this to print other regions
features=[k[0:7] for k in og.feature_dict.keys()]
print("file-params       molecule " +' '.join(features)+' ratio   CV')
for fnum,ftuple in enumerate(og.ftuples):
    for imol,mol in enumerate(list(og.molecules)+og.tot_species):
        outputvals=[str(np.round(odict[imol,fnum,regnum],4)) for feat,odict in og.mean_feature.items()]
        if og.mean_feature['baseline'][imol,fnum,regnum]>1e-5:
            outputvals.append(str(np.round(og.mean_feature['amplitude'][imol,fnum,regnum]/og.mean_feature['baseline'][imol,fnum,regnum],2)).rjust(8))
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
    if showplot==3 and data.maxvols>1 and len(data.spinelist)==0:
        fig2,col_inc,scale=pu5.plot_setup(data.molecules,og,len(data.region_dict),3)
        for regnum,reg in enumerate(data.region_dict):
            fig2[regnum].suptitle(figtitle+' '+reg)
        pu5.plotregions(data.molecules,og,fig2,col_inc,scale,data.region_dict,textsize=textsize)
    #also plot the totaled molecule forms
    if len(tot_species):
        figtot=pu5.plot_total_mol(tot_species,og,figtitle,col_inc,textsize=textsize,regions=['sa1[0]','dendsub'])   
    for feat in feature_list:
        pu5.plot_features(og,feat,figtitle)
    if spatial_bins and data.maxvols>1:
        pu5.spatial_plot(data,og)
    if len(mol_pairs):
        pu5.pairs(og,mol_pairs,pairs_timeframe,interest_region)

if len(signature):
     og.norm_sig(signature,thresh,min_max)
     figsig=pu5.plot_signature(og,thresh,figtitle,col_inc,textsize=textsize)    #plot some feature values
     for feature in og.sig_features.keys():
         print('FEATURE:',feature)
         for key in og.sig_features[feature].keys():
             print (key,':', og.sig_features[feature][key])
     if write_output:
        if params.write_trials and len(interest_region):
            og.write_sig(interest_region)
        else:
            og.write_sig()
def ZOOM_fig (figs,zoom,name):
    if not isinstance(figs, list):
        figs = [figs]  # Convert single Figure object to a list
    for f in figs:
        axis=f.axes
        for Z in zoom:
            axis[0].set_xlim(Z[0],Z[1])
            f.savefig(og.savedir+figtitle+'_'+name+'_zoom'+ str(Z[0]) + '_' + str(Z[1]) +'.png')    
if save_fig==True:
    for f,sp in zip(fig,data.spinelist):
        f.set_size_inches((13,13))
        f.savefig(og.savedir+figtitle+'_'+sp+'.png')
    name = ''
    ZOOM_fig(fig,zoom,name)         
    if showplot==3 and data.maxvols>1 and len(data.spinelist)==0:
        for f,reg in zip(fig2,data.region_dict):
            f.savefig(og.savedir+figtitle+'_'+reg+'.png')
        name = 'sp' 
        ZOOM_fig(fig,zoom,name)
    if len(tot_species):
        figtot.savefig(og.savedir+figtitle+'tot.png')
        name = 'tot' 
        ZOOM_fig(figtot,zoom,name)
    if len(signature):
        figsig.savefig(og.savedir+figtitle+'sig.png')
        name = 'sig' 
        ZOOM_fig(figsig,zoom,name)
'''
2. Possibly bring in signature code from sig.py or sig2.py and eliminate one or both of those.
    pu5.plot3D is used for signatures in sig.py.  How does this differ from spatial_plot?

3. possibly calculate some feature value relative to a control group (e.g. auc_ratio)

4. possibly try to figure out how to extract number of stimuli
inj_keyword='<onset>'
#lookup - how to decode b'<onset>46000.0</onset>'
#injtime=
#duration=b'<duration>1000.0</duration>'
#stop = onset+duration
#sort by injtime
#check whether injtime<=stop - if so, separate pulse?  if not continuous
######### Test signature using weight either in total_species or signature
'''
