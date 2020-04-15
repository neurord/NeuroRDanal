#neurdh5_anal.py
#in python, type ARGS="subdir/fileroot,par1 par2,mol1 mol2,sstart ssend,rows" then execfile('path/to/file/nrdh5_anal.py')
#DO NOT PUT ANY SPACES NEXT TO THE COMMAS, DO NOT USE TABS, rows is optional
#mol1 mol2, etc are the names of molecles to process
#par1 and optionally par2 are specifications of parameter variations, as follows:
#The filenames to read in are constructed as "subdir/fileroot"+"-"+par1+"*"-"+par2+"*"
#DO NOT use hyphens in filenames except for preceding parameter name
#if no parameters specified, then fileroot needs to be full filename (excluding the extension)
#e.g. ARGS="../Repo/plc/Model_PLCassay,Ca GaqGTP,Ca GaqGTP Ip3,15 20" time units are sec
#e.g. ARGS="plc/Model_PLCassay_Ca1,Ca Gaq,GTP IP3"
#if mol ommitted, then all molecules processed.  if sstart ssend are ommitted, then calculates basal from 7.5-10% of runtime
#in the first set of parameters below, change outputavg from 0 to 1 to generate output files of region averages for plotting
#from outside python, type python neurordh5_analysis "subdir/fileroot [par1 par2] [mol1 mol2]"
#Assumes that molecule outputs are integers
#Can process multiple parameter variations, but all files must use same morphology, meshfile, and simulation time/sampling rate.
#Can process multiple trials for each parameter variation
#It will provide region averages (each spine, dendrite submembrane, cytosol)
#If only a single file, will plot multiple trials; if multiple trials, plots the mean over trials

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
#indicate the name of submembrane region for totaling molecules that are exclusively submembrane
#only relevant for tot_species calculation. this should be name of structure concatenated with sub
submembname='sub'
dend="dend"
spinehead="head"
window_size=1  #number of seconds on either side of peak value to average for maximum
#Spatial average (=1 to process) only includes the structure dend, and subdivides into bins:
spatialaverage=0
bins=10
#how much info to print
showss=0
showplateau=0
showpairs=0
show_inject=0
print_head_stats=0
#outputavg determines whether output files are written
outputavg=0
outputauc=1
showplot=0    #2 indicates plot the head conc, 0 means no plots
stimspine='sa1[0]' #"name" of (stimulated) spine
auc_mol='2ag'
endtime=2000 #time to stop calculating AUC - make shorter than entire duration if simulations are dropping below basal, set to -1 for entire duration 
textsize=20 #for plots.  Make bigger for presentations
basestarttime=2200

#Example of how to total some molecule forms; turn off with tot_species={}
#No need to specify subspecies if uniquely determined by string
sub_species={'ras':['rasGap','RasGTPGap'], 'rap':['rap1Gap', 'Rap1GTPGap'],'Ras': ['pShcGrb2SosRas', 'CamCa4GRFRas', 'Raf1Ras', 'dRaf1Ras','dRaf1RasMEK', 'dRaf1RaspMEK','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK', 'bRafRas', 'bRafRasMEK','bRafRaspMEK', 'RasGTP', 'RasGDP', 'RasSynGap', 'RasGTPGap', 'RaspSynGap'],'Rap1GTP':['bRafRap1MEK', 'bRafRap1pMEK', 'bRafRap1', 'Raf1Rap1', 'Rap1GTP','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK'],'PKA':['PKA', 'PKAcAMP2', 'PKAcAMP4', 'PKAr'], 'erk':['ppERK','pERK'], 'RasGTP':['Raf1Ras', 'dRaf1Ras', 'dRaf1RasMEK', 'dRaf1RaspMEK', 'bRafRas', 'bRafRasMEK', 'bRafRaspMEK', 'RasGTP','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK'], 'RasSyn':['RasSynGap', 'RaspSynGap'], 'Rap1Syn':['Rap1SynGap', 'Rap1pSynGap'],'cAMP': ['EpacAMP', 'cAMP','PDE4cAMP','PDE2cAMP', 'PDE2cAMP2', 'PKAcAMP2', 'PKAcAMP4'], 'Ca':['Ca'],'ERK':['pERK', 'ppERK', 'pERKMKP1', 'ppERKMKP1', 'ppMEKERK', 'ppMEKpERK', 'ppERKpShcGrb2Sos'], 'free_Syn':['SynGap', 'pSynGap']}

tot_species=[]

#molecules that we want to check if there is any correlation by plotting them together 
mol_pairs=[['ppERK','Ca'], ['ppERK','CKpCamCa4'],['cAMP','ppERK']]
#starting time for pairing 
plot_start_time=10
#ending time for paring
plot_end_time=300
###################################################

Avogadro=6.023e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15

try:
    args = ARGS.split(",")
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

ftuples,parlist,params=h5utils.argparse(args)
figtitle=args[0].split('/')[-1]+args[1]

try:
    data.close()
except Exception:
    pass

###################################################
parval=[]
numfiles=len(ftuples)
whole_plot_array=[]
whole_space_array=[]
whole_time_array=[]
for fnum,ftuple in enumerate(ftuples):
    data,maxvols,TotVol,trials,seeds,arraysize,p=h5utils.initialize(ftuple,numfiles,parval)
    print (ftuple, 'volume',TotVol)
    if len(p):
        params=p[0]
        parval=p[1]
        parlist=p[2]
    space_array=[]
    plot_array=[]
    time_array=[]
    #
    ##########################################################
    #   Extract region and structure voxels and volumes
    ##########################################################
    if maxvols>1 and fnum==0:
        structType=data['model']['grid'][:]['type']
        region_list,region_dict,region_struct_dict=h5utils.subvol_list(structType,data['model'])
        #Replace the following with test for whether there is more than one "group"
        try:
            head_index=list(region_dict.keys()).index(spinehead)
        except ValueError:
            head_index=-1
        if head_index>0:
            #create "group" dictionary-spinelist & spinevox, with one more for nonspine voxels
            spinelist,spinevox=h5utils.multi_spines(data['model'])
        else:
            spinelist=''
        #
        dsm_tot=np.zeros((arraysize,len(tot_species)))
        head_tot=np.zeros((arraysize,len(tot_species)))
        dsm_name=dend+submembname
        try:
            dsm_vox=list(region_struct_dict.keys()).index(dsm_name)
        except ValueError:
            dsm_vox=-1
        if spatialaverage:
            spatial_dict=h5utils.spatial_average(bins,dend,data['model']['grid'])
            vox=[x['vox'] for x in spatial_dict.values()]
            if any(v==0 for v in vox):
                print ("**********Too many bins for spatial average****************")
#
    ##### Initialization done only for first file in the list
    #
    if fnum==0:
        molecules=h5utils.decode(data['model']['species'][:])
        #initialize plot stuff, arrays for static measurements, and find molecules among the output sets
        if len(args[2].split()):
            plot_molecules=args[2].split()
        else:
            plot_molecules=molecules
        num_mols=len(plot_molecules)
        auc_array=[[[] for ii in ftuples] for jj in plot_molecules]
        if numfiles>1:
            whole_plot_array=[[] for mol in plot_molecules] 
            whole_space_array=[[] for mol in plot_molecules]
            whole_time_array=[[] for mol in plot_molecules]
        #
        ss_tot=[[[] for n in range(arraysize)] for m in range(len(tot_species))]
        ss_time_array=[[[] for n in range(arraysize)] for m in range(len(tot_species))]
        ss_tot_zero=np.zeros((arraysize,len(tot_species)))
        #
    ######################################
    #Calculate various region averages, such as soma and dend, subm vs cyt, spines
    ######################################
    #
    # do this for each file since they may have different number of samples, or different locations
    out_location,dt,rows=h5utils.get_mol_info(data,plot_molecules,maxvols)
    #
    #Which "rows" should be used for baseline value, specifed in args[3].  If different for each file then will have problems later
    sstart,ssend=h5utils.sstart_end(plot_molecules,args,3,out_location,dt,rows)
    molecule_name_issue=0
    if maxvols>1:
        for imol,molecule in enumerate(plot_molecules):
          if out_location[molecule]!=-1:
            molecule_pop,time=h5utils.get_mol_pop(data,out_location[molecule],maxvols,trials)
            time_array.append(time)
            #calculate region means
            headstruct,RegionMeans,RegMeanStd=h5utils.region_means_dict(molecule_pop,region_dict,time,molecule,trials)
            #calculate region-structure means
            headreg,RegionStructMeans,RegStructMeanStd=h5utils.region_means_dict(molecule_pop,region_struct_dict,time,molecule,trials)
            #if more than one spine, calculate individual spine means
            if len(spinelist)>0:
                spineheader,spinemeans,spineMeanStd=h5utils.region_means_dict(molecule_pop,spinevox,time,molecule,trials)
            else:
                spineheader=''
            if spatialaverage:
                spacehead,spaceMeans,spaceMeanStd=h5utils.region_means_dict(molecule_pop,spatial_dict,time,molecule,trials)
                #calculate overall mean
            OverallMean=np.zeros((len(trials),np.shape(molecule_pop)[1]))
            OverallMean[:,:]=np.sum(molecule_pop[:,:,:],axis=2)/(TotVol*mol_per_nM_u3)
            header='#time ' +headstruct+headreg+molecule+'AvgTot\n'
            #
            if showplot==2:
                spine_index=[spinelist.index(stimsp) for stimsp in stimspine.split()]
                if len(spine_index):
                    if fnum==0 and imol==0:
                        figtitle=figtitle+' '+stimspine
                else:
                    spine_index=0
                    if fnum==0 and imol==0:
                        figtitle=figtitle+' '+'nonspine'
                if numfiles>1:
                    plot_array.append(np.mean(spinemeans,axis=0)[:,spine_index])
                else:
                    plot_array.append(spinemeans[:,:,spine_index])
            else:
                if numfiles>1:
                    plot_array.append(np.mean(OverallMean,axis=0))
                     #plot_array dimensions=number of molecules x sample times
                else:
                    #dimensions of plot_array=num molecules x num trials x sample times
                    plot_array.append(OverallMean)
            if spatialaverage:
                if numfiles>1:
                        space_array.append(np.mean(spaceMeans,axis=0))
                else:
                        space_array.append(spaceMeans)
            if imol==0:
                print(parval[fnum], "voxels",maxvols, "samples", len(time), "maxtime", time[-1], "conc", np.shape(molecule_pop), np.shape(plot_array), 'time',np.shape(time_array))
            #
            ############# write averages to separate files #######################3
            if outputavg==1:
                outfname=ftuple[0][0:-3]+'_'+molecule+'_avg.txt'
                if len(params)==1:
                    param_name=params[0]+str(parval[fnum])
                if len(params)==2:
                    param_name=params[0]+str(parval[fnum][0])+params[1]+str(parval[fnum][1])
                print('output file:', outfname,  ' param_name:', param_name)
                newheader, newheaderstd=h5utils.new_head(header,param_name)
                if len(trials)>1:
                    nonspine_out=np.column_stack((RegMeanStd['mean'],RegStructMeanStd['mean'],np.mean(OverallMean,axis=0),RegMeanStd['std'],RegStructMeanStd['std'],np.std(OverallMean,axis=0)))
                else: 
                    newheaderstd=''
                    nonspine_out=np.column_stack((RegionMeans[0,:,:],RegionStructMeans[0,:,:],OverallMean[0,:]))
                if len(spinelist)>1:
                    newspinehead, newspineheadstd=h5utils.new_head(spineheader,param_name)
                    if len(trials)>1:
                        wholeheader=newheader+newheaderstd+newspinehead+newspineheadstd+'\n'
                        outdata=np.column_stack((time,nonspine_out,spineMeanStd['mean'],spineMeanStd['std']))
                    else:
                        wholeheader=newheader+newspinehead+'\n'
                        outdata=np.column_stack((time,nonspine_out,spinemeans[0,:,:]))
                else:
                    wholeheader=newheader+newheaderstd+'\n'
                    outdata=np.column_stack((time,nonspine_out))
                f=open(outfname, 'w')
                f.write(wholeheader)
                np.savetxt(f, outdata, fmt='%.4f', delimiter=' ')
                f.close()


            if print_head_stats:
                print(molecule.rjust(14), end=' ')
                if head_index>-1:
                    if len(spinelist)>1:
                        stimspinenum=list(spinelist).index(stimspine)
                        headmean=np.mean(np.mean(spinemeans[:,sstart[imol]:ssend[imol],stimspinenum],axis=0),axis=0)
                        headmax=np.mean(spinemeans[:,sstart[imol]:ssend[imol],stimspinenum],axis=0).max()
                    else:
                        headmean=np.mean(RegionMeans[:,sstart[imol]:ssend[imol],head_index])
                        tempmax=np.max(RegionMeans[:,ssend[imol]:,head_index],axis=1)
                        headmax=np.mean(tempmax)
                    print("head ss:%8.4f pk %8.4f " % (headmean, headmax), end=' ')
                if dsm_vox>-1:
                    dsm_max=np.max(RegionStructMeans[:,ssend[imol]:,dsm_vox],axis=1)
                    print("dend sm %8.4f pk %8.4f" %((RegionStructMeans[:,sstart[imol]:ssend[imol],dsm_vox].mean()*region_struct_dict[dsm_name]['depth']),
                                                     (np.mean(dsm_max)*region_struct_dict[dsm_name]['depth'])))
          else:
              if fnum==0 and molecule_name_issue==0:
                  print("Choose molecules from:", molecules)
                  molecule_name_issue=1
              time_array.append([])
              plot_array.append([])
              #
    else:
        ######################################
        #minimal processing needed if only a single voxel.
        #Just extract, calculate ss, and plot specified molecules
        #might want to create output files with mean and stdev
        ######################################
        voxel=0
        for i,mol in enumerate (plot_molecules):
          if out_location[mol]!=-1:
            outset = list(out_location[mol]['location'].keys())[0]
            imol=out_location[mol]['location'][outset]['mol_index']
            tempConc=np.zeros((len(trials),out_location[mol]['samples']))
            time_array.append(data[trials[0]]['output'][outset]['times'][:]/1000)
            #generate output files for these cases
            for trialnum,trial in enumerate(trials):
               tempConc[trialnum]=data[trial]['output'][outset]['population'][:,voxel,imol]/TotVol/mol_per_nM_u3
            auc_array[i][fnum]=tempConc
            if numfiles>1:
                 plot_array.append(np.mean(tempConc,axis=0))
                 #plot_array dimensions=numfiles x number of molecules x sample times
            else:
                #plot_array dimensions=number of molecules (x number of trials) x sample times
                plot_array.append(tempConc)
          else:
              if fnum==0 and molecule_name_issue==0:
                  print("Choose molecules from:", molecules)
                  molecule_name_issue=1
              time_array.append([])
              plot_array.append([])
          if outputavg==1:
            #This output is needed to extract txt file from h5 file for plotting
            outfname=ftuple[0][0:-3]+'_'+mol+'_avg.txt'
            if len(params)==1:
                param_name=params[0]+str(parval[fnum])
            if len(params)==2:
                param_name=params[0]+str(parval[fnum][0])+params[1]+str(parval[fnum][1])
            print('output file:', outfname,  ' param_name:', param_name, np.shape(time_array[-1]),np.shape(plot_array[-1]))
            f=open(outfname, 'w')
            f.write('time    '+os.path.basename(ftuple[0]).split('.')[0]+'_'+mol+'\n')
            np.savetxt(f, np.row_stack((time_array[-1],plot_array[-1])).T, fmt='%.4f', delimiter=' ')
            f.close()
    ######################################
    #Whether 1 voxel or multi-voxel, create plotting array of means for all molecules, all files, all trials
    ##########################################
    if numfiles>1:
        #plot_array dimensions=num molecules x sample times
        #whole_plot_array dimension  = num molecules*num files*sample time
        for mol in range(num_mols):
            whole_plot_array[mol].append(plot_array[mol])
            whole_time_array[mol].append(time_array[mol])
            if spatialaverage:
                whole_space_array[mol].append(space_array[mol])
    else:
        #dimensions of plot_array=num molecules x num trials x sample times
        whole_plot_array=plot_array
        whole_space_array=space_array
        whole_time_array=[[time_array[imol] for trial in trials] for imol in range(len(plot_molecules))]
    if 'event_statistics' in data['trial0']['output'].keys() and show_inject:
        print ("seeds", seeds," injection stats:")
        if len(trials)>1:
            for each in trials:
                each_data=data[each]['output']['event_statistics'][0]
                for inject_sp,inject_num in zip(data['model']['event_statistics'][:],each_data):
                    #trials_avg=np.zeros((len(trials),np.shape(inject_sp.split())))
                    #trials_avg=np.mean(each_data[:],axis=0)
                    print (inject_sp.split()[-1].rjust(20),inject_num[:])#,trials_avg)
        else:
            for inject_sp,inject_num in zip(data['model']['event_statistics'][:],data['trial0']['output']['event_statistics'][0]):
                print (inject_sp.split()[-1].rjust(20),inject_num[:])
    #
    ###################################################
    #   in both cases (single voxel and multi-voxel):
    #   total some molecule forms, to verify initial conditions
    ###################################################
    #
    outset="__main__"
    for imol,mol in enumerate(tot_species):
        mol_set=[]
        #first set up arrays of all species (sub_species) that are a form of the molecule
        if mol in sub_species.keys():
            mol_set=sub_species[mol]
        else:
            mol_set=[sub for sub in h5utils.decode(data['model']['output']['__main__']['species'][:]) if mol in sub]
        ss_time=data[trials[0]]['output'][outset]['times'][:]/1000
        ss_tmp=np.zeros(len(ss_time))
        #second, find molecule index of the sub_species and total them
        for subspecie in mol_set:
            mol_index=h5utils.get_mol_index(data,outset,subspecie)
            mol_pop=data['trial0']['output'][outset]['population'][:,:,mol_index]
            ss_tmp+=np.sum(mol_pop,axis=1)/TotVol/mol_per_nM_u3
            if maxvols>1:
                if dsm_vox>-1:
                    dsm_tot[fnum,imol]+=mol_pop[region_struct_dict[dsm_name]['vox']].sum()/region_struct_dict[dsm_name]['vol']*region_struct_dict[dsm_name]['depth']/mol_per_nM_u3
                else:
                    dsm_tot[fnum,imol]+=-1
                if head_index>-1:
                    head_tot[fnum,imol]+=mol_pop[region_dict[spinehead]['vox']].sum()/region_dict[spinehead]['vol']/mol_per_nM_u3
                else:
                    head_tot[fnum,imol]+=-1
        ss_tot[imol][fnum]=ss_tmp
        ss_time_array[imol][fnum]=ss_time
        ss_tot_zero[fnum,imol]=ss_tmp[0]
        print("Total",mol, end=' ')
        if fnum==0:
            print(mol_set, end=' ')
        print(ss_tot_zero[fnum,imol],"nM")
        if maxvols>1:
            print(" or head:",head_tot[fnum,imol],"nM, or dsm:", dsm_tot[fnum,imol], "picoSD")
#
#####################################################################
#after main processing, extract a few characteristics of molecule trajectory
#####################################################################

dur_thresh=np.zeros((num_mols))
amplitude=np.zeros((num_mols,arraysize))
slope=np.zeros((num_mols,arraysize))
peaktime=np.zeros((num_mols,arraysize))
baseline=np.zeros((num_mols,arraysize))
baseline_auc=np.zeros((num_mols,arraysize,len(trials)))
auc=np.zeros((num_mols,arraysize,len(trials)))
auc_mean=np.zeros((num_mols,arraysize))
auc_std=np.zeros((num_mols,arraysize))
peak=np.zeros((num_mols,arraysize,len(trials)))
peak_mean=np.zeros((num_mols,arraysize))
peak_std=np.zeros((num_mols,arraysize))
baseline_std=np.zeros((num_mols,arraysize,len(trials)))
auc_thresh=np.zeros((num_mols,arraysize,len(trials)))
peakval=np.zeros((num_mols,arraysize))
dur_plateau=np.zeros((num_mols,arraysize))
lowval=np.zeros((num_mols,arraysize))
half_max=np.zeros((num_mols,arraysize))
end_dur=np.zeros((num_mols,arraysize))
start_dur=np.zeros((num_mols,arraysize))
for pnum in range(arraysize):
    print(params, parval[pnum])
    print("        molecule  baseline  peakval   ptime    slope      min     ratio")
    for imol,mol in enumerate(plot_molecules):
      if out_location[mol]!=-1:
        if endtime==-1:
              endpt=len(whole_plot_array[imol][pnum])
        else:
              endpt=int(endtime/dt[imol])
        window=int(window_size/dt[imol])
        baseline[imol,pnum]=whole_plot_array[imol][pnum][sstart[imol]:ssend[imol]].mean()
        peakpt=whole_plot_array[imol][pnum][ssend[imol]:].argmax()+ssend[imol]
        peaktime[imol,pnum]=peakpt*dt[imol]
        peakval[imol,pnum]=whole_plot_array[imol][pnum][peakpt-window:peakpt+window].mean()
        lowpt=whole_plot_array[imol][pnum][ssend[imol]:].argmin()+ssend[imol]
        lowval[imol,pnum]=whole_plot_array[imol][pnum][lowpt-10:lowpt+10].mean()
        amplitude[imol,pnum]=np.round(peakval[imol][pnum]-baseline[imol][pnum],1)
        begin_slopeval=0.2*(amplitude[imol][pnum])+baseline[imol][pnum] 
        end_slopeval=0.8*(amplitude[imol][pnum])+baseline[imol][pnum]
        exceedsthresh=np.where(whole_plot_array[imol][pnum][ssend[imol]:]>begin_slopeval)[0]+ssend[imol]
        begin_slopept=0
        end_slopept=0
        found=0
        #
        t=int((parval[pnum]*3)//dt[imol]+ssend[imol])
        end_auc=np.zeros(len(trials))
        tempauc=np.zeros(len(trials))

        
        for trialnum,trial in enumerate(trials):
            basestart=int(basestarttime/dt[imol])
            baseline_auc[imol,pnum,trialnum]=auc_array[imol][pnum][trialnum][basestart:].mean()
            baseline_std[imol,pnum,trialnum]=auc_array[imol][pnum][trialnum][basestart:].std()
            auc_thresh[imol,pnum,trialnum]=baseline_auc[imol][pnum][trialnum]+0*baseline_std[imol][pnum][trialnum]####
            peakpt_auc=auc_array[imol][pnum][trialnum][ssend[imol]:].argmax()+ssend[imol]
            peakpt_t=auc_array[imol][pnum][trialnum][t:].argmax()+t
            peak[imol,pnum,trialnum]=auc_array[imol][pnum][trialnum][peakpt_auc-window:peakpt_auc+window].mean()
            peak_mean[imol,pnum]=peak[imol,pnum].mean()
            peak_std[imol,pnum]=peak[imol,pnum].std()
            belowthresh_auc=np.where(auc_array[imol][pnum][trialnum][peakpt_t:]<auc_thresh[imol][pnum][trialnum])[0]+peakpt_t
            if len(belowthresh_auc):
                end_auc[trialnum]=np.min(belowthresh_auc)
            else:
                print ('********* ERK is not returning to basal, raise your threshold by 2**********')
            tempauc[trialnum]=np.sum(auc_array[imol][pnum][trialnum][ssend[imol]:int(end_auc[trialnum])]-baseline_auc[imol][pnum][trialnum])*dt[imol]
        auc[imol,pnum]=tempauc
        auc_mean[imol,pnum]=tempauc.mean()
        auc_std[imol,pnum]=tempauc.std()
        
        #
        half_max[imol,pnum]=0.2*(amplitude[imol][pnum])+baseline[imol][pnum]
        belowthresh=np.where(whole_plot_array[imol][pnum][ssend[imol]:peakpt]<half_max[imol][pnum])[0]+ssend[imol]
        if len(belowthresh):
            start_dur[imol,pnum]=(np.max(belowthresh))*dt[imol]
        belowthresh=np.where(whole_plot_array[imol][pnum][peakpt:]<half_max[imol,pnum])[0]+peakpt
        if len(belowthresh):
            end_dur[imol,pnum]=(np.min(belowthresh))*dt[imol]
           
        #
        if len(exceedsthresh):
            begin_slopept=np.min(exceedsthresh)
            begin_durtime=(np.min(exceedsthresh))
            dt[imol]
            end_durtime=(np.max(exceedsthresh[0]))*dt[imol]
            #print(np.round(begin_durtime,1),np.round( end_durtime,1))
            found=1
            exceedsthresh2=np.where(whole_plot_array[imol][pnum][begin_slopept:]>end_slopeval)[0]+begin_slopept
            if len(exceedsthresh2):
                end_slopept=np.min(exceedsthresh2)
            else:
                found=0
        if found and len(whole_plot_array[imol][pnum][begin_slopept:end_slopept])>1:
                slope[imol,pnum]=(peakval[imol,pnum]-baseline[imol][pnum])/((end_slopept-begin_slopept)*dt[imol])
        else:
                slope[imol,pnum]=-9999
        print(mol.rjust(16),"%8.2f" % baseline[imol,pnum],"%8.2f" %peakval[imol,pnum], end=' ')
        print("%8.2f" % peaktime[imol,pnum], "%8.3f" %slope[imol,pnum], "%8.2f" %lowval[imol,pnum], end=' ') 
        if baseline[imol,pnum]>1e-5:
            print("%8.2f" %(peakval[imol,pnum]/baseline[imol,pnum]))
        else:
            print("   inf")
#extract prolong/plateau data 
duration=end_dur-start_dur

if outputauc==1:
    if params==2:
        outfname=args[0]+'-'+'analysis'+'-'+args[1].split(' ')[0]+'-'+args[1].split(' ')[1]+'.txt'
    else:
       outfname=args[0]+'-'+'analysis'+'-'+args[1]+'.txt'         
    if len(params)<2:
        header='parval  ' #added speced at the end to have spac in between
        if len(params)==0:
            p=args[0]
        else:
            p=parval
    else:
        header='parval1  parlval2 '
        p=parval
    header+='           '.join(plot_molecules)+'_'+'auc_mean  '+'           '.join(plot_molecules)+'_'+'auc_std  '+'           '.join(plot_molecules)+'_'+'peak_mean '+'           '.join(plot_molecules)+'_'+'peak_std '+'\n'
    #header+='           '.join(plot_molecules)+'_'+'peak'+'\n'
    f=open(outfname, 'w')
    f.write(header)
    np.savetxt(f,np.column_stack((p,np.round(auc_mean[0]/1000,3),np.round(auc_std[0]/1000,3),np.round(peak_mean[0]/1000,3))),fmt='%1s', delimiter='  ')
    f.close() 

#####################################################################
#Now plot some of these molcules, either single voxel or overall average if multi-voxel
#####################################################################
#
if showplot:
    fig,col_inc,scale=pu5.plot_setup(plot_molecules,parlist,params,len(stimspine.split()),showplot)
    #need fnames
    fig.canvas.set_window_title(figtitle)
    pu5.plottrace(plot_molecules,whole_time_array,whole_plot_array,parval,fig,col_inc,scale,parlist,textsize,stimspine.split(),showplot)
    #
if spatialaverage:
    pu5.space_avg(plot_molecules,whole_space_array,whole_time_array,parval,spatial_dict)
#fig.suptitle('', fontweight='bold',fontsize=40)
#pyplot.xlabel('Time(sec)', fontweight='bold')
#pyplot.ylabel('ppERK (nM)', fontweight='bold')

#plot plateau
def plateau ():
    if len(parlist[0])>len(parlist[1]):
        par_index=0
    else:
        par_index=1
    for imol,mol in enumerate(plot_molecules):
        pyplot.figure(figtitle)
        pyplot.plot(parlist,duration[imol], label=mol)
        pyplot.xlabel('Inj_dur')
        pyplot.ylabel('duration'+'_'+ mol)
    pyplot.legend()

if showplateau==1:
    plateau()


#This code is very specific for the Uchi sims where there are two parameters: dhpg and duration
#it will work with other parameters, as long as there are two of them. Just change the auc_mol
if auc_mol and auc_mol in plot_molecules and 'dhpg' in params:
    newauc=np.zeros((arraysize,num_mols))
    molnum=plot_molecules.index(auc_mol)
    newbaseline=baseline[:,molnum].mean()
    for pnum in range(arraysize):
        for imol,mol in enumerate(plot_molecules):
            if out_location[mol]!=-1:
                newauc[pnum,imol]=np.sum(whole_plot_array[imol][pnum][ssend[imol]:endpt]-newbaseline)*dt[imol]
    dhpg0index=np.zeros(len(parlist[0]))
    for i,dhpg in enumerate(np.sort(parlist[1])):
        for j,dur in enumerate(np.sort(parlist[0])):
            pnum=parval.index((dur,dhpg))
            if i==0:
                dhpg0index[j]=pnum
            if i==0 and j==0:
                print ('{} dur dhpg auc   ratio  new_auc ratio'.format(args[3]))
            print('{0:8} {1:4} {2:.2f} {3:.3f} {4:.2f} {5:.3f}'.format(dur, dhpg, auc[pnum][molnum], auc[pnum][molnum]/auc[int(dhpg0index[j])][molnum], newauc[pnum][molnum], newauc[pnum][molnum]/newauc[int(dhpg0index[j])][molnum]))

#then plot the steady state versus parameter value for each molecule
#Needs to be fixed so that it works with non numeric parameter values
#is ss the baseline?  Or measuring at some other time point?
if len(params)>1:
    xval=[]
    for i,pv in enumerate(parval):
        if len(parlist[0])>len(parlist[1]):
            xval.append(pv[0])
        else:
            xval.append(pv[1])
    print(xval)
else:
    xval=parval
if showss:
    #also plot the totaled molecule forms
    fig,col_inc,scale=pu5.plot_setup(tot_species,parlist,params,len(stimspine.split()),showplot)
    fig.canvas.set_window_title(figtitle)
    pu5.plottrace(tot_species,ss_time_array,ss_tot,parval,fig,col_inc,scale,parlist,textsize,stimspine.split(),showplot)
    #pu5.plotss(plot_molecules,xval,baseline)



#plot pairs
def pairs ():
    plot_start=int(plot_start_time/dt[0])
    plot_end=int(plot_end_time/dt[0])
    for pair in mol_pairs:
        print(pair)
        do_plot=True
        if pair[0] in plot_molecules:
            molY=plot_molecules.index(pair[0])
        else:
            do_plot=False
        if pair[1] in plot_molecules:
            molX=plot_molecules.index(pair[1])
        else:
             do_plot=False
        if do_plot:
            pyplot.figure()
            pyplot.title('---'.join(pair))
            for pnum in range(arraysize):
                X=whole_plot_array[molX][pnum]
                Y=whole_plot_array[molY][pnum]
                time_vectorY=np.linspace(0,whole_time_array[0][0][-1],len(Y))
                time_vectorX=np.linspace(0,whole_time_array[0][0][-1],len(X))
                # check if molX & moly same length
                if len(X)==len(Y):
                    pyplot.plot(X[plot_start:plot_end],Y[plot_start:plot_end], label=xval[pnum], linestyle='--')
                if len(X)>len(Y):
                    molX_interp=np.interp(time_vectorY,time_vectorX,X)
                    pyplot.plot(molX_interp[plot_start:plot_end],Y[plot_start:plot_end], label=xval[pnum], linestyle='--')
                if len(Y)>len(X):
                    molY_interp=np.interp(time_vectorX,time_vectorY,Y)
                    pyplot.plot(X[plot_start:plot_end],molY_interp[plot_start:plot_end], label=xval[pnum], linestyle='--')
            pyplot.legend()
            pyplot.xlabel(pair[1])
            pyplot.ylabel(pair[0])
        else:
            print('*********************Molecule not in ARGS****************', pair)
if showpairs==1:
    pairs()
    
    #####################################################################

