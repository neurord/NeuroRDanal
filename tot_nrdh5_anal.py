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
from matplotlib import pyplot as plt
import sys  
import glob
from NeuroRDanal import h5utils
from NeuroRDanal import plot_h5 as pu5
import h5py as h5
import csv

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
normalized=0
csvfile=0 ######to save data to be compare.. please use only one tot_species for now
show_inject=0
print_head_stats=0
#outputavg determines whether output files are written
outputavg=0
showplot=0    #2 indicates plot the head conc, 0 means no plots
stimspine='sa1[0]' #"name" of (stimulated) spine
auc_mol='2ag'
endtime=150 #time to stop calculating AUC - make shorter than entire duration if simulations are dropping below basal
textsize=35 #for plots.  Make bigger for presentations

#Example of how to total some molecule forms; turn off with tot_species={}
#No need to specify subspecies if uniquely determined by string
sub_species={'ras':['rasGap','RasGTPGap'], 'rap':['rap1Gap', 'Rap1GTPGap'],'Ras': ['pShcGrb2SosRas', 'CamCa4GRFRas', 'Raf1Ras', 'dRaf1Ras','dRaf1RasMEK', 'dRaf1RaspMEK','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK', 'bRafRas', 'bRafRasMEK','bRafRaspMEK', 'RasGTP', 'RasGDP', 'RasSynGap', 'RasGTPGap', 'RaspSynGap'],'Rap1GTP':['bRafRap1MEK', 'bRafRap1pMEK', 'bRafRap1', 'Raf1Rap1', 'Rap1GTP','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK'],'PKA':['PKA', 'PKAcAMP2', 'PKAcAMP4', 'PKAr'], 'erk':['ppERK','pERK'], 'RasGTP':['Raf1Ras', 'dRaf1Ras', 'dRaf1RasMEK', 'dRaf1RaspMEK', 'bRafRas', 'bRafRasMEK', 'bRafRaspMEK', 'RasGTP','dRaf1bRaf','dRaf1bRafMEK','dRaf1bRafpMEK'], 'RasSyn':['RasSynGap', 'RaspSynGap'], 'Rap1Syn':['Rap1SynGap', 'Rap1pSynGap'], 'Ca':['Ca'], 'cAMP':['cAMP'],'ERK':['ppERK']}

tot_species=['ppERK']


#tot_species=['erk','ERK', 'MKP1','MEK','PP2A','bRaf', 'Raf1', 'Ras', 'Rap1', 'Syn', 'ras', 'rap', 'GRF', 'PKA', 'Cam','CK','cAMP','Ca','Gbg','Src','Epac','Sos','Cbl']

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
        if numfiles>1:
            whole_plot_array=[[] for mol in plot_molecules]
            whole_space_array=[[] for mol in plot_molecules]
            whole_time_array=[[] for mol in plot_molecules]
        #
        #ss_tot=np.zeros((arraysize,len(tot_species)))
        ss_tot=[[[] for i in range (arraysize)] for j in tot_species]
        #ss_time_array=[[[] for n in range(arraysize)] for m in range(len(tot_species))]
        slope=np.zeros((arraysize,num_mols))
        peaktime=np.zeros((arraysize,num_mols))
        baseline=np.zeros((arraysize,num_mols))
        auc=np.zeros((arraysize,num_mols))
        peakval=np.zeros((arraysize,num_mols))
        lowval=np.zeros((arraysize,num_mols))
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
                    param_name=params[0]+parval[fnum]
                if len(params)==2:
                    param_name=params[0]+parval[fnum][0]+params[1]+parval[fnum][1]
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
        for mol in plot_molecules:
          if out_location[mol]!=-1:
            outset = list(out_location[mol]['location'].keys())[0]
            imol=out_location[mol]['location'][outset]['mol_index']
            tempConc=np.zeros((len(trials),out_location[mol]['samples']))
            time_array.append(data[trials[0]]['output'][outset]['times'][:]/1000)
            #generate output files for these cases
            for trialnum,trial in enumerate(trials):
                tempConc[trialnum]=data[trial]['output'][outset]['population'][:,voxel,imol]/TotVol/mol_per_nM_u3
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
                param_name=params[0]+parval[fnum]
            if len(params)==2:
                param_name=params[0]+parval[fnum][0]+params[1]+parval[fnum][1]
            print('output file:', outfname,  ' param_name:', param_name)
            f=open(outfname, 'w')
            f.write('time    '+os.path.basename(ftuple[0]).split('.')[0]+'_'+mol+'\n')
            np.savetxt(f, np.row_stack((time_array,plot_array)).T, fmt='%.4f', delimiter=' ')
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
            for subspecie in molecules:
                if mol in subspecie:
                    mol_set.append(subspecie)
        #second, find molecule index of the sub_species and total them
        #print('mol_set',mol_set)
        time=data[trials[0]]['output'][outset]['times'][:]/1000
        tot_pop=np.zeros((len(time)))
        for subspecie in mol_set:
            mol_index=h5utils.get_mol_index(data,outset,subspecie)
            #mol_pop=data['trial0']['output'][outset]['population'][0,:,mol_index]
            mol_pop=data['trial0']['output'][outset]['population'][:,:,mol_index]  #read everything
                             
            #ss_tot[fnum,imol]+=mol_pop.sum()/TotVol/mol_per_nM_u3
            tot_pop+=np.sum(mol_pop, axis=1)/TotVol/mol_per_nM_u3  #add accross voxel
            
            
            if maxvols>1:
                if dsm_vox>-1:
                    dsm_tot[fnum,imol]+=mol_pop[region_struct_dict[dsm_name]['vox']].sum()/region_struct_dict[dsm_name]['vol']*region_struct_dict[dsm_name]['depth']/mol_per_nM_u3
                else:
                    dsm_tot[fnum,imol]+=-1
                if head_index>-1:
                    head_tot[fnum,imol]+=mol_pop[region_dict[spinehead]['vox']].sum()/region_dict[spinehead]['vol']/mol_per_nM_u3
                else:
                    head_tot[fnum,imol]+=-1
        ss_tot[imol][fnum]=tot_pop
        print("Total",mol, end=' ')
        if fnum==0:
            print(mol_set, end=' ')
        print(ss_tot[imol][fnum][0],"nM")
        if maxvols>1:
            print(" or head:",head_tot[fnum,imol],"nM, or dsm:", dsm_tot[fnum,imol], "picoSD")
#
#####################################################################
#after main processing, extract a few characteristics of molecule trajectory
#####################################################################
endpt=int(endtime/dt[0])
for pnum in range(arraysize):
    print(params, parval[pnum])
    print("        molecule  baseline  peakval   ptime    slope      min     ratio")
    for imol,mol in enumerate(plot_molecules):
      if out_location[mol]!=-1:
        window=int(window_size/dt[imol])
        baseline[pnum,imol]=whole_plot_array[imol][pnum][sstart[imol]:ssend[imol]].mean()
        peakpt=whole_plot_array[imol][pnum][ssend[imol]:].argmax()+ssend[imol]
        auc[pnum,imol]=np.sum(whole_plot_array[imol][pnum][ssend[imol]:endpt]-baseline[pnum,imol])*dt[imol]
        peaktime[pnum,imol]=peakpt*dt[imol]
        peakval[pnum,imol]=whole_plot_array[imol][pnum][peakpt-window:peakpt+window].mean()
        lowpt=whole_plot_array[imol][pnum][ssend[imol]:].argmin()+ssend[imol]
        lowval[pnum,imol]=whole_plot_array[imol][pnum][lowpt-10:lowpt+10].mean()
        begin_slopeval=0.2*(peakval[pnum,imol]-baseline[pnum,imol])+baseline[pnum,imol]
        end_slopeval=0.8*(peakval[pnum,imol]-baseline[pnum,imol])+baseline[pnum,imol]
        exceedsthresh=np.where(whole_plot_array[imol][pnum][ssend[imol]:]>begin_slopeval)
        begin_slopept=0
        end_slopept=0
        found=0
        if len(exceedsthresh[0]):
            begin_slopept=np.min(exceedsthresh[0])+ssend[imol]
            found=1
            exceedsthresh=np.where(whole_plot_array[imol][pnum][begin_slopept:]>end_slopeval)
            if len(exceedsthresh[0]):
                end_slopept=np.min(exceedsthresh[0])+begin_slopept
            else:
                found=0
        if found and len(whole_plot_array[imol][pnum][begin_slopept:end_slopept])>1:
                slope[pnum,imol]=(peakval[pnum,imol]-baseline[pnum,imol])/((end_slopept-begin_slopept)*dt[imol])
        else:
                slope[pnum,imol]=-9999
        print(mol.rjust(16),"%8.2f" % baseline[pnum,imol],"%8.2f" %peakval[pnum,imol], end=' ')
        print("%8.2f" % peaktime[pnum,imol], "%8.3f" %slope[pnum,imol], "%8.2f" %lowval[pnum,imol], end=' ')
        if baseline[pnum,imol]>1e-5:
            print("%8.2f" %(peakval[pnum,imol]/baseline[pnum,imol]))
        else:
            print("   inf")
#
#####################################################################
#Now plot some of these molcules, either single voxel or overall average if multi-voxel
#####################################################################
#
if showplot:
    fig,col_inc,scale=pu5.plot_setup(plot_molecules,parlist,params,len(stimspine.split()),showplot)
    #need fnames
    fig.canvas.set_window_title(figtitle)
    pu5.plottrace(plot_molecules,whole_time_array,whole_plot_array,parval,fig,col_inc,scale,parlist,textsize,stimspine.split(),showplot,)
    #
if spatialaverage:
    pu5.space_avg(plot_molecules,whole_space_array,whole_time_array,parval,spatial_dict)
#
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
ss=baseline
if len(params)>1:
        #print(np.column_stack((parval,ss)))
        xval=[]
        for i,pv in enumerate(parval):
                if len(parlist[0])>len(parlist[1]):
                        xval.append(pv[0])
                else:
                        xval.append(pv[1])
        print(xval)
        if showss:
                pu5.plotss(plot_molecules,xval,ss)
else:
    if showss:
        #also plot the totaled molecule forms
       fig,col_inc,scale=pu5.plot_setup(tot_species,parlist,params,len(stimspine.split()),showplot)
       fig.canvas.set_window_title(figtitle)
       pu5.plottrace(tot_species,ss_time_array,ss_tot,parval,fig,col_inc,scale,parlist,textsize,stimspine.split(),showplot) 
        #if len(tot_species):
         #      fig=plt.figure() 
        #for imol,mol in enumerate (tot_species):
         #   plt.plot(time,ss_tot[imol][fnum],label=mol)
        #fig.canvas.set_window_title(figtitle)
        #plt.xlabel('Time(sec)')
        #plt.ylabel('Conc(nM)')
        #plt.legend()
    #else:
     #   pu5.plotss(plot_molecules,parval,ss)


     ###normalized data to baseline based on arraysize
#for imol, mol in enumerate(tot_species):
#    fold_change=[[[] for i in range (arraysize)] for j in tot_species]
#for fnum in range(arraysize):
#    init_val=np.mean(ss_tot[imol][fnum][sstart[imol]:ssend[imol]]) #get true ss based, time before stim
#    fold_change[imol][fnum]=ss_tot[imol][fnum]/init_val

    
if normalized: 
#        fig,col_inc,scale=pu5.plot_setup(tot_species,parlist,params,len(stimspine.split()),showplot)
#        fig.canvas.set_window_title(figtitle)
#        pu5.plottrace(tot_species,ss_time_array,fold_change,parval,fig,col_inc,scale,parlist,textsize,stimspine.split(),showplot)
        
    for imol, mol in enumerate(tot_species):
        fig_norm=plt.figure()
        fold_change=[[[] for i in range (arraysize)] for j in tot_species]
        ###normalized data to baseline
        for fnum in range(arraysize):
           init_val=np.mean(ss_tot[imol][fnum][sstart[imol]:ssend[imol]])
           fold_change[imol][fnum]=ss_tot[imol][fnum]/init_val
           #plt.plot(time,fold_change[imol][fnum], label=xval)
           fig,col_inc,scale=pu5.plot_setup(fold_change,parlist,params,len(stimspine.split()),showplot)
           fig.canvas.set_window_title(figtitle)
           fig.suptitle('PKA pathway induces higher ppERK fold_change amplitude', fontweight='bold',fontsize=30) #title is very specifc-change to fit yours
           pu5.plottrace(tot_species,ss_time_array,fold_change,parval,fig,col_inc,scale,parlist,textsize,stimspine.split(),showplot)
           plt.xlabel('Time(sec)',fontweight='bold')
           plt.ylabel('fold_change'+mol,fontweight='bold')
           #plt.legend(parlist,parval)


           if csvfile:
               ##save data as a csv file to compare with other data using the graphing.py file
               myData=np.column_stack([time,fold_change[imol][fnum]])
               np.savetxt(args[0].split('/')[-1]+'-'+args[1]+'_'+mol+'.csv',myData, header='time(sec), '+mol, delimiter=',') ###replace arg[0] with file name that differ with diff files see how Dr Blackwell did
 
 

