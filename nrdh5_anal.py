#neurdh5_anal.py
#in python, type ARGS="subdir/fileroot,par1 par2,mol1 mol2,sstart ssend,rows" then execfile('neurdh5_anal.py')
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
#Future improvements:
#   spatial average along dendrite

from __future__ import print_function
from __future__ import division

import os
import numpy as np
from matplotlib import pyplot
from string import *
import sys  
import glob
import h5utils
import plot_h5 as pu5
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
show_inject=0
print_head_stats=0
#outputavg determines whether output files are written
outputavg=0
showplot=1    #2 indicates plot the head conc, 0 means no plots
stimspine='sa1[0]' #"name" of (stimulated) spine
auc_mol='2ag'
endtime=110 #time to stop calculating AUC
textsize=10 #for plots.  Make bigger for presentations

#Example of how to total some molecule forms; turn off with tot_species={}
#No need to specify subspecies if uniquely determined by string
sub_species={"PI": ["Ip3","Ip3degrad","Ip3degPIk","Pip2","PlcCaPip2","PlcCaGqaPip2"],
        "PKA":["PKA", "PKAcAMP2", "PKAcAMP4", "PKAr"]}
tot_species=["D1R","m4R", "m1R","Gi", "Gs", "Gq", "Plc", "AC5", "AC1", "PI", "PKA","D32","PP2A", "PP2B", "PP1", "Cam", "CK", "Pkc", "Dgl","PDE4", "PDE10", "PDE2"]
#tot_species=["Calbin", "CaM", "ncx", "pmca", "CaOut"]
tot_species=[]

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
whole_time_array=[]
for fnum,ftuple in enumerate(ftuples):
    data,maxvols,TotVol,trials,seeds,arraysize,p=h5utils.initialize(ftuple,numfiles,parval)
    if len(p):
        params=p[0]
        parval=p[1]
        parlist=p[2]
    plot_array=[]
    time_array=[]
    #
    ##########################################################
    #   Extract region and structure voxels and volumes
    ##########################################################
    if maxvols>1 and fnum==0:
        molecules=data['model']['species'][:]
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
    #
    ##### Initialization done only for first file in the list
    #
    if fnum==0:
        #initialize plot stuff, arrays for static measurements, and find molecules among the output sets
        if len(args[2].split()):
            plot_molecules=args[2].split()
        else:
            plot_molecules=molecules
        num_mols=len(plot_molecules)
        if numfiles>1:
            whole_plot_array=[[] for mol in plot_molecules]
            whole_time_array=[[] for mol in plot_molecules]
        #
        ss_tot=np.zeros((arraysize,len(tot_species)))
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
        print(params, "=",parval[fnum], "voxels=",maxvols)
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
            #calculate overall mean
            OverallMean=np.zeros((len(trials),np.shape(molecule_pop)[1]))
            OverallMean[:,:]=np.sum(molecule_pop[:,:,:],axis=2)/(TotVol*mol_per_nM_u3)
            header='#time ' +headstruct+headreg+molecule+'AvgTot\n'
            #
            if showplot==2:
                if stimspine in spinelist:
                    spine_index=spinelist.index(stimspine)
                    if fnum==0 and imol==0:
                        figtitle=figtitle+' '+stimspine
                else:
                    spine_index=0
                    if fnum==0 and imol==0:
                        figtitle=figtitle+' '+'nonspine'
                #TEST THIS PART FOR MULTIPLE SPINES/TRIALS
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
            if imol==0:
                print("samples", len(time), "maxtime", time[-1], "conc", np.shape(molecule_pop), np.shape(plot_array), 'time',np.shape(time_array))
            #
            ############# write averages to separate files #######################3
            if outputavg:
                 if molecule in plot_molecules:
                    outfname=fname[0:-3]+'_'+molecule+'_avg.txt'
                    if len(params)==1:
                        param_name=params[0]+parval[fnum]
                    if len(params)==2:
                        param_name=params[0]+parval[fnum][0]+params[1]+parval[fnum][1]
                    print('output file:', outfname,  ' param_name:', param_name)
                    newheader=''
                    newheaderstd=''
                    for item in header.split():
                        if item.startswith('#'):
                            newheader=newheader+item+' '
                        else:
                            newheader=newheader+param_name+'_'+item+' '
                        if not item.startswith('#'):
                            newheaderstd=newheaderstd+param_name+'_'+item+'std '
                    if len(trials)>1:
                        nonspine_out=np.column_stack((RegMeanStd['mean'],RegStructMeanStd['mean'],np.mean(OverallMean,axis=0),RegMeanStd['std'],RegStructMeanStd['std'],np.std(OverallMean,axis=0)))
                    else:
                        newheaderstd=''
                        nonspine_out=np.column_stack((RegionMeans[0,:,:],RegionStructMeans[0,:,:],OverallMean[0,:]))
                    if len(spinelist)>1:
                        newspinehead=''
                        newspineheadstd=''
                        for item in spineheader.split():
                            newspinehead=newspinehead+param_name+'_'+item+' '
                            newspineheadstd=newspineheadstd+param_name+'_'+item+'std '
                        if len(trials)>1:
                            wholeheader=newheader+newheaderstd+newspinehead+newspineheadstd+'\n'
                            outdata=np.column_stack((time,nonspine_out,spineMeanStd['mean'],spineMeanStd['std']))
                        else:
                            wholeheader=newheader+newspinehead+'\n'
                            outdata=np.column_stack((time,nonspine_out,spinemeans[0,:,:]))
                    else:
                        wholeheader=newheader+newheaderstd+'\n'
                        outdata=nonspine_out
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
              time_array.append(time)
              plot_array.append(np.zeros(len(time)))
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
            outset = out_location[mol]['location'].keys()[0]
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
              time_array.append(time)
              plot_array.append(np.zeros(len(time)))
    ######################################
    #Whether 1 voxel or multi-voxel, create plotting array of means for all molecules, all files, all trials
    ##########################################
    if numfiles>1:
        #plot_array dimensions=num molecules x sample times
        #whole_plot_array dimension  = num molecules*num files*sample time
        for mol in range(num_mols):
            whole_plot_array[mol].append(plot_array[mol])
            whole_time_array[mol].append(time_array[mol])
    else:
        #dimensions of plot_array=num molecules x num trials x sample times
        whole_plot_array=plot_array
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
            for subspecie in data['model']['output']['__main__']['species'][:]:
                if mol in subspecie:
                    mol_set.append(subspecie)
        #second, find molecule index of the sub_species and total them
        for subspecie in mol_set:
            mol_index=h5utils.get_mol_index(data,outset,subspecie)
            mol_pop=data['trial0']['output'][outset]['population'][0,:,mol_index]
            ss_tot[fnum,imol]+=mol_pop.sum()/TotVol/mol_per_nM_u3
            if maxvols>1:
                if dsm_vox>-1:
                    dsm_tot[fnum,imol]+=mol_pop[region_struct_dict[dsm_name]['vox']].sum()/region_struct_dict[dsm_name]['vol']*region_struct_dict[dsm_name]['depth']/mol_per_nM_u3
                else:
                    dsm_tot[fnum,imol]+=-1
                if head_index>-1:
                    head_tot[fnum,imol]+=mol_pop[region_dict[spinehead]['vox']].sum()/region_dict[spinehead]['vol']/mol_per_nM_u3
                else:
                    head_tot[fnum,imol]+=-1
        print("Total",mol, end=' ')
        if fnum==0:
            print(mol_set, end=' ')
        print(ss_tot[fnum,imol],"nM")
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
    fig,col_inc,scale=pu5.plot_setup(plot_molecules,parlist,params)
    #need fnames
    fig.canvas.set_window_title(figtitle)
    pu5.plottrace(plot_molecules,whole_time_array,whole_plot_array,parval,fig,col_inc,scale,parlist,textsize)
    #
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
            print('dur=', dur, 'dhpg=', dhpg, 'auc',np.round(auc[pnum][molnum],2), 'ratio', np.round(auc[pnum][molnum]/auc[dhpg0index[j]][molnum]), 'new auc',np.round(newauc[pnum][molnum]), 'ratio', np.round(newauc[pnum][molnum]/newauc[dhpg0index[j]][molnum],2))

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
        if len(tot_species.keys()):
                pu5.plotss(plot_molecules+tot_species.keys(),parval,np.hstack((ss,ss_tot)))
        else:
                pu5.plotss(plot_molecules,parval,ss)

