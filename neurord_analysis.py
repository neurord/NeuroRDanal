from __future__ import print_function
from __future__ import division
#Python version, i.e. alternative of NRDpostAB
#in python, type ARGS="subdir/fileroot,par1 par2,mol1 mol2,sstart ssend,rows" then execfile('neurord_analysis.py')
#DO NOT PUT ANY SPACES NEXT TO THE COMMAS, DO NOT USE TABS, rows is optional
#if mol ommitted, then all molecules processed
#e.g. ARGS="Ca GaqGTP,Ca GaqGTP Ip3,../Repo/plc/Model_PLCassay,15 20" time units are sec
#from outside python, type python neurord_analysis [par1 par2] [mol1 mol2]
#Assumes that molecule outputs are integers, and the hypens used ONLY for parameters
#Can process multiple parameter variations, but all files must use same morphology, and meshfile.  
#It will provide region averages (each spine, dendrite submembrane, cytosol) and if spatialaverage=1,
#will calculate an average of n segments along the dendrite, 
#or whatever structure name is specified in dend variable

import os
import numpy as np
from matplotlib import pyplot
from string import *
import sys  
import glob
from NeuroRDanal import header_parse as hparse
from NeuroRDanal import plot_utils as pu

#######################################################
#indicate the name of the injection spines if you want to exclude them
fake = 'FAKE'
#indicate the name of submembrane region for totaling molecules that are exclusively submembrane
#only relevant for tot_species calculation. this should be name of structure concatenated with sub
submembname='sub'
#Spatial average (=1 to process) only includes the structure "dend", and subdivides into bins:
spatialaverage=0
dend="dend"
spinehead="head"
bins=10
#how much info to print
prnvox=1
prnheader=0
prninfo=0
showss=0
#outputavg determines whether output files are written
outputavg=0
showplot=2
#showplot=1 plots overallmean, showplot=2 allows plotregion to specify region number in the region_list
plotregion=0
##change these endings depending on whether using neurord3.x:
meshend="*mesh.txt.out"
concend='-conc.txt.out'
## or neurord2.x (uncomment these)
meshend="*mesh.txt"
concend='*conc.txt'
#Example of how to total some molecule forms; turn off with tot_species={}
#tot_species={
#        "PKAtot":["PKA", "PKAcAMP2", "PKAcAMP4", "PKAr"],
#        "D1Rtot":["D1R","DaD1R", "GsD1R","DaD1RGs", "pDaD1RGs", "PKAcDaD1RGs"],
#        "pde10tot":["PDE10","pPDE10", "PDE10cAMP","pPDE10cAMP","PKAcPDE10", "PKAcPDE10cAMP"],
#        "Gitot":["Giabg","AChm4RGi","Gim4R", "GaiGTP", "GaiGDP", "ACGai", "ACGasGai", "ACGasGaiATP"],
#        "m4Rtot":["AChm4RGi","Gim4R", "m4R", "AChm4R"]}
tot_species={}
###################################################

Avogadro=6.023e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15

def sortorder(ftuple):
    ans = ftuple[1]
    #print 'sort', ftuple, '->', ans
    return ans

try:
	args = ARGS.split(",")
	print("ARGS =", ARGS, "commandline=", args)
 	do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
	args = sys.argv[1:]
	print("commandline =", args)
	do_exit = True

#1st and 2nd arguements used to construct pattern for reading in multiple files
pattern=args[0]+'*'
if len(args[1]):
        params=args[1].split(" ")
        for par in params:
                pattern=pattern+'-'+par+'*'
else:
        params=[]
whole_pattern=pattern+concend
#A single mesh file means that all files in your list must use the same morphology
meshname=pattern.split('-')[0]+meshend
lastslash=rfind(pattern,'/')
subdir=pattern[0:lastslash]

try:
    data.close()
except Exception:
   pass

###################################################

fnames = glob.glob(whole_pattern)
print("NUM FILES:", len(fnames), "CURRENT DIRECTORY:", os.getcwd(), ", Target directory:", subdir)
if len(fnames)==0:
    print("MESHFILES:", os.listdir(subdir+'/'+meshend))
ss_tot=np.zeros((len(fnames),len(tot_species.keys())))
parlist=[]
if len(args[1]):
        ftuples,parlist=pu.file_tuple(fnames,params)
        ftuples = sorted(ftuples, key=lambda x:x[1])
else:
        ftuples=[(fnames[0],1)]

#Read mesh file to determine how many voxels
if len(fnames)>0:
        meshfile=glob.glob(meshname)[0]
else:
        print("********** no meshfile **************")
maxvols,vox_volume,xloc,yloc,TotVol,deltaY=hparse.read_mesh(meshfile,prninfo)

parval=[]
for fnum,ftuple in enumerate(ftuples):
    fname=ftuple[0]
    parval.append(ftuple[1])
    if fnum == 0:
        f = open(fname, 'r+')
        #read and then parse the header to determine identity/structure of voxels and molecules
        #will not be needed once using hdf5 file, since the region and structure info is part of the mesh information
        data=f.readline()
        if (prnheader==1):
            print("header",data)
        else:
            print("header not printed")
        #UPDATE maxvols, or number of voxels in this function
        regionID,structType,molecules,volnums,maxvols=hparse.header_parse(data,maxvols,prninfo)
        if prninfo:
            print("in neurord_analysis: vox#", volnums)
            print("      regions",regionID)
            print("      structures",structType)
        print("      molecules",molecules)
        f.close()
        #prepare to plot stuff (instead of calculating averages)
        #plot_molecules determines what is plotted
        if len(args[2].split()):
            plot_molecules=args[2].split()
        else:
            plot_molecules=molecules
        if showplot:
            fig,axes,col_inc,scale,numpar=pu.plot_setup(plot_molecules,parlist,params)
            fig.suptitle(pattern.split('/')[-1])
        ss=np.zeros((len(fnames),len(plot_molecules)))
        slope=np.zeros((len(fnames),len(plot_molecules)))
        peaktime=np.zeros((len(fnames),len(plot_molecules)))
        baseline=np.zeros((len(fnames),len(plot_molecules)))
        peakval=np.zeros((len(fnames),len(plot_molecules)))
        lowval=np.zeros((len(fnames),len(plot_molecules)))
        #
        #all voxels should be read in now with labels
        #extract number of unique regions (e.g. dendrite, or sa1[0]), 
        #and create list of subvolumes which contribute to that region
        #will be simpler once using hdf5 file, since the region and structure info is part of the mesh information
        if maxvols>1:
                region_list,region_vox,region_col,region_struct_list,region_struct_vox,region_struct_col=hparse.subvol_list(structType,regionID,volnums,fake)
                dsm_vox=region_struct_list.index(dend+submembname)
                try:
                    head_vox=region_list.index(spinehead)
                except ValueError:
                    head_vox=-1
                RegVol=hparse.region_volume(region_list,region_vox,vox_volume,prnvox)
                RegStructVol=hparse.region_volume(region_struct_list,region_struct_vox,vox_volume,prnvox)
                submembVol=0
                for region in region_list:
                        smname=region+submembname
                        if smname in region_struct_list:
                                submembVol+=RegStructVol[region_struct_list.index(smname)]
                #
                if spatialaverage:
                        hparse.spatial_average(xloc,yloc,bins,regionID,structType,volnums)
     #
    #Read in the data, reshape so that each molecule in separate array
    if len(args)>4:
            time,molecule_array,rows=hparse.readdata(fname,maxvols,molecules,int(args[4]))
    else:
            time,molecule_array,rows=hparse.readdata(fname,maxvols,molecules)
    #
    plot_array=np.zeros((rows,len(plot_molecules)))
    dt=time[1]#/1000
    if len(args)>3:
            sstart = int(float(args[3].split(" ")[0]) // dt)
            ssend = int(float(args[3].split(" ")[1]) // dt)
            if ssend>0.5*rows:
                print("WARNING*******. Possible SS time issue: only", rows, "rows, end time=", time[-1])
            if ssend>rows:
                ssend=int(0.1*rows)
                sstart=int(0.075*rows)
    else:
            sstart=int(0.075*rows)
            ssend=int(0.1*rows)
    ######################################
    #Calculate various region averages, such as soma and dend, subm vs cyt, 
    #use the above lists and volume of each region, and each region-structure
    ######################################
    if maxvols>1:
        data=np.zeros((rows,maxvols),dtype=int)
        for imol in range(len(molecules)):
           if molecules[imol] in plot_molecules:
                data=molecule_array[:,imol,:]
                #calculate region means
                header,RegionMeans=hparse.region_means(data,region_list,region_col,RegVol,time,molecules[imol])
                #calculate region-structure menas
                header2,RegionStructMeans=hparse.region_means(data,region_struct_list,region_struct_col,RegStructVol,time,molecules[imol])
                #calculate overall mean
                OverallMean=np.zeros(len(time))
                for itime in range(len(time)):
                       for k in range(maxvols):
                                OverallMean[itime]+=data[itime,k]
                if (data[:,np.where(structType=='cyt')[0]].all==0):
                        OverallMean[:] /= submembVol*mol_per_nM_u3
                        OverallDensity=OverallMean*deltaY[0]
                else:
                        OverallMean[:] /= TotVol*mol_per_nM_u3
                header='#time ' +header+header2+molecules[imol]+'AvgTot\n'
                #
                plot_index=plot_molecules.index(molecules[imol])
                if showplot==1:
                    plot_array[:,plot_index]=OverallMean
                elif showplot==2:
                    plot_array[:,plot_index]=RegionMeans[:,plotregion]
                ss[fnum,plot_index]=plot_array[sstart:ssend,plot_index].mean()
                #
                #Repeat for spatial averages if specified
                if spatialaverage:
                        spaceheader,SpatialMeans=hparse.region_means(data,range(bins),bincolumns,SpatialVol,time,molecules)
                #
                #write averages to separate files
                if outputavg:
                        outfname=fname[0:-8]+molecules[imol]+'_avg.txt'
                        if molecules[imol] in plot_molecules:
                                print('output file: ', outfname, np.mean(RegionMeans[sstart:ssend],0))
                        outdata=np.column_stack((time,RegionMeans,RegionStructMeans,OverallMean))
                        f=open(outfname, 'w')
                        f.write(header)
                        np.savetxt(f, outdata, fmt='%.4f', delimiter=' ')
                        f.close()
                else:
                    print(molecules[imol].rjust(14), end=' ')
                    if head_vox>-1:
                        print("head ss:%8.4f pk %8.4f " % (RegionMeans[sstart:ssend,head_vox].mean(), RegionMeans[ssend:,head_vox].max()), end=' ')
                    print("dend sm %8.4f pk %8.4f" %((RegionStructMeans[sstart:ssend,dsm_vox].mean()*deltaY[0]), (RegionStructMeans[ssend:,dsm_vox].max()*deltaY[0])))
               #
                #write space
                if spatialaverage:
                        outnamespace=fname[0:-8]+'-'+molecules[imol]+'_space.txt'
                        outdata=np.column_stack((time,SpatialMeans))
                        f=open(outnamespace, 'w')
                        f.write("#time "+spaceheader+'\n')
                        np.savetxt(f, outdata, fmt='%.4f', delimiter=' ')
                        f.close()
    else:
        #no processing needed if only a single voxel.  Just extract, calculate ss, and plot specified molecules
        #0 in 3 index of molecule_array indicates that for 1 voxel structures 0th array has total
        for imol,mol in enumerate(plot_molecules):
                plot_array[:,imol]=molecule_array[:,molecules.index(mol),0]/TotVol/mol_per_nM_u3
                ss[fnum,imol]=plot_array[int(sstart):int(ssend),imol].mean()
    #
    #in both cases (single voxel and multi-voxel):
    #########total some molecule forms - specified by hand above for now
    for imol,mol in enumerate(tot_species.keys()):
           for subspecies in tot_species[mol]:
                   mol_sum=molecule_array[0,molecules.index(subspecies),:].sum()
                   #print imol,mol,subspecies,molecule_array[0,molecules.index(subspecies),:],mol_sum
                   ss_tot[fnum,imol]+=mol_sum/TotVol/mol_per_nM_u3
           print(imol,mol,ss_tot[fnum,imol],"nM, or in picoSD:", ss_tot[fnum,imol]*(TotVol/submembVol)*deltaY[0])
    #####################################################################
    #after main processing, extract a few characteristics of molecule trajectory
    #####################################################################
    print(params, parval[fnum])
    print("      molecule  baseline  peakval  ptime   slope     min     ratio")
    for imol,mol in enumerate(plot_molecules):
        baseline[fnum,imol]=plot_array[sstart:ssend,imol].mean()
        peakpt=plot_array[ssend:,imol].argmax()+ssend
        peaktime[fnum,imol]=peakpt*dt
        peakval[fnum,imol]=plot_array[peakpt-10:peakpt+10,imol].mean()
        lowpt=plot_array[ssend:,imol].argmin()+ssend
        lowval[fnum,imol]=plot_array[lowpt-10:lowpt+10,imol].mean()
        begin_slopeval=0.2*(peakval[fnum,imol]-baseline[fnum,imol])+baseline[fnum,imol]
        end_slopeval=0.8*(peakval[fnum,imol]-baseline[fnum,imol])+baseline[fnum,imol]
        exceedsthresh=np.where(plot_array[ssend:,imol]>begin_slopeval)
        begin_slopept=0
        end_slopept=0
        found=0
        if len(exceedsthresh[0]):
            begin_slopept=np.min(exceedsthresh[0])+ssend
            found=1
            exceedsthresh=np.where(plot_array[begin_slopept:,imol]>end_slopeval)
            if len(exceedsthresh[0]):
                end_slopept=np.min(exceedsthresh[0])+begin_slopept
            else:
                found=0
        if found and len(plot_array[begin_slopept:end_slopept,imol])>1:
                slope[fnum,imol]=(peakval[fnum,imol]-baseline[fnum,imol])/((end_slopept-begin_slopept)*dt)
        else:
                slope[fnum,imol]=-9999
        print(mol.rjust(16),"%8.2f" % baseline[fnum,imol],"%8.2f" %peakval[fnum,imol], end=' ')
        print("%8.2f" % peaktime[fnum,imol], "%8.3f" %slope[fnum,imol], end=' ')  
        print("%8.2f" %lowval[fnum,imol], "%8.2f" %(peakval[fnum,imol]/baseline[fnum,imol]))
    #
    #Now plot some of these molcules, either single voxel or overall average if multi-voxel
    #
    if showplot:
        pu.plottrace(plot_molecules,time,plot_array,parval[fnum],axes,fig,col_inc,scale,parlist)
    #
#then plot the steady state versus parameter value for each molecule
#Needs to be fixed so that it works with non numeric parameter values
if len(params)>1:
        print(np.column_stack((parval,ss)))
        xval=[]#np.zeros(len(parval))
        for i,pv in enumerate(parval):
                if len(parlist[0])>len(parlist[1]):
                        xval.append(pv[0])
                else:
                        xval.append(pv[1])
        print(xval)
        if showss:
                pu.plotss(plot_molecules,xval,ss)
else:
    if showss:
        #also plot the totaled molecule forms
        if len(tot_species.keys()):
                pu.plotss(plot_molecules+tot_species.keys(),parval,np.hstack((ss,ss_tot)))
        else:
                pu.plotss(plot_molecules,parval,ss)

