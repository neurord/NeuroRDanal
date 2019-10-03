from __future__ import print_function
from __future__ import division
import numpy as np
from string import *
import glob
import os
import h5py as h5
from NeuroRDanal import plot_h5 as pu5
from collections import OrderedDict
from orderedmultidict import omdict

#shows how to use tables instead of h5py:  Also how to test for rotation of population table
#f = tables.open_file('...', 'r')
#pop = f.root.trial0.output.__main__.population
#elements_last = getattr(pop.attrs, 'elements_last', 0)
#data.root.model.species will list all species without [:]
#data.root.model.grid[:]['type'] - will list all  voxels - note the switch to [''] after the [:]
#to make it work in python 2 & 3, use np.core.defchararray.decode(text stuff)

def decode(table):
    return np.array([s.decode('utf-8') for s in table])

Avogadro=6.02214179e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15 #0.6022 = PUVC

def new_head(header,param_name):
    newheader=''
    newheaderstd=''
    for item in header.split():
        if item.startswith('#'):
            newheader=newheader+item+' '
        else:
            newheader=newheader+param_name+'_'+item+' '
        if not item.startswith('#'):
            newheaderstd=newheaderstd+param_name+'_'+item+'std '
    return newheader,newheaderstd

def join_params(parval,params):
    if len(params)>1:
        label=join(parval)
    else:
        label=parval
    return label

def initialize(ftuple,numfiles,parval):
    fname=ftuple[0]
    parval.append(ftuple[1])
    data = h5.File(fname,"r")
    maxvols=len(data['model']['grid'])
    TotVol=data['model']['grid'][:]['volume'].sum()
    trials=[a for a in data.keys() if 'trial' in a]
    try:
        seeds=[data[trial].attrs['simulation_seed'] for trial in trials]
    except KeyError:
        seeds=[data[trial]['simulation_seed'][:] for trial in trials]
    outputsets=data[trials[0]]['output'].keys()
    if numfiles==1:
        arraysize=len(trials)
        params=['trial']
        parval=[str(x) for x in range(len(trials))]
        parlist=[parval,[]]
        p=[params,parval,parlist]
    else:
        arraysize=numfiles
        p=[]
    return data,maxvols,TotVol,trials,seeds,arraysize,p

def sstart_end(molecule_list, args, arg_num, out_location,dt,rows):
    num_mols=len(molecule_list)
    sstart=np.zeros((num_mols),dtype=np.int)
    ssend=np.zeros((num_mols),dtype=np.int)
    if len(args)>arg_num:
        for imol,molecule in enumerate(molecule_list):
            if out_location[molecule]!=-1:
                sstart[imol] = float(args[arg_num].split(" ")[0]) // dt[imol]
                ssend[imol] = float(args[arg_num].split(" ")[1]) // dt[imol]
                if ssend[imol]>0.5*rows[imol]:
                    print("WARNING*******. Possible SS time issue: only", rows, "rows")
                if ssend[imol]>rows[imol]:
                    ssend[imol]=0.1*rows[imol]
                    sstart[imol]=0.075*rows[imol]
                    print ("WARNING *****. ssend exceeds sim time, reassigning to ", ssend[imol]*dt)
    else:
        for imol,molecule in enumerate(molecule_list):
            if out_location[molecule]!=-1:
                sstart[imol]=int(0.075*rows[imol])
                ssend[imol]=int(0.1*rows[imol])
    return sstart,ssend

####### FIX/IMPROVE THIS by going back from last outputset, only using outset __main__ if no voxels?
def get_mol_info(simData,plot_molecules,gridpoints):
    outputsets=list(simData['model']['output'].keys())
    dt=np.zeros((len(plot_molecules)))
    samples=np.zeros((len(plot_molecules)),dtype=int)
    out_location={}
    for imol,molecule in enumerate(plot_molecules):
        temp_dict={}
        tot_voxels=0
        for outset in outputsets[1:]:           #better to go backward from last set, and then go to 0 set if mol not found
        #for each in outputsets[-1::-1] and while tot_voxels>grid_points:
            mol_index=get_mol_index(simData,outset,molecule)
            if mol_index>-1:
                samples[imol]=len(simData['trial0']['output'][outset]['times'])
                dt[imol]=simData['trial0']['output'][outset]['times'][1]/1000 #convert msec to sec
                tot_voxels=tot_voxels+len(simData['model']['output'][outset]['elements'])
                temp_dict[outset]={'mol_index':mol_index,'elements':simData['model']['output'][outset]['elements'][:]}
        if len(temp_dict)>0:
            out_location[molecule]={'samples':samples[imol],'dt':dt[imol],'voxels': tot_voxels,'location': temp_dict}
        else:
            outset=outputsets[0]
            print("************* MOLECULE",molecule, " NOT IN REGULAR OUTPUT SETS !")
            mol_index=get_mol_index(simData,outset,molecule)
            if mol_index>-1:
                samples[imol]=len(simData['trial0']['output'][outset]['times'])
                dt[imol]=simData['trial0']['output'][outset]['times'][1]/1000 #convert msec to sec
                temp_dict[outset]={'mol_index':mol_index,'elements':simData['model']['output'][outset]['elements'][:]}
                out_location[molecule]={'samples':samples[imol],'dt':dt[imol],'voxels': gridpoints,'location': temp_dict}
            else:
                out_location[molecule]=-1
                print("** Even Worse: MOLECULE",molecule, " DOES NOT EXIST !!!!!!!!!!!!!")
    return out_location,dt,samples
         
def get_mol_index(simData,outputset,molecule):
    species = decode(simData['model']['output'][outputset]['species'])
    indices=np.where(species == molecule)[0]
    if len(indices) == 1:
        return indices[0]
    else:
        return -1

def get_mol_pop(simData, out_location,gridpoints,trials):
    samples=out_location['samples']
    conc=np.zeros((len(trials),samples,gridpoints))
    for outset in out_location['location'].keys():
        elements=out_location['location'][outset]['elements']
        time=simData[trials[0]]['output'][outset]['times'][:]/1000     #Convert msec to sec
        for trialnum,trial in enumerate(trials):
            tempConc=simData[trial]['output'][outset]['population'][:,:,out_location['location'][outset]['mol_index']]
            #if simulation is still running, array sizes may not be the same. 
            #print('samples', samples, 'tempConc',np.shape(tempConc))
            if np.shape(tempConc)[0]>samples:
                trialConc=np.resize(tempConc,(samples,len(elements)))
                #print ('tempconc',np.shape(trialConc),'too big')
            elif np.shape(tempConc)[0]<samples:
                extrazeros=np.zeros((samples-np.shape(tempConc)[0],len(elements)))
                #print('zeros',np.shape(extrazeros))
                trialConc=np.vstack((tempConc,extrazeros))
                #print('tempconc too small, add',np.shape(extrazeros), 'for', np.shape(trialConc))
            else:
                trialConc=tempConc
            #transpose required to undo the transpose automatically done by python when specifying elements as 3d index
            conc[trialnum,:,elements]=trialConc.T
    return conc,time

def argparse(args):

    #def sortorder(ftuple):
    #    ans = ftuple[1]
    #    #print 'sort', ftuple, '->', ans
    #    return ans
    def sort_paramNum(ftuples,parlist,par):
        print('**********************1:',parlist[1],parlist[0], 'sort_paraNum')
        if np.all([i in '0123456789.' for item in parlist[1] for i in item ]):
            parlist[1]=[float(item) for item in parlist[1]]
            if len(par)>1:
                newftuples=[(tup[0],(tup[1][0],float(tup[1][1]))) for tup in ftuples]
            else:
                newftuples=[(tup[0],float(tup[1])) for tup in ftuples]
            ftuples=sorted(newftuples,key=lambda x:x[1])
        parlist[1]=sorted(parlist[1],key=lambda x:x)
        parlist[0]=sorted(parlist[0],key=lambda x:x)
        return ftuples,parlist

    #sort par
parfloat=[0 for p in params]
for i in range (len(par)):
    if np.all([i in '0123456789.' for item in parlist[i]for i in item ]):
        parlist[i]=[float(item) for item in parlist[i]]
        parfloat[i]=True
    if len(par)==1 and parfloat[i]=True:
        newftuples=[(tup[0],float(tup[1])) for tup in ftuples]
    #1st and 2nd arguements used to construct pattern for reading in multiple files
    pattern=args[0]
    if len(args[1]):
        params=args[1].split(" ")
        for par in params:
            pattern=pattern+'-'+par+'*'
    else:
        params=[]
    whole_pattern=pattern+'.h5'
    print("pattern:", whole_pattern)

    subdir=os.path.dirname(pattern)
    if len(subdir)==0:
        subdir='.'
    fnames = glob.glob(whole_pattern)

    print("files:", fnames)
    print("NUM FILES:", len(fnames), "CURRENT DIRECTORY:", os.getcwd(), ", Target directory:", subdir)
    if len(fnames)==0:
        print("FILES:", *glob.glob(subdir+'/'+'*.h5'), sep='\n')
        raise IOError("no files found")

    parlist=[]
    if len(args[1]):
        ftuples,parlist=pu5.file_tuple(fnames,params)
        ftuples = sorted(ftuples, key=lambda x:x[1])
        print('**********************2:',parlist[1],parlist[0], ftuples, par)
        if len(parlist[1]):
            ftuples,parlist=sort_paramNum(ftuples,parlist,par)
    else:
        star=str.find(pattern,'*')
        if star>-1:
            dash=str.rfind(pattern,'-',0,star)
            params=[pattern[dash+1:star]]
            ftuples,parlist=pu5.file_tuple(fnames,params)
        else:
            ftuples=[(fnames[0],'1')]
    return ftuples,parlist,params

def subvol_list(structType,model):
    #use dictionaries to store voxels corresponding to regions, region_classes (e.g. head) or regions/structures
    region_list=decode(model['regions'])
    region_dict=OrderedDict()
    region_struct_dict=OrderedDict()
    #create dictionary of voxels and volumes for each region
    reg_voxel=omdict(( zip(model['grid'][:]['region'],range(len(model['grid'])) ) ))
    reg_voxel_vol=omdict(( zip(model['grid'][:]['region'],model['grid'][:]['volume'] ) ))
    for regnum in reg_voxel.keys():
        region_dict[region_list[regnum]]={'vox': reg_voxel.allvalues(regnum), 'vol': sum(reg_voxel_vol.allvalues(regnum))}
        # for regions of more than one type, create dictionary of voxels and volumes for each type of each region
        if len(np.unique(model['grid'][reg_voxel.allvalues(regnum)]['type']))>1:
            types = decode(model['grid'][reg_voxel.allvalues(regnum)]['type'])
            struct_voxels=omdict(( zip(types,reg_voxel.allvalues(regnum)) ))
            struct_vox_vol=omdict(( zip(types,reg_voxel_vol.allvalues(regnum)) ))
            for struct in struct_voxels.keys():
                depth=model['grid'][struct_voxels.allvalues(struct)]['y0']-model['grid'][struct_voxels.allvalues(struct)]['y2']
                #Depth is an array.  For submemb, only a single value, for cyt - different values.  Presently only storing one of the values
                key = region_list[regnum] + struct[0:3]
                region_struct_dict[key]={'vox': struct_voxels.allvalues(struct),
                                         'depth':depth[0],
                                         'vol': sum(struct_vox_vol.allvalues(struct))}
    return region_list,region_dict,region_struct_dict

def multi_spines(model):
    spine_dict=OrderedDict()
    #create list of spine voxels
    #first identify all spine voxels and spine labels
    groups=model['grid'][:]['group']
    newgroups=list()
    for n,i in enumerate(groups):
        #in python3 must explicitly decode from byte to string
        if type(i) is np.bytes_:
            groups[n]=i.decode('UTF-8')
            newgroups.append(i.decode('UTF-8'))
        if newgroups[n] =='':
            newgroups[n]='nonspine'
    groups=newgroups
    spine_voxel=omdict((zip(groups,range(len(model['grid'])) ) ))
    spine_voxel_vol=omdict(( zip(groups,model['grid'][:]['volume']) ))
    #create a unique set of spine labels
    newspinelist=spine_voxel.keys()
    for spinenum,spine in enumerate(newspinelist):
        spine_dict[spine]={'vox':spine_voxel.allvalues(spine), 'vol': sum(spine_voxel_vol.allvalues(spine))}
    return newspinelist,spine_dict

def region_means_dict(data,regionDict,time,molecule,trials):
    samples=np.shape(data)[1]
    RegionMeans=np.zeros((len(trials),samples,len(regionDict)))
    header=''       #Header for output file
    for j,item in enumerate(regionDict):
        RegionMeans[:,:,j]=np.sum(data[:,:,regionDict[item]['vox']],axis=2)/(regionDict[item]['vol']*mol_per_nM_u3)
        header=header+molecule+'_'+item+' '       #Header for output file
    RegionMeanStd={}
    if len(trials)>1:
        RegionMeanStd['mean']=np.mean(RegionMeans,axis=0)
        RegionMeanStd['std']=np.std(RegionMeans,axis=0)
    return header,RegionMeans,RegionMeanStd

def spatial_average(bins,region,grid):
    #may want to modify this to use of group instead of region_name
    xloc=[row[0] for row in grid if row['region_name'] == region]
    xdim=np.max(xloc)-np.min(xloc)
    yloc=[row[1] for row in grid if row['region_name'] == region]
    ydim=np.max(yloc)-np.min(yloc)
    if (xdim >= ydim):    #xdim is larger
        loc=xloc
        coord='x0'
        spaceheader='#time, x='
        extraloc=[row[3] for row in grid if row['region_name'] == region]
    else:                #ydim is larger
        loc=yloc
        coord='y0'
        spaceheader='#time, y='
        extraloc=[row[7] for row in grid if row['region_name'] == region]
    minloc=min(np.min(loc),np.min(extraloc))
    maxloc=max(np.max(loc),np.max(extraloc))
    bininc=(maxloc-minloc)/bins
    binmin=[minloc+j*bininc for j in range(bins)]
    binmin.append(maxloc)
    binvoxels=[[] for j in range(bins)]
    binvolumes=[[] for j in range(bins)]
    spatial_dict=OrderedDict()
    for j in range(bins):
        binvoxels[j]=[i for i,row in enumerate(grid) if (row[coord]>=binmin[j] and row[coord]< binmin[j+1]) and row['region_name'] == region]
        binvolumes[j]=[row['volume'] for row in grid if (row[coord]>=binmin[j] and row[coord]< binmin[j+1]) and row['region_name'] == region]
        spatial_dict[str(np.round(binmin[j],3))]={'vox': binvoxels[j], 'vol': sum(binvolumes[j])}
    return spatial_dict

def rolling(indices,dur):
    adjacent=np.diff(indices)==1
    contig_above=[]
    for i,index in enumerate(indices[0:-dur]):
        if adjacent[i:i+dur].all():
            contig_above.append(index)
    return contig_above

def calc_sig(all_sig_array,all_peaks,all_molecules,ltp_molecules,trials,fname_roots,time,col_num):
    sig_ltp=np.zeros((len(fname_roots),len(time),len(col_num)))
    sig_ltd=np.zeros((len(fname_roots),len(time),len(col_num)))

    num_regions=int(len(col_num)/trials)
    #for each file, calculate signature:
    #   sum normalized molecules separately for dend and spine, separately for LTP and LTD
    #        apply to each trial, 1st normalize by baseline subtraction (and peak factor)
    for molnum,mol in enumerate(all_molecules):
        region_peak=np.max(all_peaks[mol],axis=0)
        #print('rg1',mol,region_peak)
        for reg in range(num_regions):
            colset=[tr+reg*trials for tr in range(trials)]
            region_peak[colset]=np.max(region_peak[colset])
        print ('region peak',mol,region_peak,colset)
        for fnum,fname in enumerate(fname_roots):
            if mol in ltp_molecules:
                sig_ltp[fnum]=sig_ltp[fnum]+all_sig_array[mol][fnum]/region_peak
            else:
                sig_ltd[fnum]=sig_ltd[fnum]+all_sig_array[mol][fnum]/region_peak
    return sig_ltp,sig_ltd

def sign_file(col_num,trials,fname_roots,all_sig_array,samplepoints,dt,outfname):
    f=open(outfname, 'w')
    header='filename trial '
    #This for tr in range(trials) and colset=[col*trials+tr for col in range(num_regions)] is specific to sig.py output format
    num_regions=int(len(col_num)/trials)
    for fnum,fname in enumerate(fname_roots):
        for tr in range(trials):
            colset=[col*trials+tr for col in range(num_regions)]
            samples=[]
            for mol in sorted(all_sig_array.keys()):
                #print(fname,mol)
                if fnum==0 and tr==0:
                    for col in range(num_regions):
                        header=header+str(mol)+'_'+str(int((samp_point[0]+samp_point[1])*0.5*dt))+'c'+str(col)+' '
                        extracted_data=all_sig_array[mol][fnum][samp_point[0]:samp_point[1]].T
                    samples.append(list(np.mean(extracted_data[colset],axis=1)))
            if fnum==0 and tr==0:
                f.write(header+'\n')
            outputrow=[ np.round(val,2) for sublist in samples for val in sublist]
            out=str(tr)+' '+" ".join(str(e) for e in outputrow)
            f.write(fname[fname.find('-'):]+' '+out+'\n')
        #print(fname,out)
    f.close()
    return

def sig_outfile_multi(time_thresh,dt,samp_times,pstrt,pend,fname_roots,fname_suffix,col_num,trials,all_sig_array):
    for dur in time_thresh:  
        samplepoints=[(pstrt,pend)] 
        win=int(float(dur)/dt/2) 
        for st in samp_times: 
            t=int(float(st)/dt)+pstrt-win
            samplepoints.append((t,t+2*win)) #alternative 1 - if multiple sample points, create single file with multiple points
	#
        sample_strings=[str(int((point[0]-pstrt+win)*dt)) for point in samplepoints[1:]]
        sample_str='_'+'and'.join(sample_strings)+'_'
        outfname=fname_roots[0].split('-')[0]+sample_str+str(dur)+fname_suffix+'.txt'
        f=open(outfname, 'w')
        #
        sign_file(col_num,trials,fname_roots,all_sig_array,samplepoints,dt,outfname)
    return

def sig_outfile_single(time_thresh,time,samp_times,pstart,pend,fname_roots,fname_suffix,col_num,trials,all_sig_array):
    for dur in time_thresh:  
        samplepoints=[(pstrt,pend)] 
        win=int(float(dur)/time[1]/2) 
        for st in samp_times: 
            t=int(float(st)/time[1])+pstrt-win
            samplepoints=[(pstrt,pend),(t,t+2*win)] #alternative 2 - if multiple sample points, create multiple files
	    #
            sample_strings=[str(int((point[0]-pstrt+win)*time[1])) for point in samplepoints[1:]]
            sample_str='_'+'and'.join(sample_strings)+'_'
            outfname=fname_roots[0].split('-')[0]+sample_str+str(dur)+fname_suffix+'.txt'
            #
            sign_file(col_num,trials,fname_roots,all_sig_array,samplepoints,time,outfname)
    return

