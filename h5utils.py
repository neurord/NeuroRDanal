from __future__ import print_function
from __future__ import division
import numpy as np
from string import *
import glob
import os
import h5py as h5
import plot_h5 as pu5
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

Avogadro=6.023e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15

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

    def sortorder(ftuple):
        ans = ftuple[1]
        #print 'sort', ftuple, '->', ans
        return ans

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

    lastslash=str.rfind(pattern,'/')
    if lastslash > -1:
        subdir=pattern[0:lastslash]
    else:
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
    for n,i in enumerate(groups):
        if i =='':
            groups[n]='nonspine'
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
    
 
