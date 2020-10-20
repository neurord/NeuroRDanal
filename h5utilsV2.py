from __future__ import print_function
from __future__ import division
import numpy as np
from string import *
import glob
import os
import h5py as h5
from collections import OrderedDict
from orderedmultidict import omdict

def decode(table):
    return np.array([s.decode('utf-8') for s in table])

Avogadro=6.02214179e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15 #0.6022 = PUVC
ms_to_sec=1000

def join_params(parval,params):
    if len(params)>1:
        label=parval.join('-')
    else:
        label=parval
    return label

def sstart_end(molecule_list, out_location,dt,rows,args=None):
    num_mols=len(molecule_list)
    sstart={}#np.zeros((num_mols),dtype=np.int)
    ssend={}#np.zeros((num_mols),dtype=np.int)
    if args is not None:
        for imol,mol in enumerate(molecule_list):
            if out_location[mol]!=-1:
                sstart[mol] = int(np.round(float(args.split(" ")[0]) / dt[mol]))
                ssend[mol] = int(np.round(float(args.split(" ")[1]) / dt[mol]))
                if ssend[mol]>0.5*rows[imol]:
                    print("WARNING*******. Possible SS time issue: only", rows, "rows")
                if ssend[mol]>rows[imol]:
                    ssend[mol]=0.1*rows[imol]
                    sstart[mol]=0.075*rows[imol]
                    print ("WARNING *****. ssend exceeds sim time, reassigning to ", ssend[mol]*dt)
    else:
        for imol,mol in enumerate(molecule_list):
            if out_location[molecule]!=-1:
                sstart[mol]=int(0.075*rows[imol])
                ssend[mol]=int(0.1*rows[imol])
    return sstart,ssend

def get_mol_info(simData,plot_molecules):
    outputsets=list(simData['model']['output'].keys())
    dt={}#np.zeros((len(plot_molecules)))
    samples=np.zeros((len(plot_molecules)),dtype=int)
    num_outsets=len(outputsets)
    out_location={}
    for imol,molecule in enumerate(plot_molecules):
        temp_dict={}
        tot_voxels=0
        for setnum,outset in enumerate(reversed(outputsets)): #will be problem if molecule in several sets that have different dt
            if setnum<num_outsets-1 or len(temp_dict)==0:
                mol_index=get_mol_index(simData,outset,molecule)
                if mol_index>-1:
                    samples[imol]=len(simData['trial0']['output'][outset]['times'])
                    dt[molecule]=simData['trial0']['output'][outset]['times'][1]/ms_to_sec #convert msec to sec
                    tot_voxels=tot_voxels+len(simData['model']['output'][outset]['elements'])
                    temp_dict[outset]={'mol_index':mol_index,'elements':simData['model']['output'][outset]['elements'][:]}
        if len(temp_dict)>0:
            out_location[molecule]={'samples':samples[imol],'dt':dt[molecule],'voxels': tot_voxels,'location': temp_dict}
        else:
            out_location[molecule]=-1
            print("** MOLECULE", molecule, " DOES NOT EXIST !!!!!!!!!!!!!")
    return out_location,dt,samples
         
def get_mol_index(simData,outputset,molecule):
    species = decode(simData['model']['output'][outputset]['species'])
    indices=np.where(species == molecule)[0]
    if len(indices) == 1:
        return indices[0]
    else:
        return -1

def get_mol_pop(Data,mol):
    samples=Data.out_location[mol]['samples']
    counts=np.zeros((len(Data.trials),samples,Data.maxvols))
    for outname,outset in Data.out_location[mol]['location'].items():
        time=Data.data[Data.trials[0]]['output'][outname]['times'][:]/ms_to_sec     #Convert msec to sec
        for trialnum,trial in enumerate(Data.trials):
            tempcounts=Data.data[trial]['output'][outname]['population'][:,:,outset['mol_index']]
            #if simulation is still running, array sizes may not be the same. 
            if np.shape(tempcounts)[0]>samples:
                trialcounts=np.resize(tempcounts,(samples,len(outset['elements'])))
            elif np.shape(tempcounts)[0]<samples:
                extrazeros=np.zeros((samples-np.shape(tempcounts)[0],len(outset['elements'])))
                trialcounts=np.vstack((tempcounts,extrazeros))
            else:
                trialcounts=tempcounts
            #transpose required to undo the transpose automatically done by python when specifying elements as 3d index
            counts[trialnum,:,outset['elements']]=trialcounts.T
    return counts,time

def argparse(args):

    def check_for_float(parlist):
        testset=[i for item in parlist for i in item ]
        if np.all([i in '0123456789.' for i in testset]) and len(testset):
            parlist=[float(item) for item in parlist]
            parlist=sorted(parlist,key=lambda x:x)
            par_is_float=True
        else:
            par_is_float=False
        return parlist,par_is_float
        
    def sort_paramNum(ftuples,parlist):
        parlist[0],par0_is_float=check_for_float(parlist[0])
        if len(parlist[1])>0:
            parlist[1],par1_is_float=check_for_float(parlist[1])
            if par1_is_float and not par0_is_float:
                newftuples=[(tup[0],(tup[1][0],float(tup[1][1]))) for tup in ftuples]
            if par0_is_float and not par1_is_float:  
                newftuples=[(tup[0],(float(tup[1][0]),tup[1][1])) for tup in ftuples]
            if par0_is_float and par1_is_float: 
                newftuples=[(tup[0],(float(tup[1][0]),float(tup[1][1]))) for tup in ftuples]
			if not par0_is_float and not par1_is_float:
                newftuples=[(tup[0],(tup[1][0],tup[1][1])) for tup in ftuples]
        elif par0_is_float:
            newftuples=[(tup[0],tuple((float(tup[1]),))) for tup in ftuples]
            par1_is_float=False
        else:
            par1_is_float=False
        if par0_is_float or par1_is_float:
		ftuples=sorted(newftuples,key=lambda x:x[1])
	else:
		ftuples=newftuples
        return ftuples,parlist

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
        ftuples,parlist=file_tuple(fnames,params)
        ftuples = sorted(ftuples, key=lambda x:x[1])
        if len(parlist[1]) or len(parlist[0]):
            ftuples,parlist=sort_paramNum(ftuples,parlist)
    else:
        star=str.find(pattern,'*')
        if star>-1:
            dash=str.rfind(pattern,'-',0,star)
            params=[pattern[dash+1:star]]
            ftuples,parlist=file_tuple(fnames,params)
            ftuples = sorted(ftuples, key=lambda x:x[1])
        else:
            ftuples=[(fnames[0],('1'))]
    return ftuples,parlist,params #list of filenames with params, list of just params

def file_tuple(fnames,params):
    #Extract value of parameter from filename and put into parameter list
    ftuple=[]
    par0list=[]
    par1list=[]
    print('pu: ',params, fnames)
    for fname in fnames:
          dotloc=fname.rfind('.')
          if params[0]=='*':
               split_text='-'
          else:
               split_text='-'+params[0]
          part_fname=fname[0:dotloc].split(split_text)[-1] 
          hyphen=part_fname.find('-')
          if hyphen>-1:
               parval0=part_fname[0:hyphen]
          else:
               parval0=part_fname 
          if (parval0 not in par0list):
               par0list.append(parval0)
          if len(params)>1:
               parval1=part_fname.split('-'+params[1])[-1] 
               ftuple.append((fname,(parval0,parval1)))
               if (parval1 not in par1list):
                    par1list.append(parval1)
               print('ft: fname: {}, par0: {}, par1:{}'.format(fname,par0list, par1list))
          else:
               ftuple.append((fname,parval0))
               print('ft: fname: {}, par0: {}'.format(fname,par0list))
    return ftuple,[par0list,par1list] #list of filenames with params, list of just params

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

def region_means_dict(Data,molecule,regionDict):
    samples=len(Data.time[molecule])
    RegionMeans=np.zeros((len(Data.trials),samples,len(regionDict)))
    header=''       #Header for output file
    pars=[str(q) for q in Data.parval] #list of strings representing parameters
    for j,item in enumerate(regionDict.keys()):
        #print(item,regionDict[item]['vox'],regionDict[item]['vol'],np.shape(Data.counts[molecule]))
        RegionMeans[:,:,j]=np.sum(Data.counts[molecule][:,:,regionDict[item]['vox']],axis=2)/(regionDict[item]['vol']*mol_per_nM_u3)
        header=header+molecule+'_'+'_'.join(pars)+'_'+item+' '       #Header for output file
    return header,RegionMeans

def spatial_average(bins,region,grid):
    #may want to modify this to use of group instead of region_name
    xloc=[row[0] for row in grid if row['region_name'].decode('utf-8') == region]
    xdim=np.max(xloc)-np.min(xloc)
    yloc=[row[1] for row in grid if row['region_name'].decode('utf-8') == region]
    ydim=np.max(yloc)-np.min(yloc)
    if (xdim >= ydim):    #xdim is larger
        loc=xloc
        coord='x0'
        spaceheader='#time, x='
        extraloc=[row[3] for row in grid if row['region_name'].decode('utf-8') == region]
    else:                #ydim is larger
        loc=yloc
        coord='y0'
        spaceheader='#time, y='
        extraloc=[row[7] for row in grid if row['region_name'].decode('utf-8') == region]
    minloc=min(np.min(loc),np.min(extraloc))
    maxloc=max(np.max(loc),np.max(extraloc))
    bininc=(maxloc-minloc)/bins
    binmin=[minloc+j*bininc for j in range(bins)]
    binmin.append(maxloc)
    binvoxels=[[] for j in range(bins)]
    binvolumes=[[] for j in range(bins)]
    spatial_dict=OrderedDict()
    for j in range(bins):
        binvoxels[j]=[i for i,row in enumerate(grid) if (row[coord]>=binmin[j] and row[coord]< binmin[j+1]) and row['region_name'].decode('utf-8') == region]
        binvolumes[j]=[row['volume'] for row in grid if (row[coord]>=binmin[j] and row[coord]< binmin[j+1]) and row['region_name'].decode('utf-8') == region]
        spatial_dict[str(np.round(binmin[j],3))]={'vox': binvoxels[j], 'vol': sum(binvolumes[j])}
    return spatial_dict

def rolling(indices,dur):
    adjacent=np.diff(indices)==1
    contig_above=[]
    for i,index in enumerate(indices[0:-dur]):
        if adjacent[i:i+dur].all():
            contig_above.append(index)
    return contig_above

