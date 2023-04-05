
'''ARGS="model -start stime etime -IC icfile -Rxn rxnfile"
## model is name of h5 file from which you are going to take initial conditions
## stime etime specify the timeframe from which to obtain molecule concentrations for IC
## IC is the name of the IC file to be updated
## Rxn file is the reaction file used for the simulation
###exec(open('/path/to/directory/UpdateIC_basal_spatial.py').read())
# from terminal window,
python3 UpdateIC_basal_spatial.py h5file '' ''  'stime etime' ICfile Rxnfile
do not include file extensions in the filenames
e.g.
python3 -i UpdateIC_basal_spatial.py Model_SPNspineAChm4R_Gshydr5_GapD-nostim -start 450 500 -IC IC_SPNspine_noATP -Rxn Rxn_SPNspine_noATP
'''
from lxml import etree
from xml.etree import ElementTree as ET
import os
import sys
import numpy as np
import h5py as h5
from NeuroRDanal import h5utilsV2
from NeuroRDanal.nrd_output import nrdh5_output
from NeuroRDanal.nrd_group import nrdh5_group


def AllSpecies (root):
    all_species ={}
    cytosol_species=[]
    submembrane_species=[]
    for specie_sub in root.findall('.//PicoSD'):
        submembrane_species.append(specie_sub.get("specieID"))
    all_species['submemb']=np.unique(submembrane_species)
    for specie_cyt in root.findall('.//NanoMolarity'):
        cytosol_species.append(specie_cyt.get("specieID"))
    all_species['cytosol']=np.unique(cytosol_species)# no repet mol name 
    return all_species

def DiffuseSpecies (Rxn_filename):
    species_diff=[];species_nondiff=[]
    tree_rxn=ET.parse(Rxn_filename+'.xml')
    root_rxn=tree_rxn.getroot()
    for elem_rxn in root_rxn:
        if elem_rxn.tag== 'Specie':
            if float(elem_rxn.get('kdiff'))>0:
                species_diff.append(elem_rxn.get('id'))
            else:
                species_nondiff.append(elem_rxn.get('id'))
    return species_diff ,species_nondiff

try:
    args = ARGS.split(" ")
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

#### to update IC values from a simulation with stimulation, specify -start beg end
#where end is just prior to stimulation, and beg is a few ms prior to stimulation
#if there are two regions, e.g. dend and soma, this code  may not work
dendname="dend" #from morph file
spinehead="head" #from morph file
submembname=dendname+'sub'
decimals=4 #how many decimals to use in the updated IC file values

#load the data 
params=h5utilsV2.parse_args(args,do_exit)#h5utilsV2.argparse(args)
if len(args[1]):
    all_species=params.mol
else:
    all_species=None

tot_species,weight,sub_species=h5utilsV2.get_tot(params)

og=nrdh5_group(params.fileroot,params.par,tot_species) 
for fnum,ftuple in enumerate(og.ftuples):
    data=nrdh5_output(ftuple)
    data.rows_and_columns(all_species,params.start,params.end)
    data.molecule_population()
    #print(data.data['model']['grid'][:])
    if data.maxvols>1:
        data.region_structures(dendname,submembname,spinehead)#,stimspine) #stimspine is optional
        data.average_over_voxels()
    data.total_subspecies(tot_species,sub_species,params.start,weights=weight)

#mol concentration for both spine and submem
mole_conc_ic={M:{} for M in data.molecules}
for mol in data.molecules:
    mole_conc_ic[mol]['general']=np.mean(data.OverallMean[mol][0,data.sstart[mol]:data.ssend[mol]])# 0 indicates first trial
   
if data.maxvols>1:
    mylist=['region','struct']
    for ii, regions_params in enumerate (mylist):
        for mol in data.molecules:
            for jj,key in enumerate(data.output_labels[regions_params][mol].split()):
                region=key.split('_')[-1]
                mole_conc_ic[mol][region]=np.mean(data.means[regions_params][mol][0,data.sstart[mol]:data.ssend[mol],jj])# 0 indicates first trial

#read in initial conditions file
IC_filename=params.IC 
print(IC_filename)
tree=ET.parse(IC_filename+'.xml')
root=tree.getroot()
tags=list(np.unique([rt.tag for rt in root]))

IC_types={t:{} for t in tags} #dictionary to map region names in IC file to region names in mole_conc_ic
for elem_ic in root:
    if elem_ic.tag=='ConcentrationSet':
        IC_types[elem_ic.tag][elem_ic.get('region')]=elem_ic.get('region')
    if elem_ic.tag=='SurfaceDensitySet':
        IC_types[elem_ic.tag][elem_ic.get('region')]=elem_ic.get('region')+'sub'

## read in reaction file
Rxn_filename=params.Rxn
## determine which species diffuse
diffuse_species,non_diffuse_mol=DiffuseSpecies(Rxn_filename)

########### Initialize some lists and dictionaries 
IC_molecules=[]
all_species=AllSpecies(root)
region_conc={m:{} for sp in all_species.values() for m in sp }
all_species=data.molecules #Use this? This could be more molecules than from AllSpecies, but it includes molecules not in IC file
region_conc={sp:{} for sp in all_species }

vol={r:data.region_dict[r]['vol'] for r in data.region_dict.keys()}
for r in data.region_struct_dict.keys():
    vol[r]=data.region_struct_dict[r]['vol']
########### Update concentration in IC file #####################
for rt in root:
    subtree=ET.ElementTree(rt)
    if len(rt.attrib):
        region=IC_types[rt.tag][rt.attrib['region']]
    else:
        region='general'
    if rt.tag=='ConcentrationSet' and region=='general':
        re_do_root=rt
        elems=subtree.findall('.//NanoMolarity')
        for e in elems:
            mol=e.attrib['specieID']
            if mol in diffuse_species:
                e.attrib['value']=str(np.round(mole_conc_ic[mol][region],decimals)) 
                print(mol,'updated elems',region,e.attrib)
                IC_molecules.append(mol)
                region_conc[mol][region]=mole_conc_ic[mol][region]
            else:
                print(mol,'NOT updating ', region,' for non-diffusible', e.attrib) #get rid of this, since all are updated on second pass?
    else:
        if rt.tag=='ConcentrationSet':
            elems=subtree.findall('.//NanoMolarity')
            height=1
        elif rt.tag=='SurfaceDensitySet':
            elems=tree.findall('.//PicoSD')
            height=data.region_struct_dict[region]['depth']
        for e in elems:
            mol=e.attrib['specieID']
            e.attrib['value']=str(np.round(mole_conc_ic[mol][region]*height,decimals))
            IC_molecules.append(mol)
            region_conc[mol][region]=mole_conc_ic[mol][region]
            print(mol,' SD or Region - updated elems',rt.tag,region,e.attrib,region_conc[mol])
IC_molecules=list(np.unique(IC_molecules))
### Some molecules are not diffusible, but are not specified in region or SurfaceDensitySet
### Those should be specified using general concentration set #######
### Some molecules are not diffusible, and specificed in region or surfaceDensitySet, and should also be non-zero in cytosol
print(u'◢◤◢◤◢◤◢◤ second pass ◢◤◢◤◢◤◢◤')
subtree=ET.ElementTree(re_do_root) #re_do_root gets subtree for general concentration set only
region='dendcyt' ###update general with this value for non-diffusible molecules to allow a dendcyt value different from dendsub
elems=subtree.findall('.//NanoMolarity')
for e in elems:
    mol=e.attrib['specieID']
    if mol not in IC_molecules:
        e.attrib['value']=str(np.round(mole_conc_ic[mol]['general'],decimals))
        print('Updated conc set for non-diffusible molecule (on 2nd pass)',region,mol,e.attrib)
        IC_molecules.append(mol)
        region_conc[mol]['general']=mole_conc_ic[mol]['general']
    elif  mol in non_diffuse_mol and e.attrib['value']!=0: # if not explicitly zero, update general to add molecules to cytosol
        e.attrib['value']=str(np.round(mole_conc_ic[mol][region],decimals))
        region_conc[mol]['dend'+'cyt']=mole_conc_ic[mol][region]
        print('>>>>> Updated (on 2nd pass) general conc set, using ',region,'with region overrides',region,mol,', IC=',e.attrib,region_conc[mol])
    '''else: #this stuff should be commented about after debugging is finished
        print('&&&&&&&&&&&& NOT updating!',mol,region,', total (moles)=',region_conc[mol],', IC (conc)=',e.attrib)
        if  mol in non_diffuse_mol:
            print('    mol in non_diffuse_mol and conc zero')
        elif e.attrib['value']!=0:
            print('      e.attrib is non-zero in IC, molecule diffusible')
        else:
            print('      e.attrib is zero in IC and molecule diffusible')
    '''    
########### check if all mol are present in the IC files ##############
IC_file_mols=set(IC_molecules)
h5_mols=set(mole_conc_ic.keys())
both_mols=IC_file_mols & h5_mols

if len(h5_mols-both_mols)>0:
    for mol in h5_mols-both_mols:
        if np.max(list(mole_conc_ic[mol].values()))>0:
            print ('******* ', mol,'missing from IC file, conc non-zero in h5:', mole_conc_ic[mol])
if len(IC_file_mols-both_mols)>0:
    print('************** molecules',IC_file_mols-both_mols, 'in conc file but not in IC file')

############ Calculate total concentration specified in IC file and compare to h5 ##############
total_conc={'general':{mol:0 for mol in tot_species}, 'Sum_specific':{mol:0 for mol in tot_species}}
for imol,mol in enumerate(tot_species):
    mol_set=[];total_subsp={}
    #### first set up arrays of all species (sub_species) that are a form of the molecule
    if mol in sub_species.keys():
        mol_set=sub_species[mol]
    else:
        mol_set=[subsp for subsp in IC_molecules if mol in subsp] #h5utils.decode(self.data['model']['output'][outset]['species'][:])
    for subspecie in mol_set:
        total_conc['general'][mol]+=mole_conc_ic[subspecie]['general']
        if list(region_conc[subspecie].keys())==['general']:
            total_subsp[subspecie]=region_conc[subspecie]['general']
            print (subspecie,'one general IC spec:', region_conc[subspecie].keys(),', tot conc subsp', total_subsp[subspecie])
        elif 'general' not in region_conc[subspecie].keys():
            total_subsp[subspecie]=np.sum([conc*vol[reg] for reg,conc in region_conc[subspecie].items()])/data.TotVol
            print (subspecie,'only region specific IC specs:', region_conc[subspecie].keys(),', tot conc subsp',total_subsp[subspecie])
        elif len(region_conc[subspecie].keys()):
            total_subsp[subspecie]=np.sum([conc*vol[reg] for reg,conc in region_conc[subspecie].items() if reg != 'general'])
            vol_remain=data.TotVol-np.sum([vol[reg] for reg in region_conc[subspecie].keys() if reg != 'general'])
            total_subsp[subspecie]+=region_conc[subspecie]['general']*vol_remain #add in molecules in additional regions
            total_subsp[subspecie]/=data.TotVol
            print (subspecie,'region specific AND general IC specs', region_conc[subspecie].keys(),', tot conc subsp',total_subsp[subspecie])
        else:
            print(subspecie,' !!!! no keys !!!!')
        total_conc['Sum_specific'][mol]+=total_subsp[subspecie]
        #print('>>> running total for:', subspecie,total_conc['Sum_specific'][mol])
    print('######## Total for',mol,', IC general=',round(total_conc['general'][mol],4),
          ', IC specific=',round(total_conc['Sum_specific'][mol],4),
          'total in h5 (beg,end)=',round(data.total_trace['Overall'][mol][0][0],4),round(data.total_trace['Overall'][mol][0][-1],4))
        
############# write the new IC file ###################
outfile=IC_filename+'h5_update.xml'
with open(outfile, 'wb') as out:
    out.write(ET.tostring(root))            
 
