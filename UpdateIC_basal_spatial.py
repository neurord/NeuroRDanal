
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
    species_diff=[]
    tree_rxn=ET.parse(Rxn_filename+'.xml')
    root_rxn=tree_rxn.getroot()
    for elem_rxn in root_rxn:
        if elem_rxn.tag== 'Specie':
            if float(elem_rxn.get('kdiff'))>0:
                species_diff.append(elem_rxn.get('id'))
    return species_diff 

try:
    args = ARGS.split(" ")
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

dendname="dend" #from morph file
spinehead="head" #from morph file
submembname=dendname+'sub' 

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
    mole_conc_ic[mol]['general']=np.mean(np.mean(data.OverallMean[mol]))
if data.maxvols>1:
    mylist=['region','struct']
    for ii, regions_params in enumerate (mylist):
        for mol in data.molecules:
            for jj,key in enumerate(data.output_labels[regions_params][mol].split()):
                region=key.split('_')[-1]
                mole_conc_ic[mol][region]=np.mean(np.mean(data.means[regions_params][mol][:,-10:,jj],axis=1))

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
diffuse_species=DiffuseSpecies(Rxn_filename)
all_species=AllSpecies(root) 

########### Initialize some lists and dictionaries 
IC_molecules=[];non_diffuse_mol=[]
total={m:{} for sp in all_species.values() for m in sp }

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
                e.attrib['value']=str(np.round(mole_conc_ic[mol][region],3))
                print('updated elems',region,mol,e.attrib)
                IC_molecules.append(mol)
                total[mol][region]=mole_conc_ic[mol][region]*data.TotVol
            else:
                print('NOT updating ', region,' for non-diffusible', mol,e.attrib)
                non_diffuse_mol.append(mol)
    else:
        if rt.tag=='ConcentrationSet':
            elems=subtree.findall('.//NanoMolarity')
            height=1
            volume=data.region_dict[region]['vol']
        elif rt.tag=='SurfaceDensitySet':
            elems=tree.findall('.//PicoSD')
            height=data.region_struct_dict[region]['depth']
            if region != 'general':
                volume=data.region_struct_dict[region]['vol']
        for e in elems:
            mol=e.attrib['specieID']
            e.attrib['value']=str(np.round(mole_conc_ic[mol][region]*height,3))
            print(' SD or Region - updated elems',rt.tag,region,e.attrib)
            IC_molecules.append(mol)
            total[mol][region]=mole_conc_ic[mol][region]*volume
            print('**********',mol,region,'vol=',round(volume,5),'total=',total[mol])

### Some molecules are not diffusible, but are not specified in region or SurfaceDensitySet
### Those should be specified using general concentration set #######
### Some molecules are not diffusible, and specificed in region or surfaceDensitySet, and should also be non-zero in cytosol
subtree=ET.ElementTree(re_do_root)
region='general'
elems=subtree.findall('.//NanoMolarity')
for e in elems:
    mol=e.attrib['specieID']
    if mol not in IC_molecules:
        e.attrib['value']=str(np.round(mole_conc_ic[mol]['general'],3))
        print('Updated (on 2nd pass)',region,mol,e.attrib)
        IC_molecules.append(mol)
        total[mol]['general']=mole_conc_ic[mol]['general']*data.TotVol
    elif  mol in non_diffuse_mol and e.attrib['value']!=0: # if not explicitly zero, update with dendcyt value
        e.attrib['value']=str(np.round(mole_conc_ic[mol][region],3))
        print('Updated (on 2nd pass) non-zero dend cyt',region,mol,e.attrib)
        volume=data.region_struct_dict['dend'+'cyt']['vol'] #FIXME  'dend' is a region
        total[mol]['dend'+'cyt']=mole_conc_ic[mol][region]*volume
        print('**********',mol,region,'vol=',volume,'total=',total[mol])

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
        if 'general' in total[subspecie].keys():
            print (subspecie,'one total key?', total[subspecie].keys())
            total_subsp[subspecie]=total[subspecie]['general']/data.TotVol
        else:
            total_subsp[subspecie]=np.sum([t for m,t in total[subspecie].items() if m != 'general'])/data.TotVol
            if len(total[subspecie].keys()):
                print (subspecie,'many total keys?', total[subspecie].keys(),'sum',total_subsp[subspecie])
            else:
                print(subspecie,'no keys')
        total_conc['Sum_specific'][mol]+=total_subsp[subspecie]
        #print('>>> running total for:', subspecie,total_conc['Sum_specific'][mol])
    print('######## Total for',mol,', IC general=',round(total_conc['general'][mol],4),
          ', IC specific=',round(total_conc['Sum_specific'][mol],4),
          'total in h5 (beg,end)=',round(data.total_trace['Overall'][mol][0][0],4),round(data.total_trace['Overall'][mol][0][-1],4))
        
############# write the new IC file ###################
outfile=IC_filename+'h5_update.xml'
with open(outfile, 'wb') as out:
    out.write(ET.tostring(root))            
 
