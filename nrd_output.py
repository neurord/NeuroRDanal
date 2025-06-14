# -*- coding: utf-8 -*-
"""
Created on Sun Apr 26 14:41:21 2020

@author: kblackw1
"""
from tokenize import group
import numpy as np
import h5py as h5
import h5utilsV2 as h5utils

Avogadro=6.02214179e14 #molecules per nanoMole
mol_per_nM_u3=Avogadro*1e-15 #0.6022 = PUVC, 1e-15 is liters per micron cubed
mol_per_pM_u2=Avogadro*1e-15

class nrdh5_output(object):
    def __init__(self,ftuple):
        self.fname=ftuple[0]   
        self.data = h5.File(self.fname,"r") 
        self.parval=ftuple[1] 
        self.TotVol=self.data['model']['grid'][:]['volume'].sum() 
        self.maxvols=len(self.data['model']['grid'])  
        self.trials=[a for a in self.data.keys() if 'trial' in a] 
        self.seeds=[self.data[trial].attrs['simulation_seed'] for trial in self.trials]      
        self.outputsets=self.data[self.trials[0]]['output'].keys()
        self.head=None
        self.spinelist=[]
        self.dsm_vox=None
        self.spatial_dict=None
        print('file:',self.fname,'parameters:',self.parval,'volume',self.TotVol,'voxels',self.maxvols,'trials',self.trials)
        
    def region_structures(self,dendname,submembname,spinehead,stimspine=None): #spinehead may not need to be parameter
        structType=self.data['model']['grid'][:]['type']
        self.region_list,self.region_dict,self.region_struct_dict=h5utils.subvol_list(structType,self.data['model'])
        #Replace the following with test for whether there is more than one "group"
        if spinehead in self.region_dict.keys():
            self.head=spinehead
            self.spinelist,self.spinevox=h5utils.multi_spines(self.data['model'])
            if stimspine is not None:
                self.stim_spine_index=[self.spinelist.index(stimsp) for stimsp in stimspine]        #
            else:
                self.stim_spine_index=-1
            #create "group" dictionary-spinelist & spinevox, with one more for nonspine voxels
        self.dsm_name=dendname+submembname
        if self.dsm_name in self.region_struct_dict.keys():
            self.dsm_vox=self.region_struct_dict[self.dsm_name]
            self.dsm_index=list(self.region_struct_dict.keys()).index(self.dsm_name)
            
    def spatial_structures(self,bins,dendname):
        self.spatial_dict=h5utils.spatial_average(bins,dendname,self.data['model']['grid'])
        vol=[x['vol'] for x in self.spatial_dict.values()]
        if any(v==0 for v in vol):
            print ("**********Too many bins for spatial average, dropping these bins ****************")
            delete_keys=[] #identify keys, but don't delete
            for key,vox_vol in self.spatial_dict.items():
                if len(vox_vol['vox'])==0:
                    print(key,vox_vol)
                    delete_keys.append(key)
            for key in delete_keys: #now, delete keys of ordered dict while looping over delete_keys
                del self.spatial_dict[key]
        print('final dict:',self.spatial_dict)
#
    def rows_and_columns(self,molecules,start_end,endtime=None):
        #which columns in the data set contain counts of molecules of interest
        self.all_molecules=h5utils.decode(self.data['model']['species'][:])
        if molecules is not None:
             self.molecules = [mol for mol in molecules if mol in self.all_molecules]
             missing_molecules=[mol for mol in molecules if mol not in self.all_molecules]
             if len(missing_molecules):
                 print("** MOLECULE", missing_molecules, " DOES NOT EXIST !!!!!!!!!!!!!")
                 print('choose from',self.all_molecules)
        else:
            self.molecules=self.all_molecules 
        out_location,dt,rows=h5utils.get_mol_info(self.data, self.molecules,endtime)
        self.out_location=out_location
        self.dt=dt
        self.rows=rows
        #Which "rows" should be used for baseline value. If different for each file then will have problems later
        sstart,ssend=h5utils.sstart_end(self.molecules,self.out_location,self.dt,self.rows,start_end)
        self.sstart=sstart
        self.ssend=ssend
        
    def molecule_population(self):
        self.counts={}
        self.time={}
        self.OverallMean={}
        for imol,molecule in enumerate(self.molecules):
              counts,time=h5utils.get_mol_pop(self,molecule)
              self.counts[molecule]=counts
              self.time[molecule]=time
              self.OverallMean[molecule]=np.sum(counts[:,:,:],axis=2)/(self.TotVol*mol_per_nM_u3)
        
    def average_over_voxels(self):
        #calculate various averages.  These will be used for plotting and output
        self.output_labels={'struct':{},'region':{}}
        self.means={'struct':{},'region':{}}
        if self.spinelist:
            self.output_labels['spines']={};self.means['spines']={}
        if self.spatial_dict:
            self.output_labels['space']={}; self.means['space']={}
        for imol,mol in enumerate(self.molecules):
            #calculate region-structure means, such as dendrite submembrane and dendrite cytosol
            #dimensions: number of trials x time samples x number of regions
            #labels hold the name of the region / structure, such as dendrite submembrane and dendrite cytosol
            self.output_labels['struct'][mol],self.means['struct'][mol]=h5utils.region_means_dict(self,mol,self.region_struct_dict)            
            #regions, such as dendrite, soma, spine head
            self.output_labels['region'][mol],self.means['region'][mol]=h5utils.region_means_dict(self,mol,self.region_dict)
            #if more than one spine, calculate individual spine means
            if self.spinelist:
                self.output_labels['spines'][mol],self.means['spines'][mol]=h5utils.region_means_dict(self,mol,self.spinevox)
            if self.spatial_dict:
                self.output_labels['space'][mol],self.means['space'][mol]=h5utils.region_means_dict(self,mol,self.spatial_dict)

    def write_average(self,savedir): #one file per molecule, mean and std of all regions
        import os ########### NEED TO PUT SPATIAL STUFF IN SEPARATE FILE?
        for mol_set in [self.molecules,self.tot_species.keys()]:
            for mol in mol_set:
                outfilename=savedir+os.path.splitext(os.path.basename(self.fname))[0]+mol+'_avg.txt'
                col_name='_'.join([str(q) for q in self.parval])
                mean_header=mol+col_name+'_All ' #first non-time column of header
                if mol in self.molecules: 
                    output_means=np.mean(self.OverallMean[mol],axis=0) #average over trials
                    output_std=np.std(self.OverallMean[mol],axis=0)
                    time=self.time[mol]
                    if self.maxvols>1:
                        for mean_dict in self.means.values():
                            output_means=np.column_stack((output_means,np.mean(mean_dict[mol],axis=0)))
                            output_std=np.column_stack((output_std,np.std(mean_dict[mol],axis=0)))
                        for label in self.output_labels.values():
                            mean_header=mean_header+label[mol]
                elif mol in self.tot_species.keys():
                    output_means=np.mean(self.total_trace['Overall'][mol],axis=0)
                    output_std=np.std(self.total_trace['Overall'][mol],axis=0)
                    time=np.linspace(0,self.endtime[mol], len(output_means))
                    if self.maxvols>1:
                        for reg in list(self.total_trace.keys())[1:]:
                               output_means=np.column_stack((output_means,np.mean(self.total_trace[reg][mol],axis=0)))
                               output_std=np.column_stack((output_std,np.std(self.total_trace[reg][mol],axis=0)))
                               mean_header=mean_header+' '+mol+col_name+'_'+reg
                std_header='_std '.join(mean_header.split())
                output_header='Time '+ mean_header+' '+std_header+'_std\n'
                #print(outfilename,output_header)
                np.savetxt(outfilename, np.column_stack((time,output_means,output_std)), fmt='%.4f', delimiter=' ', header=output_header)            
        
    #
    #FOR SIGNATURE: use this code, but specify which molecules to total for the signature
    #add option to baseline subtract - separate function?
    def total_subspecies(self,tot_species,sub_species,start_end,outset='__main__',weights={}):
        print_trial=-1 #set to 0 for debugging output
        self.tot_species={t:[] for t in tot_species}
        self.endtime={t:[] for t in tot_species}
        self.total_trace={'Overall':{}}
        if self.maxvols > 1:
            for xdict in [self.region_dict, self.region_struct_dict]:
                for reg in xdict.keys():
                    self.total_trace[reg]={}
        if self.spinelist:
            for sp in self.spinelist:
                self.total_trace[sp]={}
        for imol,mol in enumerate(tot_species):
            mol_set=[];tot_wt=0
            #### first set up arrays of all species (sub_species) that are a form of the molecule
            if mol in sub_species.keys():
                mol_set=sub_species[mol]
            else:
                mol_set=[subsp for subsp in self.all_molecules if mol in subsp] #h5utils.decode(self.data['model']['output'][outset]['species'][:])
            self.tot_species[mol]=mol_set
            out_location,dt,rows=h5utils.get_mol_info(self.data, mol_set)
            sstart,ssend=h5utils.sstart_end(mol_set,out_location,dt,rows,start_end)
            if len(np.unique(list(dt.values())))>1:
                print('WARNING, subspecies have',mol,' have different dt',dt,'using largest and interpolating',rows)
            dt_index=np.argmax(list(dt.values()))
            self.dt[mol]=list(dt.values())[dt_index] #use largest dt, and then subsample subspecies with smaller dt
            self.endtime[mol]=(rows[dt_index]-1)*self.dt[mol] #FIXME, use rows corresponding to correct dt, possible rows[1]?
            self.total_trace['Overall'][mol]=np.zeros((len(self.trials),rows[dt_index]))
            if self.maxvols > 1:
                for xdict in [self.region_dict, self.region_struct_dict]:
                    for reg, reg_dict in xdict.items():
                        self.total_trace[reg][mol]=np.zeros((len(self.trials),rows[dt_index]))       
            #if self.dsm_vox:
                #self.total_trace['dsm'][mol]=np.zeros((len(self.trials),rows[dt_index]))
            if self.spinelist: 
                for sp in self.spinelist:
                    self.total_trace[sp][mol]=np.zeros((len(self.trials),rows[dt_index]))
            #### second, find molecule index of the sub_species and total them
            for subspecie in mol_set:
                wt=weights[subspecie] if subspecie in weights.keys() else 1
                for trialnum,trial in enumerate(self.trials):
                    if subspecie not in self.molecules:
                        for outname,outset in out_location[subspecie]['location'].items():
                            ###### trials,times,voxels ########
                            mol_pop=self.data[trial]['output'][outname]['population'][:,:,outset['mol_index']] # FIXME won't work if > 1 outset
                    else:
                        mol_pop=self.counts[subspecie][trialnum]
                    pop_sum=np.sum(mol_pop,axis=1) #sum over voxels
                    #print('tot_species',mol,subspecie,np.shape(mol_pop),np.shape(self.total_trace['Overall'][mol][trialnum]))
                    if dt[subspecie] == self.dt[mol]:
                        self.sstart[mol]=sstart[subspecie]
                        self.ssend[mol]=ssend[subspecie]
                    else:
                        #interpolate
                        target_t=np.arange(len(self.total_trace['Overall'][mol][trialnum]))*self.dt[mol]
                        subspecie_t=np.arange(len(mol_pop))*dt[subspecie]
                        if trialnum==print_trial:
                            print('mol=',mol,subspecie,'trial',trial,'len tot=',len(target_t),'len subsp=',len(subspecie_t),'tot time=',target_t[-1],'subsp time=', subspecie_t[-1])
                        pop_sum=np.interp(target_t,subspecie_t, pop_sum)
                        if trialnum==print_trial:
                            print('after interp, len interp=',len(pop_sum), 'len tot=', len(self.total_trace['Overall'][mol][trialnum]))
                    self.total_trace['Overall'][mol][trialnum]+=wt*pop_sum/self.TotVol/mol_per_nM_u3
                    #then total sub_species in submembrane and spine head, if they exist  
                    if self.maxvols > 1:
                        for ydict in [self.region_dict, self.region_struct_dict]:
                            for reg,val in ydict.items():   
                                pop_sum=np.sum(mol_pop[:,val['vox'] ],axis=1)
                                if dt[subspecie] != self.dt[mol]: #assume that subspecie_t is same for all regions and structures 
                                    pop_sum=np.interp(target_t,subspecie_t, pop_sum)  
                                if reg == self.dsm_name:
                                    self.total_trace[reg][mol][trialnum]+=wt*(pop_sum/val['vol'])*val['depth']/mol_per_pM_u2
                                else:
                                    self.total_trace[reg][mol][trialnum]+=wt*(pop_sum/val['vol'])/mol_per_nM_u3
                    if self.spinelist:
                        for spnum,spname in enumerate(self.spinevox):
                            pop_sum=np.sum(mol_pop[:,self.spinevox[spname]['vox']],axis=1)
                            if dt[subspecie] != self.dt[mol]: 
                                if trialnum==print_trial:
                                    print(spname,subspecie,' dt subsp=',dt[subspecie],'dt mol=',self.dt[mol])
                                pop_sum=np.interp(target_t,subspecie_t, pop_sum)
                            self.total_trace[spname][mol][trialnum]+=wt*pop_sum/self.spinevox[spname]['vol']/mol_per_nM_u3
            outputline='####  '+str(self.parval)+'  >>>'+mol+' TOTAL: '+str(np.round(self.total_trace['Overall'][mol][0,0],3))+' nM'
            if self.spinelist:
                for spname in self.spinevox:
                    outputline +=', '+spname+': '+str(np.round(self.total_trace[spname][mol][0,0],3))+' nM'
            if self.dsm_vox:
                outputline +=',  dsm: '+str(np.round(self.total_trace[self.dsm_name][mol][0,0],3))+' pSD'
            print(outputline)
  
    def print_head_stats(self):
        for imol,mol in enumerate(self.molecules):
            outputline=mol.rjust(14)
            if self.spinelist:
                if len(self.spinelist)>1 and self.stim_spine_index!=-1: 
                    headmean=[np.mean(np.mean(self.means['spines'][mol][:,self.sstart[mol]:self.ssend[mol],sp],axis=1),axis=0) for sp in self.stim_spine_index]
                    headmax=[np.mean(np.max(self.means['spines'][mol][:,self.ssend[mol]:,sp],axis=1),axis=0) for sp in self.stim_spine_index]
                else:
                    headmean=[np.mean(np.mean(self.means['region'][mol][:,self.sstart[mol]:self.ssend[mol],self.head],axis=1),axis=0)]
                    headmax=[np.mean(np.max(self.means['region'][mol][:,self.ssend[mol]:,self.head],axis=1),axis=0)]
                outputline+="   head ss: "+' '.join([str(np.round(h,3)) for h in headmean])+' pk: '+' '.join([str(np.round(h,3)) for h in headmax])
            if self.dsm_vox:
                dsm_mean=np.mean(np.mean(self.means['struct'][mol][:,self.sstart[mol]:self.ssend[mol],self.dsm_index],axis=1),axis=0)
                dsm_max=np.mean(np.max(self.means['struct'][mol][:,self.ssend[mol]:,self.dsm_index],axis=1),axis=0)
                outputline+="   dend sm: %8.4f pk %8.4f" %(dsm_mean,dsm_max)
            print(outputline)
