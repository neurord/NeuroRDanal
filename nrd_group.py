# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:59:11 2020

@author: kblackw1
"""
import numpy as np
import os
import h5utilsV2 as h5utils

ms_to_sec=1000
class nrdh5_group(object): 
    def __init__(self,fileroot,parameters,tot_species=[],savedir=None):
        self.ftuples,self.parlist,self.params=h5utils.create_filenames(fileroot,parameters)
        self.par_keys=[tuple([str(r)+str(p) for r,p in zip(self.params,par[1])]) for par in self.ftuples]
        self.file_set_conc={k:{'Overall':{}} for k in self.par_keys}
        self.time_set={k:{} for k in self.par_keys}
        self.means={k:{reg:{} for reg in ['space']} for k in self.par_keys} 
        if savedir:
            self.savedir=savedir+'/'
        elif os.path.dirname(fileroot):
            self.savedir=os.path.dirname(fileroot)+'/'
        else:
            self.savedir=''

        if len(tot_species):
            self.tot_species=tot_species
            self.endtime={k:{sp:[] for sp in self.tot_species} for k in self.par_keys}
            self.file_set_tot={k:{'Overall':{sp:[] for sp in self.tot_species}} for k in self.par_keys}
        else:
            self.tot_species=[]
            self.file_set_tot={k:{'Overall':{}} for k in self.par_keys}

    def conc_arrays(self,data):
        self.molecules=data.molecules
        self.trials=data.trials
        #These are overwritten with each data file, and must be the same for each data file
        self.sstart={mol:data.sstart[mol] for mol in data.molecules}
        self.ssend={mol:data.ssend[mol] for mol in data.molecules}
        self.dt={mol:data.dt[mol] for mol in data.molecules}
        self.all_regions=['Overall']
        parval=tuple([str(r)+str(p) for r,p in zip(self.params,data.parval)])
        if data.maxvols>1:
            for reg_dict in [data.region_dict,data.region_struct_dict]:
                for regnum,reg in enumerate(reg_dict.keys()):
                    if reg not in self.file_set_conc[parval].keys():
                        self.file_set_conc[parval][reg]={}
                        self.all_regions.append(reg)
            if data.spinelist:
                for regnum,reg in enumerate(data.spinelist):
                    if reg not in self.file_set_conc[parval].keys():
                        self.file_set_conc[parval][reg]={}
                        self.all_regions.append(reg)
        for imol,molecule in enumerate(data.molecules):
            self.time_set[parval][molecule]=data.time[molecule]
            self.file_set_conc[parval]['Overall'][molecule]=data.OverallMean[molecule]
            if data.maxvols>1:
                for reg_type,reg_dict in zip(['region','struct'],[data.region_dict,data.region_struct_dict]):
                    for regnum,reg in enumerate(reg_dict.keys()):
                        self.file_set_conc[parval][reg][molecule]=data.means[reg_type][molecule][:,:,regnum]
                if data.spinelist:
                     for regnum,reg in enumerate(data.spinelist):
                        self.file_set_conc[parval][reg][molecule]=data.means['spines'][molecule][:,:,regnum]
                if data.spatial_dict:
                    self.means[parval]['space'][molecule]=data.means['space'][molecule]
            else:
                self.spatial_data=None
        if len(self.tot_species):
            for reg in data.total_trace.keys(): #is order of regions in total_trace same as file_set_conc , file_set_tot, and all_regions?
                if reg not in self.file_set_tot[parval].keys():
                    self.file_set_tot[parval][reg]={sp:[] for sp in self.tot_species} 
                for imol,sp in enumerate(self.file_set_tot[parval][reg].keys()):
                    self.file_set_tot[parval][reg][sp]=data.total_trace[reg][sp][:,:]
                    self.sstart[sp]=data.sstart[sp]
                    self.ssend[sp]=data.ssend[sp]
                    self.dt[sp]=data.dt[sp]
                    self.endtime[parval][sp]=data.endtime[sp]
        
    def trace_features(self,window_size,lo_thresh_factor=0.2,hi_thresh_factor=0.8,std_factor=2,numstim=1,end_baseline_start=0,filt_length=5,aucend=None,iti=None):
        import operator
        self.feature_list=['baseline','basestd','peakval','peaktime','amplitude','duration','slope','minval','auc','auc_thresh','start_plateau', 'end_plateau']
        self.feature_scale={'baseline':1,'basestd':1,'peakval':1,'peaktime':ms_to_sec,'amplitude':1,'duration':ms_to_sec,'slope':1,'minval':1,'auc':ms_to_sec,'auc_thresh':1,'start_plateau':ms_to_sec, 'end_plateau':ms_to_sec}
        self.feature_dict={feat:np.zeros((len(self.molecules)+len(self.tot_species),len(self.ftuples), len(self.all_regions),len(self.trials))) for feat in self.feature_list}
        self.mean_feature={}
        self.std_feature={}

        molecules=self.molecules #use individual molecules for self.file_set_conc and self.means
        ii = 0
        for parnum,par in enumerate(self.par_keys):
            #print('*************',region,regnum,traces.keys())
            self.tf_base_function(window_size, molecules, self.file_set_conc[par], ii, parnum, par)

        molecules=self.tot_species
        ii=1
        for parnum,par in enumerate(self.par_keys):
            #print('****** nrd_group line 84 *******',region,regnum,traces.keys())
            self.tf_base_function(window_size, molecules, self.file_set_tot[par], ii, parnum, par) 

        for feat in self.feature_dict.keys():
            self.mean_feature[feat]=np.nanmean(self.feature_dict[feat],axis=-1)
            self.std_feature[feat]=np.nanstd(self.feature_dict[feat],axis=-1)

    def tf_base_function(self,window_size, molecules, traces, ii, parnum, par, lo_thresh_factor=0.2,hi_thresh_factor=0.8,std_factor=2,numstim=1,end_baseline_start=0,filt_length=5,aucend=None,iti=None):
        import operator
        for regnum, region in enumerate(traces.keys()): 
            for jmol,mol in enumerate(molecules):
                #print(par,mol,np.shape(traces[par][mol]))
                imol=jmol+ii*len(self.molecules)
                window=int(window_size/self.dt[mol]) #FIXME self.dt[mol],sstart[mol],ssend[mol]
                if window==0:
                    print('trace_features_loop, window size too small for',par,mol,'dt',self.dt[mol],'window size',window_size,'using',self.dt[mol])
                    window=1
                self.feature_dict['baseline'][imol,parnum,regnum,:]= np.mean(traces[region][mol][:,self.sstart[mol]:self.ssend[mol]],axis=1)
                self.feature_dict['basestd'][imol,parnum,regnum,:]=np.std(traces[region][mol][:,self.sstart[mol]:self.ssend[mol]],axis=1)
                peakpt=np.argmax(traces[region][mol][:,self.ssend[mol]:],axis=1)+self.ssend[mol]
                self.feature_dict['peaktime'][imol,parnum,regnum,:]=peakpt*self.dt[mol]/self.feature_scale['peaktime']
                self.feature_dict['peakval'][imol,parnum,regnum,:]=[np.mean(traces[region][mol][i,p-window:p+window]) 
                                                             for i,p in enumerate(peakpt)]
                lowpt=np.argmin(traces[region][mol][:,self.ssend[mol]:],axis=1)+self.ssend[mol]
                #DEBUGGING: peak time and peakval from mean of trace
                p=np.argmax(np.mean(traces[region][mol][:,self.ssend[mol]:],axis=0))+self.ssend[mol]
                pt=p*self.dt[mol]
                pval=np.mean(np.mean(traces[region][mol],axis=0)[p-window:p+window])
                #print('peaktime {} & peakval {} from mean trace'.format(pt,pval))
                #print('lowpt',mol,lowpt,'window',window,'lowval mean',np.mean([traces[par][mol][i,p-window:p+window] for i,p in enumerate(lowpt)]))
                #end DEBUGGING
                self.feature_dict['minval'][imol,parnum,regnum,:]=[np.mean(traces[region][mol][i,p-window:p+window]) 
                                                            for i,p in enumerate(lowpt)]
                self.feature_dict['amplitude'][imol,parnum,regnum,:]=self.feature_dict['peakval'][imol,parnum,regnum,:]-self.feature_dict['baseline'][imol,parnum,regnum,:]
                ####################

                self.slope(traces[region][mol], imol, parnum, mol, regnum)
                self.plateau_duration(traces[region][mol], imol, parnum, mol, peakpt, filt_length, regnum)
                self.auc(traces[region][mol], imol, parnum, par, mol, numstim, std_factor, aucend, iti, end_baseline_start, filt_length, regnum)
                
    def slope(self, traces, imol, parnum, mol, regnum,lo_thresh_factor=0.2,hi_thresh_factor=0.8):           
        
        import operator        
        #FIND SLOPE OF INCREASE - Use thresholds defined by lo_thresh, and hi_thresh, e.g. 20 and 80%
        lo_thresh=lo_thresh_factor*(self.feature_dict['amplitude'][imol,parnum,regnum,:])+self.feature_dict['baseline'][imol,parnum,regnum,:] #get the 5% above the max value
        hi_thresh=hi_thresh_factor*(self.feature_dict['amplitude'][imol,parnum,regnum,:])+self.feature_dict['baseline'][imol,parnum,regnum,:]
        self.ssend_list=[self.ssend[mol] for t in self.trials]
        #
        begin_slope=exceeds_thresh_points(traces, self.ssend_list, lo_thresh,operator.gt)
        end_slope=exceeds_thresh_points(traces, begin_slope,hi_thresh,operator.gt)
        for i,(end,beg) in enumerate(zip(end_slope,begin_slope)):
            if end-beg>1 and ~np.isnan(end) and ~np.isnan(beg): #FIX - check for end=beg - making zero slope, or even end_slope.beginslope=1
                self.feature_dict['slope'][imol,parnum,regnum,i]=(self.feature_dict['peakval'][imol,parnum,regnum,i]-
                                                        self.feature_dict['baseline'][imol,parnum,regnum,i])/((end-beg)*self.dt[mol])
            else:
                self.feature_dict['slope'][imol,parnum,regnum,:]=np.nan                    
        ####################
                        
    def plateau_duration(self, traces, imol, parnum, mol, peakpt, filt_length, regnum): 

        import operator          
        # FIND PLATEAU DURATION - USE thresholds defined by midpoints, and two different time periods
        #could also use thresholds defined by lo_thresh or hi_thresh
        midpoints=0.5*(self.feature_dict['amplitude'][imol,parnum,regnum,:])+self.feature_dict['baseline'][imol,parnum,regnum,:]
        start_platpt=exceeds_thresh_points(traces, self.ssend_list, midpoints,operator.gt) #earliest point that trace exceeds threshold
        end_platpt=exceeds_thresh_points(traces,peakpt,midpoints,operator.lt,filt_length=filt_length) #earliest point that trace exceeds threshold AFTER the peak.  shold this be .lt
        self.feature_dict['start_plateau'][imol,parnum,regnum,:]=[platpt*self.dt[mol]/self.feature_scale['start_plateau'] for platpt in start_platpt]
        print_end = self.feature_dict['end_plateau'][imol,parnum,regnum,:]=[platpt*self.dt[mol]/self.feature_scale['start_plateau'] for platpt in end_platpt]
        #print('DURATION: mol',mol,'param',parnum,'start pt',start_platpt,'end',[round(p,2) for p in print_end])
        self.feature_dict['duration'][imol,parnum,regnum,:]=[(end-start)*self.dt[mol]/self.feature_scale['start_plateau']
                                                    for end,start in zip(end_platpt,start_platpt)]
        ####################
    def auc(self, traces, imol, parnum, par, mol, numstim, std_factor, aucend, iti, end_baseline_start, filt_length, regnum):                            
        import operator
        # CALCULATE AUC, using baseline+stdev as threshold - possibly use this for plateau?
        #also use the latest stimulation time if specified
        stim_time=self.ssend[mol] #default.  Overwrite under some conditions
        if numstim>1:
            if iti:
                stim_time=int((float(iti)*(numstim-1))/self.dt[mol]+self.ssend[mol])
            elif len(par):
                testset=[i for i in par[-1] ]
                if np.all([i in '0123456789.' for i in testset]):
                    stim_time=int((float(par[-1])*(numstim-1))/self.dt[mol]+self.ssend[mol])#get previous to last stimuation time
        if end_baseline_start:
            basestart=int(end_baseline_start/self.dt[mol])
            end_baseline=np.mean(traces[:,basestart:],axis=1)
            end_basestd=np.std(traces[:,basestart:],axis=1)
            self.feature_dict['auc_thresh'][imol,parnum,regnum,:]=end_baseline+std_factor*end_basestd
            baseline=end_baseline
        else:
            self.feature_dict['auc_thresh'][imol,parnum,regnum,:]=(self.feature_dict['baseline'][imol,parnum,regnum,:]+
                                                            std_factor*self.feature_dict['basestd'][imol,parnum,regnum,:])####
            baseline=self.feature_dict['baseline'][imol,parnum,regnum]
        #find latest point prior to molecule increase; find earliest point (after lastest stim) that molecule conc dips below basal
        peakpt_stim=np.argmax(traces[:,stim_time:],axis=1)+stim_time
        #could start auc from 1st point not belowthresh, but currently not using begin_auc
        #need to pass np.max into function to use begin_auc
        begin_auc=exceeds_thresh_points(traces, self.ssend_list,
                                            self.feature_dict['auc_thresh'][imol,parnum,regnum,:],
                                            operator.gt,peakpt_stim) #earliest time that trace is above threshold, 
        end_auc=exceeds_thresh_points(traces,peakpt_stim,
                                            self.feature_dict['auc_thresh'][imol,parnum,regnum,:],
                                            operator.lt,filt_length=filt_length) #earliest time that trace drops below threshold
        #print('AUC: mol=',mol,'param=',par,'start pt=',begin_auc,'peakpt_stim=',peakpt_stim,'end=',end_auc)
        if aucend is not None:
            self.feature_dict['auc'][imol,parnum,regnum,:]=[np.sum(traces[i,self.ssend[mol]:int(aucend/self.dt[mol])]-b)*self.dt[mol]/self.feature_scale['auc'] for i,b in enumerate(baseline)]
            print('end_auc',end_auc,'specified end auc',aucend)
        elif np.any(np.isnan(end_auc)):
            self.feature_dict['auc'][imol,parnum,regnum,:]=[np.sum(traces[i,self.ssend[mol]:]-b)*self.dt[mol] for i,b in enumerate(baseline)]
            print ('*********',mol,' is not returning to basal, calculating AUC to end of simulation, possibly raise your threshold **********')
        else:
            self.feature_dict['auc'][imol,parnum,regnum,:]=[np.sum(traces[i,self.ssend[mol]:end]-
                                                            baseline[i])*self.dt[mol] for i,end in enumerate(end_auc)]
            #FIXME: TEST using begin_auc
            #self.feature_dict['auc'][imol,parnum,regnum,:]=[np.sum(traces[i,begin:end]-
            #                                                       baseline[i])*self.dt[mol] for i,(begin,end) in enumerate(zip([begin_auc,end_auc]))]
 
    def write_features(self,feature_list,fileroot,regions,write_trials=False): #arg0=fileroot
        outfnam=self.savedir+os.path.basename(fileroot)+'-'+'analysis'+'-'.join([i for i in self.params])
        for regnum,reg in enumerate(self.all_regions): #average, one file per region, all molecules in each file
            if reg in regions: #only output a subset of regions
                outfname=outfnam+'-'+reg+'.txt'  #test whether this works for single file
                print('in write_features',outfname, feature_list,reg,regnum)
                if len(self.ftuples)==1:
                    outputdata=fileroot
                    header='file      ' 
                else:
                    outputdata=['-'.join(par) for par in self.par_keys]
                    header='-'.join([i for i in self.params])+'  '
                for feat in feature_list:
                    outputdata=np.column_stack((outputdata,np.round(self.mean_feature[feat][:,:,regnum].T/ms_to_sec,5),np.round(self.std_feature[feat][:,:,regnum].T/ms_to_sec,5)))       
                    header+=' '.join([m+'_'+feat+'_mean ' for m in list(self.molecules)+self.tot_species ]+ \
                                    [m+'_'+feat+'_std ' for m in list(self.molecules)+self.tot_species] )
                np.savetxt(outfname,outputdata,fmt='%1s', delimiter='     ', header=header)
        #write individual trials, one file per molecule, all features on one line, each parameter and trial on separate line
        if write_trials:
            #print('writing trials')
            par_string=['-'.join(par)+'  '+trial for par in self.par_keys for trial in self.trials]
            for imol,mol in enumerate(list(self.molecules)+self.tot_species):
                outfname=outfnam+'-'+mol+'-trials.txt'
                header='param  trial '
                outputdata=par_string
                for feat in feature_list:
                    for regnum,reg in enumerate(self.all_regions):
                        if reg in regions:
                            outputdata=np.column_stack((outputdata,np.round(self.feature_dict[feat][imol,:,regnum,:]/ms_to_sec,5).flatten()))
                            header=header+'   '+'_'.join([mol,reg,feat])
                np.savetxt(outfname,outputdata,fmt='%1s', delimiter='     ', header=header)

    def norm(self,sig_molecules,regnum,region,num_denom,min_max=None):
        all_dts=[self.dt[mol] for mol in sig_molecules] 
        if len(np.unique(all_dts))>1:
            #select largest dt
            dt_index=np.argmax(all_dts)
        else:
            dt_index=0
        self.dt[num_denom]=all_dts[dt_index]
        for parnum,par in enumerate(self.par_keys):
            if sig_molecules[dt_index] in self.file_set_tot[par][region].keys():                
                trace_length=np.shape(self.file_set_tot[par][region][sig_molecules[dt_index]])[-1]
            elif sig_molecules[dt_index] in self.file_set_conc[par][region].keys():               
                trace_length=np.shape(self.file_set_conc[par][region][sig_molecules[dt_index]])[-1]
            else:
                print('nrd_group, line 240. sig molecule',sig_molecules[dt_index],'not found in file_set_tot or file_set_conc')
            self.norm_traces[par][region][num_denom]=np.zeros((len(sig_molecules),len(self.trials),trace_length))
            for jmol,mol in enumerate(sig_molecules):
                if mol in self.molecules:
                    imol=self.molecules.index(mol)
                    traces=self.file_set_conc 
                elif mol in self.tot_species:
                    imol=self.tot_species.index(mol)+len(self.molecules)
                    traces=self.file_set_tot 
                #print('NORM, line 244, par=',par,'mol=',mol,'imol=',imol,jmol,'trace len',trace_length,'dt',self.dt[mol])
                #maxVal=np.max(self.feature_dict['peakval'][imol,:,regnum,:]) #max across trials - below first takes the mean, then takes max across protocols
                #minVal=np.min(self.feature_dict['minval'][imol,:regnum,:]) #same as above
                if min_max:
                    minVal=min_max[mol][region]['min']
                    maxVal=min_max[mol][region]['max']
                else:
                    maxVal=np.max(np.mean(self.feature_dict['peakval'][imol,:,regnum,:],axis=-1)) #replace first : with parnum to normalize separately for each protocol
                    minVal=np.min(np.mean(self.feature_dict['minval'][imol,:regnum,:],axis=-1)) #same as above
                if parnum==0:
                    temp_dict={mol:{region:{'max':round(maxVal,3),'min':round(minVal,3)}}}
                    print('norm constants for', num_denom, temp_dict) #print in format needed for tot_species, minmax
                    #print('norm constants for', num_denom, 'mol=', mol,'region=',region, 'max=', maxVal,'min=', minVal)
                for t in range(len(self.trials)):
                    # constrain norm between -1 and 1: 
                    if self.dt[mol]==self.dt[num_denom]:
                        new_trace=traces[par][region][mol][t,:]
                    else:
                        #interpolate traces
                        target_t=np.arange(trace_length)*self.dt[num_denom]
                        trace_t=np.arange(len(traces[par][region][mol][t,:]))*self.dt[mol]
                        new_trace=np.interp(target_t,trace_t,traces[par][region][mol][t,:])
                    self.norm_traces[par][region][num_denom][jmol,t]=(new_trace-self.feature_dict['baseline'][imol,parnum,regnum,t])/(maxVal-minVal)
                self.sstart[num_denom]=list(self.sstart.values())[dt_index]
                self.ssend[num_denom]=list(self.ssend.values())[dt_index]  
    
    def signature(self,mol,region,regnum,thresh,window=5):
        import operator
        for parnum,par in enumerate(self.par_keys):
            #both numerator and denom will be between -1 and 1, with baseline = 0
            if mol.startswith('product'):
                numerator=np.prod(self.norm_traces[par][region]['numerator'],axis=0) #product of molecules, dimensions are trials x time
            else:
                numerator=np.mean(self.norm_traces[par][region]['numerator'],axis=0) #average over molecules, dimensions are trials x time
            if len(self.norm_traces[par][region]['denom']):
                denom=np.mean(self.norm_traces[par][region]['denom'],axis=0) #average over molecules, dimensions are trials x time
            else:
                denom=1
            self.sig[mol][par][region]=numerator/denom  #2D array - trial x time
            self.dt[mol]=self.dt['numerator']
            self.sstart[mol]=self.sstart['numerator']
            self.ssend[mol]=self.ssend['numerator']
            #
            #now calculate features
            #
            self.sig_features['basestd'][mol][par][region]=np.std(self.sig[mol][par][region][:,self.sstart[mol]:self.ssend[mol]],axis=1)
            peakpt=np.argmax(self.sig[mol][par][region][:,self.ssend[mol]:],axis=1)+self.ssend[mol]
            self.sig_features['peaktime'][mol][par][region]=peakpt*self.dt[mol]/self.sig_feat_scale['peaktime']
            #this should be 1.0 unless there is a bug
            self.sig_features['amplitude'][mol][par][region]=[np.mean(self.sig[mol][par][region][i,p-window:p+window]) for i,p in enumerate(peakpt)] 
            #thresh_val=[thresh[region]*amp for amp in self.sig_features['amplitude'][mol][par][region]]
            #use this if specify different thresholds for each key in signature
            thresh_val=[thresh[mol][region] for t in range(len(self.trials))] 
            start_platpt=exceeds_thresh_points(self.sig[mol][par][region], self.ssend_list, thresh_val,operator.gt) #earliest point that trace exceeds threshold
            end_platpt=exceeds_thresh_points(self.sig[mol][par][region],peakpt,thresh_val,operator.lt) #earliest point that trace drops below threshold AFTER the peak.
            for i,ep in enumerate(end_platpt):
                if np.isnan(ep):
                    end_platpt[i]=len(self.sig[mol][par][region][i])
            self.sig_features['duration'][mol][par][region]=[(end-start)*self.dt[mol]/self.sig_feat_scale['duration']
                                                    for end,start in zip(end_platpt,start_platpt)]
            self.sig_features['auc'][mol][par][region]=[np.sum(self.sig[mol][par][region][i,self.ssend[mol]:end]-thresh_val[i])*self.dt[mol]/self.sig_feat_scale['auc']
                                                    for i,end in enumerate(end_platpt)] #sum area above the threshold
            self.sig_features['start_plateau'][mol][par][region]= [p*self.dt[mol]/self.sig_feat_scale['start_plateau'] for p in start_platpt]  
            self.sig_features['end_plateau'][mol][par][region]= [p*self.dt[mol]/self.sig_feat_scale['end_plateau'] for p in end_platpt]  
            start_dip_pt=exceeds_thresh_points(self.sig[mol][par][region], end_platpt, [0]* len(end_platpt), operator.lt) #earliest point that trace reaches below baseline
            self.sig_features['start_dip'][mol][par][region]= [p*self.dt[mol]/self.sig_feat_scale['start_dip'] for p in start_dip_pt]  

    def norm_sig(self,signature,thresh,min_max):
        self.sig={key:{p:{} for p in self.par_keys} for key in signature.keys()}
        self.sig_feature_list=['basestd','amplitude','duration','auc','start_plateau','end_plateau','peaktime', 'start_dip']
        self.sig_feat_scale={'basestd':1,'amplitude':1,'duration':ms_to_sec,'auc':ms_to_sec,'start_plateau':ms_to_sec,'end_plateau':ms_to_sec,'peaktime':ms_to_sec, 'start_dip':ms_to_sec}
        self.sig_features={feat:{key:{p:{} for p in self.par_keys} for key in signature.keys()} for feat in self.sig_feature_list}
        for key,sig in signature.items():
            num_molecules=sig['num']
            denom_molecules=sig['denom']
            self.norm_traces={p:{reg:{'numerator':{},'denom':{}} for reg in thresh[key].keys()} for p in self.par_keys}
            for region in thresh[key].keys():
                regnum=self.all_regions.index(region)
                if len(min_max) and key in min_max:
                    self.norm(num_molecules,regnum,region,'numerator',min_max[key]['num'])
                else:
                    self.norm(num_molecules,regnum,region,'numerator')
                if len(denom_molecules):
                    if len(min_max) and key in min_max:
                        self.norm(denom_molecules,regnum,region,'denom', min_max[key]['denom'])
                    else:
                        self.norm(denom_molecules,regnum,region,'denom')
            for regnum,region in enumerate(thresh[key].keys()):
                self.signature(key,region,regnum,thresh)

    def write_sig(self, regions, params,feature_list):  #one file per signature and parameter, all regions, average across trials. 
        for key in self.sig.keys():
            for par,ftuple in zip(self.sig[key].keys(),self.ftuples):
                par_str='-'.join(par)
                #self.sig[mol][par][region] #2D array - trial x time
                outfilename=self.savedir+os.path.splitext(os.path.basename(ftuple[0]))[0]+'-'+key
                columns=['_'.join([key,par_str,reg,tp]) for reg in regions for tp in ['mean','std'] ] 
                header='Time   '+'    '.join(columns)
                output_sig=self.dt[key]*np.arange(np.shape(self.sig[key][par][regions[0]])[-1]) 
                for reg in regions:
                    output_sig=np.column_stack((output_sig,np.mean(self.sig[key][par][reg],axis=0),np.std(self.sig[key][par][reg],axis=0)))
                np.savetxt(outfilename+'-sig.txt', output_sig, fmt='%.4f', delimiter=' ', header=header) #write signature trials with different filename
                if params.write_trials:
                    trial_cols=['_'.join([key,par_str,reg,tr]) for reg in regions for tr in self.trials ]
                    trial_header='Time   '+'    '.join(trial_cols)
                    output_trials= self.dt[key]*np.arange(np.shape(self.sig[key][par][regions[0]])[-1]) 
                    for reg in regions:
                        output_trials=np.column_stack((output_trials,self.sig[key][par][reg].T)) #need to transpose self.sig??  
                np.savetxt(outfilename+'-sig-trials.txt', output_trials, fmt='%.4f', delimiter=' ', header=trial_header) #write signature trials with different filename           
        print('in write_sig - features',feature_list,regions)
        #write sig features for individual trials, one file per molecule, all features on one line, each parameter and trial on separate line
        if params.write_trials:
            #write signature features for individual trials
            out_par=['-'.join(par)+'  '+trial for par in self.par_keys for trial in self.trials]
            for mol in self.sig_features[feature_list[0]].keys(): 
                outfname=self.savedir+os.path.basename(params.fileroot)+'-'+'anal'+'-'.join([i for i in self.params])+'-'+mol+'-trials.txt'
                header='param  trial '
                trials=len(self.trials)
                nfeat=len(feature_list)
                for ifeat,feat in enumerate(feature_list): #each feature separate column
                    for ipar,par in enumerate(self.sig_features[feat][mol].keys()): #different files are different rows
                        #row stack the pars
                        nrows=len(self.sig_features[feat][mol])*len(self.trials) #number of files * number of trials
                        ncols=len(self.sig_features[feat][mol][par])*nfeat #number of regions
                        output_array=np.empty((nrows,ncols))
                        for ireg,reg in enumerate(self.sig_features[feat][mol][par].keys()): #each region a different column
                            #column stack the regions
                            output_array[ipar*trials:(ipar+1)*trials,ireg*nfeat:ireg*nfeat+1]=np.round(np.array(self.sig_features[feat][mol][par][reg]),2).reshape(trials,1)
                            if ipar==0:
                                 header=header+'   '+'_'.join([mol,reg,feat])
                outputdata=np.column_stack((out_par,output_array))
                np.savetxt(outfname,outputdata,fmt='%1s', delimiter='     ', header=header)


    #one file per parameter and molecule, only regions of interest.
    def write_trace_trials(self, regions,fileroot):
        import os 
        #### write trials for specified regions providing they are in file_set_conc
        for par in self.file_set_conc.keys():
            par_name='-'.join([str(q) for q in par])
            reg0=regions[0] 
            for file_set in [self.file_set_conc, self.file_set_tot]:
                for mol in file_set[par][reg0].keys(): 
                    outfilename=self.savedir+os.path.splitext(os.path.basename(fileroot))[0]+'-'+par_name+'-'+mol+'-trials.txt'
                    header='Time '
                    if mol in self.time_set[par].keys():
                        output=self.time_set[par][mol]
                    else:
                        output=np.arange(len(file_set[par][reg][mol].T))*self.dt[mol]
                    for reg in regions:  
                        output=np.column_stack((output,file_set[par][reg][mol].T))
                        header=header+' '.join(['_'.join([par_name,mol,reg,t]) for t in self.trials])+' '
                    header=header
                    np.savetxt(outfilename, output, fmt='%.4f', delimiter='     ',header=header)           

def exceeds_thresh_points(traces,startpoints,thresh,relate,endpoint=-1,filt_length=0):
    #Find earliest point when traces (from startpoint to endpoint) is over or under the threshold
    #relate is either > (operator.gt) or < (operator.lt)
    #need to replace np.min with function in case want to find latest point
    if not np.isscalar(endpoint):
        endpoint_list=endpoint
    else:
        endpoint_list=[endpoint for t in range(len(traces))]
    
    earliest_points=[np.nan for i in startpoints]
    
    #print('start',startpoints,'thresh',[round(t,2) for t in thresh],'traces',np.shape(traces),'endpoint', [round(ep,3) for ep in endpoint_list])

    for i,(sp,t,endpt) in enumerate(zip(startpoints,thresh,endpoint_list)):
        if not np.isnan(sp):
            if filt_length>0:
                mov_avg_filt=np.ones(filt_length)/filt_length
                newtraces=np.convolve(mov_avg_filt,traces[i],'same')
            else:
                newtraces=traces[i]
            pointset=np.where(relate(newtraces[sp:endpt],t))[0]+sp
            if len(pointset):
                earliest_points[i]=np.min(pointset)
    return earliest_points




