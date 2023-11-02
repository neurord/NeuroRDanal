# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:59:11 2020

@author: kblackw1
"""
import numpy as np
from NeuroRDanal import h5utilsV2 as h5utils

ms_to_sec=1000
class nrdh5_group(object): 
    def __init__(self,fileroot,parameters,tot_species=[]):
        self.ftuples,self.parlist,self.params=h5utils.create_filenames(fileroot,parameters)
        self.file_set_conc={'Overall':{k[1]:{} for k in self.ftuples}}
        self.time_set={k[1]:{} for k in self.ftuples}
        self.means={reg:{k[1]:{} for k in self.ftuples} for reg  in ['space']}

        if len(tot_species):
            self.tot_species=tot_species
            self.endtime={k[1]:{sp:[] for sp in self.tot_species} for k in self.ftuples}
            self.file_set_tot={'Overall':{k[1]:{sp:[] for sp in self.tot_species} for k in self.ftuples}}
        else:
            self.tot_species=[]
            self.file_set_tot={'Overall':{}}

    def conc_arrays(self,data):
        self.molecules=data.molecules
        self.trials=data.trials
        #These are overwritten with each data file, but must be the same for each data file
        self.sstart={mol:data.sstart[mol] for mol in data.molecules}
        self.ssend={mol:data.ssend[mol] for mol in data.molecules}
        self.dt={mol:data.dt[mol] for mol in data.molecules}
        if data.maxvols>1:
            for reg_dict in [data.region_dict,data.region_struct_dict]:
                for regnum,reg in enumerate(reg_dict.keys()):
                    if reg not in self.file_set_conc.keys():
                        self.file_set_conc[reg]={k[1]:{} for k in self.ftuples}
            if data.spinelist:
                for regnum,reg in enumerate(data.spinelist):
                    if reg not in self.file_set_conc.keys():
                        self.file_set_conc[reg]={k[1]:{} for k in self.ftuples}
        for imol,molecule in enumerate(data.molecules):
            self.time_set[data.parval][molecule]=data.time[molecule]
            self.file_set_conc['Overall'][data.parval][molecule]=data.OverallMean[molecule]
            if data.maxvols>1:
                for reg_type,reg_dict in zip(['region','struct'],[data.region_dict,data.region_struct_dict]):
                    for regnum,reg in enumerate(reg_dict.keys()):
                        self.file_set_conc[reg][data.parval][molecule]=data.means[reg_type][molecule][:,:,regnum]
                if data.spinelist:
                     for regnum,reg in enumerate(data.spinelist):
                        self.file_set_conc[reg][data.parval][molecule]=data.means['spines'][molecule][:,:,regnum]
                if data.spatial_dict:
                    self.means['space'][data.parval][molecule]=data.means['space'][molecule]
            else:
                self.spatial_data=None
        if len(self.tot_species):
            for region in data.total_trace.keys():
                if region not in self.file_set_tot.keys():
                    self.file_set_tot[region]={k[1]:{sp:[] for sp in self.tot_species} for k in self.ftuples}
                for imol,sp in enumerate(self.file_set_tot[region][data.parval].keys()):
                    self.file_set_tot[region][data.parval][sp]=data.total_trace[region][sp][:,:]
                    self.sstart[sp]=data.sstart[sp]
                    self.ssend[sp]=data.ssend[sp]
                    self.dt[sp]=data.dt[sp]
                    self.endtime[data.parval][sp]=data.endtime[sp]
        
    def trace_features(self,window_size,lo_thresh_factor=0.2,hi_thresh_factor=0.8,std_factor=2,numstim=1,end_baseline_start=0,filt_length=5,aucend=None,iti=None):
        import operator
        self.feature_list=['baseline','basestd','peakval','peaktime','amplitude','duration','slope','minval','auc','auc_thresh','start_plateau', 'end_plateau']
        self.feature_dict={feat:np.zeros((len(self.molecules)+len(self.tot_species),len(self.ftuples), len(self.file_set_conc.keys()),len(self.trials))) for feat in self.feature_list}
        self.mean_feature={}
        self.std_feature={}

        molecules=self.molecules #use individual molecules for self.file_set_conc and self.means
        ii = 0
        for regnum, (region,traces) in enumerate(self.file_set_conc.items()):
            #print('*************',region,regnum,traces.keys())
            self.tf_base_function(window_size, molecules, traces, ii, regnum)

        molecules=self.tot_species
        ii=1
        for regnum, (region,traces) in enumerate(self.file_set_tot.items()):
            #print('****** nrd_group line 84 *******',region,regnum,traces.keys())
            self.tf_base_function(window_size, molecules, traces, ii, regnum)

        for feat in self.feature_dict.keys():
            self.mean_feature[feat]=np.nanmean(self.feature_dict[feat],axis=-1)
            self.std_feature[feat]=np.nanstd(self.feature_dict[feat],axis=-1)

    def tf_base_function(self,window_size, molecules, traces, ii, regnum, lo_thresh_factor=0.2,hi_thresh_factor=0.8,std_factor=2,numstim=1,end_baseline_start=0,filt_length=5,aucend=None,iti=None):
        import operator
        for parnum,(fname,par) in enumerate(self.ftuples):
            for jmol,mol in enumerate(molecules):
                #print(par,mol,np.shape(traces[par][mol]))
                imol=jmol+ii*len(self.molecules)
                window=int(window_size/self.dt[mol]) #FIXME self.dt[mol],sstart[mol],ssend[mol]
                if window==0:
                    print('trace_features_loop, window size too small for',par,mol,'dt',self.dt[mol],'window size',window_size,'using',self.dt[mol])
                    window=1
                self.feature_dict['baseline'][imol,parnum,regnum,:]= np.mean(traces[par][mol][:,self.sstart[mol]:self.ssend[mol]],axis=1)
                self.feature_dict['basestd'][imol,parnum,regnum,:]=np.std(traces[par][mol][:,self.sstart[mol]:self.ssend[mol]],axis=1)
                peakpt=np.argmax(traces[par][mol][:,self.ssend[mol]:],axis=1)+self.ssend[mol]
                self.feature_dict['peaktime'][imol,parnum,regnum,:]=peakpt*self.dt[mol]
                self.feature_dict['peakval'][imol,parnum,regnum,:]=[np.mean(traces[par][mol][i,p-window:p+window]) 
                                                             for i,p in enumerate(peakpt)]
                lowpt=np.argmin(traces[par][mol][:,self.ssend[mol]:],axis=1)+self.ssend[mol]
                #DEBUGGING: peak time and peakval from mean of trace
                p=np.argmax(np.mean(traces[par][mol][:,self.ssend[mol]:],axis=0))+self.ssend[mol]
                pt=p*self.dt[mol]
                pval=np.mean(np.mean(traces[par][mol],axis=0)[p-window:p+window])
                #print('peaktime {} & peakval {} from mean trace'.format(pt,pval))
                #print('lowpt',mol,lowpt,'window',window,'lowval mean',np.mean([traces[par][mol][i,p-window:p+window] for i,p in enumerate(lowpt)]))
                #end DEBUGGING
                self.feature_dict['minval'][imol,parnum,regnum,:]=[np.mean(traces[par][mol][i,p-window:p+window]) 
                                                            for i,p in enumerate(lowpt)]
                self.feature_dict['amplitude'][imol,parnum,regnum,:]=self.feature_dict['peakval'][imol,parnum,regnum,:]-self.feature_dict['baseline'][imol,parnum,regnum,:]
                ####################

                self.slope(traces[par][mol], imol, parnum, mol, regnum)
                self.plateau_duration(traces[par][mol], imol, parnum, mol, peakpt, filt_length, regnum)
                self.auc(traces[par][mol], imol, parnum, par, mol, numstim, std_factor, aucend, iti, end_baseline_start, filt_length, regnum)
                
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
        self.feature_dict['start_plateau'][imol,parnum,regnum,:]=[platpt*self.dt[mol] for platpt in start_platpt]
        print_end = self.feature_dict['end_plateau'][imol,parnum,regnum,:]=[platpt*self.dt[mol] for platpt in end_platpt]
        print('DURATION: mol',mol,'param',parnum,'start pt',start_platpt,'end',[round(p,2) for p in print_end])
        self.feature_dict['duration'][imol,parnum,regnum,:]=[(end-start)*self.dt[mol]
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
        print('AUC: mol=',mol,'param=',par,'start pt=',begin_auc,'peakpt_stim=',peakpt_stim,'end=',end_auc)
        if aucend is not None:
            self.feature_dict['auc'][imol,parnum,regnum,:]=[np.sum(traces[i,self.ssend[mol]:int(aucend/self.dt[mol])]-b)*self.dt[mol] for i,b in enumerate(baseline)]
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
 
    def write_features(self,feature_list,arg0,write_trials=False): #FIXME - not writing correct values
        for regnum,reg in enumerate(self.file_set_tot.keys()):
            outfname=arg0+'-'+'analysis'+'-'.join([i for i in self.params])+'_'+reg+'.txt'  #test whether this works for single file
            print('in write_features',outfname, feature_list,reg,regnum)
            if len(self.ftuples)==1:
                outputdata=arg0
                header='file      ' 
            else:
                outputdata=['-'.join([str(p) for p in par[1]]) for par in self.ftuples]
                header='-'.join([i for i in self.params])+'  '
            header+=' '.join([m+'_'+f+'_mean ' for f in feature_list for m in list(self.molecules)+self.tot_species ]+ \
                               [m+'_'+f+'_std ' for f in feature_list for m in list(self.molecules)+self.tot_species] )
            for feat in feature_list:
                outputdata=np.column_stack((outputdata,np.round(self.mean_feature[feat][:,:,regnum].T/ms_to_sec,3),np.round(self.std_feature[feat][:,:,regnum].T/ms_to_sec,3)))
            f=open(outfname, 'w')
            f.write(header+'\n')        
            np.savetxt(f,outputdata,fmt='%1s', delimiter='     ')
            f.close() 

        #write individual trials
        if write_trials:
            #print('writing trials')
            for imol,mol in enumerate(self.molecules):
                outfname=arg0+'-'+'analysis'+'-'.join([i for i in self.params])+'-'+mol+'-trials.txt'
                params=[ftuple[1] for ftuple in self.ftuples]
                output_data=np.column_stack((params,np.round(self.feature_dict[feat][imol,:,:]/ms_to_sec,5)))
                header='param '+' '.join(['trial'+str(n) for n in range (np.shape(self.feature_dict[feat])[2])])
                f=open(outfname, 'w')
                f.write(header+'\n')        
                np.savetxt(f,output_data,fmt='%1s', delimiter='     ')
                f.close()

    def norm(self,sig_molecules,regnum,region,num_denom):
        all_dts=[self.dt[mol] for mol in sig_molecules] 
        if len(np.unique(all_dts))>1:
            #select largest dt
            dt_index=np.argmax(list(self.dt.values()))
        else:
            dt_index=0
        self.dt[num_denom]=all_dts[dt_index]
        for parnum,(fname,par) in enumerate(self.ftuples):
            if sig_molecules[dt_index] in self.file_set_tot[region][par].keys():                
                trace_length=np.shape(self.file_set_tot[region][par][sig_molecules[dt_index]])[-1]
            elif sig_molecules[dt_index] in self.file_set_conc[region][par].keys():               
                trace_length=np.shape(self.file_set_conc[region][par][sig_molecules[dt_index]])[-1]
            else:
                print('nrd_group, line 240. sig molecule',sig_molecules[dt_index],'not found in file_set_tot or file_set_conc')
            self.norm_traces[par][region][num_denom]=np.zeros((len(sig_molecules),len(self.trials),trace_length))
            for jmol,mol in enumerate(sig_molecules):
                if mol in self.molecules:
                    imol=self.molecules.index(mol)
                    traces=self.file_set_conc[region]
                elif mol in self.tot_species:
                    imol=self.tot_species.index(mol)+len(self.molecules)
                    traces=self.file_set_tot[region]
                print('NORM, line 244, par=',par,'mol=',mol,'imol=',imol,jmol,'trace len',trace_length,'dt',self.dt[mol])
                for t in range(len(self.trials)):
                    # constrain norm between -1 and 1: 
                    if self.dt[mol]==self.dt[num_denom]:
                        new_trace=traces[par][mol][t,:]
                    else:
                        #interpolate traces
                        target_t=np.arange(trace_length)*self.dt[num_denom]
                        trace_t=np.arange(len(traces[par][mol][t,:]))*self.dt[mol]
                        new_trace=np.interp(target_t,trace_t,traces[par][mol][t,:])
                    maxVal=np.max(self.feature_dict['peakval'][imol,parnum,regnum,:]) #replace parnum with : if normalized accross protcol
                    minVal=np.min(self.feature_dict['minval'][imol,parnum,regnum,:]) #same as above 
                    self.norm_traces[par][region][num_denom][jmol,t]=(new_trace-self.feature_dict['baseline'][imol,parnum,regnum,t])/(maxVal-minVal)
                self.sstart[num_denom]=list(self.sstart.values())[dt_index]
                self.ssend[num_denom]=list(self.ssend.values())[dt_index]  
    
    def signature(self,mol,region,regnum,thresh,window=5):
        import operator
        for parnum,(fname,par) in enumerate(self.ftuples):
            #both numerator and denom will be between -1 and 1, with baseline = 0
            if mol == 'product':
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
            self.sig_features['peaktime'][mol][par][region]=peakpt*self.dt[mol]
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
            self.sig_features['duration'][mol][par][region]=[(end-start)*self.dt[mol]
                                                    for end,start in zip(end_platpt,start_platpt)]
            self.sig_features['auc'][mol][par][region]=[np.sum(self.sig[mol][par][region][i,self.ssend[mol]:end]-thresh_val[i])*self.dt[mol] 
                                                    for i,end in enumerate(end_platpt)] #sum area above the threshold
            self.sig_features['start_plateau'][mol][par][region]= [p*self.dt[mol] for p in start_platpt]  
            self.sig_features['end_plateau'][mol][par][region]= [p*self.dt[mol] for p in end_platpt]  

    def norm_sig(self,signature,thresh):
        self.sig={key:{p[1]:{} for p in self.ftuples} for key in signature.keys()}
        self.sig_feature_list=['basestd','amplitude','duration','auc','start_plateau','end_plateau','peaktime']
        self.sig_features={feat:{key:{p[1]:{} for p in self.ftuples} for key in signature.keys()} for feat in self.sig_feature_list}
        for key,sig in signature.items():
            num_molecules=sig['num']
            denom_molecules=sig['denom']
            self.norm_traces={p[1]:{reg:{'numerator':{},'denom':{}} for reg in thresh[key].keys()} for p in self.ftuples}
            for region in thresh[key].keys():
                regnum=list(self.file_set_conc.keys()).index(region)
                self.norm(num_molecules,regnum,region,'numerator')
                if len(denom_molecules):
                    self.norm(denom_molecules,regnum,region,'denom')
            for regnum,region in enumerate(thresh[key].keys()):
                self.signature(key,region,regnum,thresh)

def exceeds_thresh_points(traces,startpoints,thresh,relate,endpoint=-1,filt_length=0):
    #Find earliest point when traces (from startpoint to endpoint) is over or under the threshold
    #relate is either > (operator.gt) or < (operator.lt)
    #need to replace np.min with function in case want to find latest point
    if not np.isscalar(endpoint):
        endpoint_list=endpoint
    else:
        endpoint_list=[endpoint for t in range(len(traces))]
    
    earliest_points=[np.nan for i in startpoints]
    
    print('start',startpoints,'thresh',[round(t,2) for t in thresh],'traces',np.shape(traces),'endpoint', [round(ep,3) for ep in endpoint_list])

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




