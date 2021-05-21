# -*- coding: utf-8 -*-
"""
Created on Tue Apr 28 15:59:11 2020

@author: kblackw1
"""
import numpy as np
from NeuroRDanal import h5utilsV2 as h5utils

ms_to_sec=1000
class nrdh5_group(object): 
    def __init__(self,args,tot_species=[]):
        self.ftuples,self.parlist,self.params=h5utils.argparse(args)
        self.file_set_conc={k[1]:{} for k in self.ftuples} 
        self.time_set={k[1]:{} for k in self.ftuples}
        self.spatial_means={k[1]:{} for k in self.ftuples}
        self.regions_means={k[1]:{} for k in self.ftuples}
        self.regions_structure_means={k[1]:{} for k in self.ftuples}
        self.spine_means={k[1]:{} for k in self.ftuples}
        if len(tot_species):
            #Expand this to include spine_means and other spatial information
            self.file_set_tot={k[1]:{sp:[] for sp in tot_species} for k in self.ftuples}
            self.tot_species=tot_species
            self.endtime={k[1]:{sp:[] for sp in tot_species} for k in self.ftuples}
        else:
            self.tot_species=[]
            self.file_set_tot={}

    def conc_arrays(self,data):
        self.molecules=data.molecules
        self.trials=data.trials
        #These are overwritten with each data file, but must be the same for each data file
        self.sstart={mol:data.sstart[mol] for mol in data.molecules}
        self.ssend={mol:data.ssend[mol] for mol in data.molecules}
        self.dt={mol:data.dt[mol] for mol in data.molecules}
        for imol,molecule in enumerate(data.molecules):
            self.time_set[data.parval][molecule]=data.time[molecule]
            self.file_set_conc[data.parval][molecule]=data.OverallMean[molecule]
            if data.maxvols>1:
                self.regions_means[data.parval][molecule]=data.means['struct'][molecule]
                self.regions_structure_means[data.parval][molecule]=data.means['region'][molecule]
                if data.spatial_dict:
                    self.spatial_means[data.parval][molecule]=data.means['space'][molecule]
                if data.spinelist:
                    self.spine_means[data.parval][molecule]=data.means['spines'][molecule]
            else:
                self.spatial_data=None
        if len(data.ss_tot):
            #Expand this to include spine_means and other spatial information
            for imol,sp in enumerate(self.file_set_tot[data.parval].keys()):
                self.file_set_tot[data.parval][sp]=data.ss_tot[sp][:,:]
                self.sstart[sp]=data.sstart[sp]
                self.ssend[sp]=data.ssend[sp]
                self.dt[sp]=data.dt[sp]
                self.endtime[data.parval][sp]=data.endtime[sp]
        
    def trace_features(self,trials,window_size,lo_thresh_factor=0.2,hi_thresh_factor=0.8,std_factor=1,numstim=1,end_baseline_start=0,filt_length=5,aucend=None):
        import operator
        self.feature_list=['baseline','basestd','peakval','peaktime','amplitude','duration','slope','minval','auc','auc_thresh','start_plateau']
        self.feature_dict={feat:np.zeros((len(self.molecules)+len(self.tot_species),len(self.ftuples),len(trials))) for feat in self.feature_list}
        self.mean_feature={feat:np.zeros((len(self.molecules)+len(self.tot_species),len(self.ftuples))) for feat in self.feature_list}
        self.std_feature={feat:np.zeros((len(self.molecules)+len(self.tot_species),len(self.ftuples))) for feat in self.feature_list}

        def exceeds_thresh_points(traces,startpoints,thresh,relate,endpoints=[-1 for t in trials],filt_length=0):
            #Find earliest point when traces (from startpoint to endpoint) is over or under the threshold
            #relate is either > (operator.gt) or < (operator.lt)
            #need to replace np.min with function in case want to find latest point
            #print('start',startpoints,'thresh',thresh,'traces',np.shape(traces))
            earliest_points=[np.nan for i in startpoints]
            for i,(sp,t,endpt) in enumerate(zip(startpoints,thresh,endpoints)):
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
        for ii,(molecules,traces) in enumerate(zip([self.molecules,self.tot_species],[self.file_set_conc,self.file_set_tot])):
            for parnum,(fname,par) in enumerate(self.ftuples):
                for jmol,mol in enumerate(molecules):
                    imol=jmol+ii*len(self.molecules)
                    window=int(window_size/self.dt[mol]) #FIXME self.dt[mol],sstart[mol],ssend[mol]
                    if window==0:
                        print('trace_features_loop, window size too small for',par,mol,'dt',self.dt[mol],'window size',window_size,'using',self.dt[mol])
                        window=1
                    self.feature_dict['baseline'][imol,parnum,:]= np.mean(traces[par][mol][:,self.sstart[mol]:self.ssend[mol]],axis=1)
                    self.feature_dict['basestd'][imol,parnum,:]=np.std(traces[par][mol][:,self.sstart[mol]:self.ssend[mol]],axis=1)
                    peakpt=np.argmax(traces[par][mol][:,self.ssend[mol]:],axis=1)+self.ssend[mol]
                    self.feature_dict['peaktime'][imol,parnum,:]=peakpt*self.dt[mol]
                    self.feature_dict['peakval'][imol,parnum,:]=[np.mean(traces[par][mol][i,p-window:p+window]) 
                                                                 for i,p in enumerate(peakpt)]
                    lowpt=np.argmin(traces[par][mol][:,self.ssend[mol]:],axis=1)+self.ssend[mol]
                    #DEBUGGING: peak time and peakval from mean of trace
                    p=np.argmax(np.mean(traces[par][mol][:,self.ssend[mol]:],axis=0))+self.ssend[mol]
                    pt=p*self.dt[mol]
                    pval=np.mean(np.mean(traces[par][mol],axis=0)[p-window:p+window])
                    #print('peaktime {} & peakval {} from mean trace'.format(pt,pval))
                    #print('lowpt',mol,lowpt,'window',window,'lowval mean',np.mean([traces[par][mol][i,p-window:p+window] for i,p in enumerate(lowpt)]))
                    #end DEBUGGING
                    self.feature_dict['minval'][imol,parnum,:]=[np.mean(traces[par][mol][i,p-window:p+window]) 
                                                                for i,p in enumerate(lowpt)]
                    self.feature_dict['amplitude'][imol,parnum,:]=self.feature_dict['peakval'][imol,parnum,:]-self.feature_dict['baseline'][imol,parnum,:]
                    ####################
                    #FIND SLOPE OF INCREASE - Use thresholds defined by lo_thresh, and hi_thresh, e.g. 20 and 80%
                    lo_thresh=lo_thresh_factor*(self.feature_dict['amplitude'][imol,parnum,:])+self.feature_dict['baseline'][imol,parnum,:] #get the 5% above the max value
                    hi_thresh=hi_thresh_factor*(self.feature_dict['amplitude'][imol,parnum,:])+self.feature_dict['baseline'][imol,parnum,:]
                    ssend_list=[self.ssend[mol] for t in trials]
                    #
                    begin_slope=exceeds_thresh_points(traces[par][mol], ssend_list, lo_thresh,operator.gt)
                    end_slope=exceeds_thresh_points(traces[par][mol], begin_slope,hi_thresh,operator.gt)
                    for i,(end,beg) in enumerate(zip(end_slope,begin_slope)):
                        if end-beg>1 and ~np.isnan(end) and ~np.isnan(beg): #FIX - check for end=beg - making zero slope, or even end_slope.beginslope=1
                            self.feature_dict['slope'][imol,parnum,i]=(self.feature_dict['peakval'][imol,parnum,i]-
                                                                       self.feature_dict['baseline'][imol,parnum,i])/((end-beg)*self.dt[mol])
                        else:
                            self.feature_dict['slope'][imol,parnum,:]=np.nan                    
                    ####################
                    # FIND PLATEAU DURATION - USE thresholds defined by midpoints, and two different time periods
                    #could also use thresholds defined by lo_thresh or hi_thresh
                    midpoints=0.5*(self.feature_dict['amplitude'][imol,parnum,:])+self.feature_dict['baseline'][imol,parnum,:]
                    start_platpt=exceeds_thresh_points(traces[par][mol],ssend_list, midpoints,operator.gt)
                    end_platpt=exceeds_thresh_points(traces[par][mol],peakpt,midpoints,operator.gt,filt_length=filt_length)
                    self.feature_dict['start_plateau'][imol,parnum,:]=[platpt*self.dt[mol] for platpt in start_platpt]
                    #print('plateau',self.feature_dict['start_plateau'][imol,parnum,:],end_platpt)
                    self.feature_dict['duration'][imol,parnum,:]=[(end-start)*self.dt[mol] 
                                                                  for end,start in zip(end_platpt,start_platpt)]
                    ####################
                    # CALCULATE AUC, using baseline+stdev as threshold - possibly use this for plateau?
                    #also use the latest stimulation time if specified
                    if numstim>1 and len(par):
                        stim_time=int((float(par[-1])*(numstim-1))/self.dt[mol]+self.ssend[mol])#get previous to last stimuation time
                    else:
                        stim_time=self.ssend[mol]
                    if end_baseline_start:
                        basestart=int(end_baseline_start/self.dt[mol])
                        end_baseline=np.mean(traces[par][mol][:,basestart:],axis=1)
                        end_basestd=np.std(traces[par][mol][:,basestart:],axis=1)
                        self.feature_dict['auc_thresh'][imol,parnum,:]=end_baseline+std_factor*end_basestd
                        baseline=end_baseline
                    else:
                        self.feature_dict['auc_thresh'][imol,parnum,:]=(self.feature_dict['baseline'][imol,parnum,:]+
                                                                        std_factor*self.feature_dict['basestd'][imol,parnum,:])####
                        baseline=self.feature_dict['baseline'][imol,parnum]
                    #find latest point prior to molecule increase; find earliest point (after lastest stim) that molecule conc dips below basal
                    peakpt_stim=np.argmax(traces[par][mol][:,stim_time:],axis=1)+stim_time
                    #could start auc from 1st point not belowthresh, but currently not using begin_auc
                    #need to pass np.max into function to use begin_auc
                    begin_auc=exceeds_thresh_points(traces[par][mol],ssend_list,
                                                          self.feature_dict['auc_thresh'][imol,parnum,:],
                                                          operator.lt,peakpt_stim)
                    end_auc=exceeds_thresh_points(traces[par][mol],peakpt_stim,
                                                          self.feature_dict['auc_thresh'][imol,parnum,:],
                                                          operator.lt,filt_length=filt_length)
                    if aucend is not None:
                         self.feature_dict['auc'][imol,parnum,:]=[np.sum(traces[par][mol][i,self.ssend[mol]:int(aucend/self.dt[mol])]-b)*self.dt[mol] for b in baseline]
                         print('end_auc',end_auc,'specified end auc',aucend)
                    elif np.any(np.isnan(end_auc)):
                         self.feature_dict['auc'][imol,parnum,:]=[np.sum(traces[par][mol][i,self.ssend[mol]:]-b)*self.dt[mol] for b in baseline]
                         print ('*********',mol,' is not returning to basal, calculating AUC to end of simulation, possibly raise your threshold **********')
                    else:
                        self.feature_dict['auc'][imol,parnum,:]=[np.sum(traces[par][mol][i,self.ssend[mol]:end]-
                                                                        baseline[i])*self.dt[mol] for i,end in enumerate(end_auc)]

                    #calculate mean and stdev of various features
                    for feat in self.feature_dict.keys():
                        self.mean_feature[feat][imol,parnum]=np.nanmean(self.feature_dict[feat][imol,parnum,:])
                        self.std_feature[feat][imol,parnum]=np.nanstd(self.feature_dict[feat][imol,parnum,:])

    def write_features(self,feature_list,arg0,write_trials=False):
        outfname=arg0+'-'+'analysis'+'-'.join([i for i in self.params])+'.txt'  #test whether this works for single file
        #print('in write_features',outfname, feature_list)
        if len(self.ftuples)==1:
            outputdata=arg0
            header='file      ' 
        else:
            outputdata=['-'.join([str(p) for p in par[1]]) for par in self.ftuples]
            header='-'.join([i for i in self.params])+'  '
        header+=' '.join([m+'_'+f+'_mean ' for m in self.molecules+self.tot_species for f in feature_list]+ [m+'_'+f+'_std ' for m in self.molecules+self.tot_species for f in feature_list])
        for feat in feature_list:
            outputdata=np.column_stack((outputdata,np.round(self.mean_feature[feat].T/ms_to_sec,3),np.round(self.std_feature[feat].T/ms_to_sec,3)))
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



