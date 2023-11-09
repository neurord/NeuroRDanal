from __future__ import print_function
from __future__ import division
import numpy as np
from matplotlib import pyplot

legtextsize=9
     
colors=pyplot.get_cmap('viridis')
#colors=pyplot.get_cmap('plasma')
colors2D=[pyplot.get_cmap('gist_heat'),pyplot.get_cmap('summer'),pyplot.get_cmap('Blues')]
offset=[0,0,63]  #avoid the light colors in low indices for the 'Blues' map
partial_scale=0.75 #avoid the very light colors.  Note that non-zero offset must be <= (1-partial_scale)*255
ms_to_sec=1000

def xval_from_params(dataset):
    params=[f[1] for f in dataset.ftuples]
    if len(dataset.params)>1:
        if len(dataset.parlist[0])*len(dataset.parlist[1])==len(dataset.ftuples):
            if len(dataset.parlist[0])>len(dataset.parlist[1]):
                xval_index=0
            else:
                xval_index=1
            label_index=(xval_index+1)%2
            xvals=dataset.parlist[xval_index]
            labels={'labels':dataset.parlist[label_index],'xindex':xval_index,'lindex':label_index}
        else:
            xvals=[str(p[0])+'-'+str(p[1]) for p in params]
            labels=[]
    else:
        xvals=[str(p[0]) for p in params]
        labels=[]
    return xvals,labels

def plot_features(dataset,feature,title):
    xvals,p=xval_from_params(dataset)
    if len(dataset.tot_species):
        feature_mol=dataset.molecules+dataset.tot_species
    else:
        feature_mol=dataset.molecules
    for imol,mol in enumerate(feature_mol):
        pyplot.figure() #new figure panel for each molecules
        pyplot.suptitle(title)
        #FIXME need to loop over region, or something that that if len(p)
        if len(p): #reshape feature values for plotting if 2 params and all combos
            labels=p['labels']; xval_index=p['xindex'];label_index=p['lindex']
            new_yvals=np.zeros((len(xvals),len(labels)))
            for fnum,(fname,param) in enumerate(dataset.ftuples):
                row=xvals.index(param[xval_index])
                col=labels.index(param[label_index])
                new_yvals[row,col]=dataset.mean_feature[feature][imol,fnum]
            for col,label in enumerate(labels):
                pyplot.scatter(xvals,new_yvals[:,col],label=dataset.params[label_index]+' '+str(label))
            pyplot.xlabel(dataset.params[xval_index])
        else: #just plot feature values vs param or param combo            
            for regnum,key in enumerate(dataset.file_set_tot.keys()):           
                pyplot.scatter(xvals,dataset.mean_feature[feature][imol,:,regnum], label=key)
            pyplot.xlabel('-'.join(dataset.params))
        pyplot.ylabel(mol+' '+feature)
        pyplot.legend()
#IndexError: too many indices for array: array is 1-dimensional, but 2 were indexed

def spatial_plot(data,dataset,plot_trials=0):
    numtrials=len(data.trials)
    if plot_trials==0:
        fig,axes=pyplot.subplots(len(data.molecules),len(dataset.ftuples), sharex=True)
        
    for par,(fname,param) in enumerate(dataset.ftuples):
        if plot_trials==1:
            fig,axes=pyplot.subplots(len(data.molecules),numtrials+1, sharex=True)
            fig.suptitle(fname)
        axes=fig.axes
        for imol,mol in enumerate(data.molecules):
            if plot_trials==1:
                for trial in range(numtrials):
                   ax=imol*(numtrials+1)+trial
                   axes[ax].imshow(dataset.spatial_means[param][mol][trial].T, aspect='auto',origin='lower',
                                   extent=[0, np.max(dataset.time_set[param][mol]), float(list(data.spatial_dict.keys())[0]), float(list(data.spatial_dict.keys())[-1])])
                   #axes[imol,trial].colorbar()
                axes[imol*(numtrials+1)].set_ylabel (mol +', location (um)')
                #to plot the mean across trials
                ax=imol*(numtrials+1)+numtrials
                axes[ax].imshow(np.mean(dataset.spatial_means[param][mol][:],axis=0).T, aspect='auto',origin='lower',
                                   extent=[0, np.max(dataset.time_set[param][mol]), float(list(data.spatial_dict.keys())[0]), float(list(data.spatial_dict.keys())[-1])])
                for trial in range(numtrials+1):
                    axes[imol*(numtrials+1)+trial].set_xlabel('time (ms)')
            else:
                ax=imol*(len(dataset.ftuples))+par
                axes[ax].imshow(np.mean(dataset.spatial_means[param][mol][:],axis=0).T, aspect='auto',origin='lower',
                            extent=[0, np.max(dataset.time_set[param][mol]), float(list(data.spatial_dict.keys())[0]), float(list(data.spatial_dict.keys())[-1])])
            axes[imol*(len(dataset.ftuples))].set_ylabel (mol +', location (um)')
        axes[imol*(len(dataset.ftuples))+par].set_xlabel('time (ms)')
        axes[par].set_title(param[0])
        
def plot_setup(plot_molecules,data,num_spines,plottype):
    pyplot.ion()
    if len(plot_molecules)>8:
          rows=int(np.round(np.sqrt(len(plot_molecules))))
          cols=int(np.ceil(len(plot_molecules)/float(rows)))
    else:
          cols=1
          rows=len(plot_molecules)
    if plottype==3:
        fig=[]
        for sp in range(num_spines):
            f,a=pyplot.subplots(rows, cols,sharex=True) #n rows,  1 column #,figsize=(4*cols,rows)
            fig.append(f)
    else:
        fig,axes=pyplot.subplots(rows, cols,sharex=True)
    col_inc=[0.0,0.0]
    scale=['lin','lin']
    for i,paramset in enumerate(data.parlist):
        if len(paramset)>1: #color indexed by paramset
            col_inc[i]=(len(colors.colors)-1)/(len(paramset)-1)
            if plottype==2 and num_spines>1:
                col_inc[i]=(len(colors.colors)-1)/(len(paramset)*num_spines)
        elif plottype==2 and num_spines>1: #color indexed by spine number, if more than one spine (that are put onto same plot)
            col_inc[i]=(len(colors.colors)-1)/num_spines
        elif len(data.trials)>1 or (plottype==2 and num_spines<2) : #if only a single file with multiple trials and at most 1 spine, color indexed by trial
            #if only a single file, check/use the number of trials
            col_inc[i]=(len(colors.colors)-1)/(len(data.trials)-1)
        else:
            col_inc[i]=0.0 
            
    return fig,col_inc,scale

def get_color_label(parlist,params,colinc,parnames):
    if len(parlist[1])<len(parlist[0]):
        par_index=0
    else:
        par_index=1
    if len(parlist[1])==0 or len(parlist[0])==0:
        map_index=0
        color_index=int(parlist[par_index].index(params[par_index])*colinc[par_index]*partial_scale)
        mycolor=colors.colors[color_index]
    else:
        list_index=(par_index+1)%2
        map_index=parlist[list_index].index(params[list_index])
        color_index=int(parlist[par_index].index(params[par_index])*colinc[par_index]*partial_scale)
        mycolor=colors2D[map_index].__call__(color_index+offset[map_index])
    plotlabel='-'.join([parnames[ii]+str(k) for ii,k in enumerate(params)])
    return mycolor,plotlabel,par_index,map_index
 
def plotregions(plotmol,dataset,fig,colinc,scale,region_dict,textsize=12,regions=None):
    if not regions:
        regions=region_dict.keys()
        #call plot_setup, but pass len(region_dict.keys()) instead of number of stimulated spines
    num_regions=len(regions)
    axis=[]
    for reg in range(num_regions):
        axis.append(fig[reg].axes)

    for (fname,param) in dataset.ftuples:
        #First, determine the color scaling
        if len(dataset.ftuples)==1: 
            mycolor=[0,0,0]
            plotlabel=''
        else:
            mycolor,plotlabel,par_index,map_index=get_color_label(dataset.parlist,param,colinc,dataset.params)
        #Second, plot each molecule
        for imol,mol in enumerate(plotmol):
            maxpoint=min(len(dataset.time_set[param][mol]),np.shape(dataset.file_set_conc[param][mol])[1])
            for regnum,reg in enumerate(regions):
                #if len(dataset.ftuples)==1:  #NEED TO TEST THIS ON MORPHOLOGY WITHOUT SPINES
                #    for t in range(len(dataset.trials)):
                #        mycolor=colors.colors[int(colinc[0]*t)]
                #        axis[regnum][imol].plot(dataset.time_set[param][mol][0:maxpoint],dataset.means['regions'][param][mol][t,0:maxpoint].T[regnum],
                #                                   label=reg+' trial'+str(t),color=mycolor)
                #else:
                axis[regnum][imol].plot(dataset.time_set[param][mol][0:maxpoint],np.mean(dataset.file_set_conc[reg][param][mol][0:maxpoint],axis=0).T[regnum],
                                       label=plotlabel,color=mycolor)
                axis[regnum][imol].set_ylabel(mol+' (nM)',fontsize=textsize)
                axis[regnum][imol].tick_params(labelsize=textsize)
                axis[regnum][imol].set_xlabel('Time (sec)',fontsize=textsize)
        for regnum,reg in enumerate(regions):
            axis[regnum][imol].legend(fontsize=legtextsize, loc='best')

def spine_name(text,char=5):
    non_list=text.split('[')[-1].split(']')[0]
    return non_list[0:char]

def plottrace(plotmol,dataset,fig,colinc,scale,spinelist,plottype,textsize=12):
    num_spines=len(spinelist)
    print("plottrace: plotmol,parlist,parval:", plotmol,dataset.parlist,[p[1] for p in dataset.ftuples])

    if plottype==3:
        axis=[]
        for sp in range(num_spines):
            axis.append(fig[sp].axes)
    else:
        axis=fig.axes  #increments col index 1st
    print('***************','shape of axis for plottrace', np.shape(axis))

    if len(dataset.ftuples)==1:
        fname,param = dataset.ftuples[0]
        for imol,mol in enumerate(plotmol):
            maxpoint=min(len(dataset.time_set[param][mol]),np.shape(dataset.file_set_conc['Overall'][param][mol])[1])
            if plottype==1: # one graph, only plotting overall conc - plot individual trials
                for t in range(len(dataset.trials)):
                    mycolor=colors.colors[int(colinc[0]*t)]
                    axis[imol].plot(dataset.time_set[param][mol][0:maxpoint],dataset.file_set_conc['Overall'][param][mol][t,0:maxpoint],label='trial'+str(t),color=mycolor)
            elif plottype==2:  #one figure
                if num_spines>2: #plot spine conc, average across trials, color indexed by spine number
                    for spnum,sp in enumerate(spinelist):
                        new_col=colors.colors[int(colinc[0]*spnum)]
                        axis[imol].plot(dataset.time_set[param][mol][0:maxpoint],np.mean(dataset.file_set_conc[sp][param][mol][:,0:maxpoint],axis=0).T, label=spine_name(sp),color=new_col)
                else:        #plot individual trials, color indexed by trial and spine num
                    for spnum,sp in enumerate(spinelist):
                        map_index=spnum #either 0 or 1
                        for t in range(len(dataset.trials)):
                            new_index=int(colinc[0]*t*partial_scale)
                            mycolor=colors2D[map_index].__call__(new_index+offset[map_index])
                            axis[imol].plot(dataset.time_set[param][mol][0:maxpoint],dataset.file_set_conc[sp][param][mol][t,0:maxpoint].T,
                                                   label=spine_name(sp)+' trial'+str(t),color=mycolor)
            elif plottype==3:   #each spine & non-spine on separate figure, color indexed by trial
                for spnum,sp in enumerate(spinelist):
                    for t in range(len(dataset.trials)):
                        mycolor=colors.colors[int(colinc[0]*t)]
                        axis[spnum][imol].plot(dataset.time_set[param][mol][0:maxpoint],dataset.file_set_conc[sp][param][mol][t,0:maxpoint].T,
                                               label='trial'+str(t),color=mycolor)
                        axis[spnum][imol].set_ylabel(mol+' (nM)',fontsize=textsize)
                        axis[spnum][imol].tick_params(labelsize=textsize)
                        axis[spnum][imol].set_xlabel('Time (sec)',fontsize=textsize)
    else:
        for (fname,param) in dataset.ftuples:
            #First, determine the color scaling
            mycolor,plotlabel,par_index,map_index=get_color_label(dataset.parlist,param,colinc,dataset.params)
            #Second, plot each molecule
            for imol,mol in enumerate(plotmol):
                #axis[imol].autoscale(enable=True,tight=False)
                maxpoint=min(len(dataset.time_set[param][mol]),np.shape(dataset.file_set_conc['Overall'][param][mol])[1])
                if num_spines>1 and plottype==2: #index color by parameter and spine number if all spines on one plot
                    for spnum,sp in enumerate(spinelist):
                        new_index=int((dataset.parlist[par_index].index(param[par_index])*num_spines+spnum)*colinc[par_index]*partial_scale)
                        new_col=colors2D[map_index].__call__(new_index+offset[map_index]) #colors.colors[new_index] #
                        axis[imol].plot(dataset.time_set[param][mol][0:maxpoint],np.mean(dataset.file_set_conc[sp][param][mol][:,0:maxpoint],axis=0).T, label=plotlabel+' '+spine_name(sp),color=new_col)
                elif plottype==3: #index color by parameter if each spine on separate plot
                    for spnum,sp in enumerate(spinelist):
                        axis[spnum][imol].plot(dataset.time_set[param][mol][0:maxpoint],np.mean(dataset.file_set_conc[sp][param][mol][:,0:maxpoint],axis=0).T,
                                            label=plotlabel,color=mycolor)
                        axis[spnum][imol].set_ylabel(mol+' (nM)',fontsize=textsize)
                        axis[spnum][imol].tick_params(labelsize=textsize)
                        axis[spnum][imol].set_xlabel('Time (sec)',fontsize=textsize)
                else: #either plottype==1 (only plot overall) or plottype==2 and only one region (also plot Overall), then use param to index color, 
                    axis[imol].plot(dataset.time_set[param][mol][0:maxpoint],np.mean(dataset.file_set_conc['Overall'][param][mol][:,0:maxpoint],axis=0),label=plotlabel,color=mycolor)
    if plottype<3:
        for imol,mol in enumerate(plotmol):
            axis[imol].set_ylabel(mol+' (nM)',fontsize=textsize)
            axis[imol].tick_params(labelsize=textsize)
            axis[imol].set_xlabel('Time (sec)',fontsize=textsize)
    if plottype==3:
        for spnum,sp in enumerate(spinelist):
            axis[spnum][imol].legend(fontsize=legtextsize, loc='best')
    else:
        axis[imol].legend(fontsize=legtextsize, loc='best')
    #pyplot.tight_layout() #KEEP??????
    #fig.canvas.draw()
    return

def plotss(plot_mol,xparval,ss):
    fig,axes=pyplot.subplots(len(plot_mol), 1,sharex=True)
    for imol,mol in enumerate(plot_mol):
        axes[imol].plot(xparval,ss[:,imol],'o',label=mol)
        axes[imol].set_ylabel('nM')
        if max(xparval)/min(xparval)>100:
            axes[imol].set_xscale("log")
        axes[imol].legend(fontsize=legtextsize, loc='best',  fontweight='bold')
    fig.canvas.draw()
    return

def plot_total_mol(tot_species,dataset,figtitle,colinc,textsize,regions=None):
    numcols=len(tot_species)
    #Will need to specify whether plotting spine (and non-spine) totals, currently, regions are overall, dsm and spine head.  need non-spine (dendrite)
    if not regions:
        regions=dataset.file_set_tot.keys()
    numrows= len(regions) 
    fig,axes=pyplot.subplots(numrows,numcols,sharex=True)
    fig.canvas.manager.set_window_title(figtitle+'Totals')
    axis=fig.axes
    for row,region in enumerate(regions):
        for i,(param,total_trace) in enumerate(dataset.file_set_tot[region].items()):
            if len(dataset.file_set_tot[region].keys())==1:
                mycolor=[0,0,0]
                plotlabel=''
            else:
                mycolor,plotlabel,par_index,map_index=get_color_label(dataset.parlist,param,colinc,dataset.params)
            for col,(mol,trace) in enumerate(total_trace.items()): 
                #print('$$$$$$$$$$ pu.ps',param,mol,np.shape(trace),mycolor,plotlabel)
                ax=col+row*numcols
                newtime = np.linspace(0,dataset.endtime[param][mol], np.shape(trace)[1]) #convert from ms to sec
                axis[ax].plot(newtime,np.mean(trace,axis=0),label=plotlabel,color=mycolor)
                axis[col].set_title(mol+' TOTAL',fontsize=textsize)
                axis[-1].set_xlabel('Time (sec)',fontsize=textsize)
                axis[ax].tick_params(labelsize=textsize)
                axis[row*numcols].set_ylabel(region+' Conc (nM)',fontsize=textsize)
        axis[0].legend(fontsize=legtextsize, loc='best')#for now put legend into panel 0
    fig.canvas.draw()

def plot_signature(dataset,thresholds,figtitle,colinc,textsize):
    numcols=len(dataset.sig) #og.sig[mol][par][region]
    key0=list(dataset.sig.keys())[0]
    numrows= len(thresholds[key0].keys())
    fig,axes=pyplot.subplots(numrows,numcols,sharex=True)
    fig.canvas.manager.set_window_title(figtitle+' Signature')
    axis=fig.axes
    if len(dataset.ftuples)==1:
        fname,param = dataset.ftuples[0]
        for col,mol in enumerate(dataset.sig.keys()):
            for row,(region,trace) in enumerate(dataset.sig[mol][param].items()): 
                for t,trial in enumerate(dataset.trials):
                    mycolor=colors.colors[int(colinc[0]*t)]
                    ax=col+row*numcols
                    newtime = np.arange(np.shape(trace)[1])*dataset.dt[mol] #convert from ms to sec
                    axis[ax].plot(newtime,trace[t,:],label=trial,color=mycolor)
                axis[row*numcols].set_ylabel(region,fontsize=textsize)
            axis[col].set_title(mol+' SIGNATURE',fontsize=textsize)
            axis[ax].set_xlabel('Time (sec)',fontsize=textsize)
        axis[0].legend(fontsize=legtextsize, loc='best')#for now put legend into panel 0
    else:
        for col,mol in enumerate(dataset.sig.keys()):
            for i,(param,total_trace) in enumerate(dataset.sig[mol].items()):
                mycolor,plotlabel,par_index,map_index=get_color_label(dataset.parlist,param,colinc,dataset.params)
                for row,(region,trace) in enumerate(total_trace.items()): 
                    #print('$$$$$$$$$$ pu.ps',param,mol,np.shape(trace),mycolor,plotlabel)
                    ax=col+row*numcols
                    newtime = np.arange(np.shape(trace)[1])*dataset.dt[mol] #convert from ms to sec
                    axis[ax].plot(newtime,np.mean(trace,axis=0),label=plotlabel,color=mycolor)
                    axis[ax].tick_params(labelsize=textsize)
                    axis[row*numcols].set_ylabel(region+' Conc (nM)',fontsize=textsize)
            axis[col].set_title(mol+' SIGNATURE',fontsize=textsize)
            axis[ax].set_xlabel('Time (sec)',fontsize=textsize)
            axis[0].legend(fontsize=legtextsize, loc='best')#for now put legend into panel 0
    if len(thresholds[key0].keys()): 
        for col,mol in enumerate(thresholds.keys()):
            for param in dataset.sig[mol].keys():
                for row,thresh_val in enumerate(thresholds[mol].values()):
                    ax=col+row*numcols
                    axis[ax].plot([0,newtime[-1]],[thresh_val,thresh_val],color='gray',linestyle= 'dashed')
    fig.canvas.draw()

def tweak_fig(fig,yrange,legendloc,legendaxis,legtextsize):
     axes=fig.axes
     for axis in axes:
          axis.set_ylim(yrange)
          axis.set_ylim(yrange)
          axes[legendaxis].legend(fontsize=legtextsize, loc=legendloc)
     fig.tight_layout()

def axis_config(ax):
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.get_xaxis().set_tick_params(direction='out', right=0, which='both')
    ax.get_yaxis().set_tick_params(direction='out', top=0, which='both')

def axlabel(ax, label):
    ax.text(-0.2, 1.05, label, transform=ax.transAxes,
            fontweight='bold', va='top', ha='right')   

def plot3D(image,parval,figtitle,molecules,xvalues,time):
     from matplotlib.ticker import MultipleLocator
     minx=float(xvalues[0])
     maxx=float(xvalues[-1])
     asp=time[-1]/(maxx-minx)/len(parval) #may depend on number of subplots! 
     fig,axes=pyplot.subplots(len(parval),1,sharex=True,sharey=True,figsize=(6,9))
     fig.canvas.manager.set_window_title(figtitle)
     fig.suptitle('+'.join(molecules))
     for par in range(len(parval)):
          #for some reason, y axes are not correct without *10 in extent
          cax=axes[par].imshow(image[par].T,extent=[0,time[-1],minx,maxx],aspect=asp,vmin=0,vmax=np.max(image),origin='lower')
          if np.min(image[par])>=0:
               fig.colorbar(cax,ax=axes[par],ticks=MultipleLocator(round(np.max(image)/4)))
          axes[par].set_ylabel(parval[par])
          axes[par].set_xlabel('Time (sec)')

def pairs (dataset,mol_pairs,timeframe):
    for pair in mol_pairs:
        do_plot=(pair[0] in dataset.molecules) and (pair[1] in dataset.molecules) 
        if do_plot:
            if (dataset.dt[pair[0]]==dataset.dt[pair[1]]):
                plot_start=int(timeframe[0]/dataset.dt[pair[0]])
                plot_end=int(timeframe[1]/dataset.dt[pair[0]])
                pyplot.figure()
                pyplot.title('---'.join(pair))
                for pnum,(param,data) in enumerate(dataset.file_set_conc.items()):
                    labl='-'.join([str(k) for k in param])
                    X=np.mean(data[pair[0]],axis=0)
                    Y=np.mean(data[pair[1]],axis=0)
                    print('pairs plot',pnum,param,data.keys(),np.shape(X),np.shape(Y))
                    # check if molX & molY same length
                    if len(X)==len(Y):
                        pyplot.plot(X[plot_start:plot_end],Y[plot_start:plot_end], label=labl, linestyle='--')
                    else:
                        time_vectorY=dataset.time_set[param][pair[1]]
                        time_vectorX=dataset.time_set[param][pair[0]]
                        if len(X)>len(Y):
                            X=np.interp(time_vectorY,time_vectorX,X)
                        if len(Y)>len(X):
                            Y=np.interp(time_vectorX,time_vectorY,Y)
                        pyplot.plot(X[plot_start:plot_end],Y[plot_start:plot_end], label=labl, linestyle='--')
                pyplot.legend()
                pyplot.xlabel(pair[1])
                pyplot.ylabel(pair[0])
            else:
                print('******* Molecules', pair, 'in ARGS but saved with different dt **********', dataset.dt[pair[0]],dataset.dt[pair[0]])               
        else:
            print('******* Molecule not in ARGS****** or molecules saved with different dt **********', pair)
    
#from matplotlib.ticker import FuncFormatter
#def other_stuff():
     #PercentFormatter = FuncFormatter(lambda x, pos: '{:.0%}'.format(x).replace('%', r'\%'))
     #plt.rc('text', usetex=1)
     #plt.rc('text.latex', unicode=1)
     #plt.rc('font', weight='bold')
     #plt.rc('xtick', labelsize=20)
     #plt.rc('ytick', labelsize=20)
     #plt.rc('axes', labelsize=25)
     #plt.rc('axes', labelweight='bold')
     #plt.rc('legend', frameon=False)
     #plt.rc('legend', fontsize=20)
     #plt.rc('figure.subplot', bottom=0.15, left=0.18, right=0.93, top=0.93)
     #plt.rc('axes', color_cycle=['r', 'g', 'b', 'c', 'm', 'k', 'y'])
     #plt.rc('legend', numpoints=1)
     #matplotlib.rc('axes.formatter', useoffset=False)
