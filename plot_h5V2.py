from __future__ import print_function
from __future__ import division
import numpy as np
from matplotlib import pyplot

legtextsize=15
     
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
    for imol,mol in enumerate(dataset.molecules+dataset.tot_species):
        pyplot.figure() #new figure panel for each molecules
        pyplot.suptitle(title)
        if len(p): #reshape feature values for plotting if 2 params and all combos
            labels=p['labels'];xval_index=p['xindex'];label_index=p['lindex']
            new_yvals=np.zeros((len(xvals),len(labels)))
            for fnum,(fname,param) in enumerate(dataset.ftuples):
                row=xvals.index(param[xval_index])
                col=labels.index(param[label_index])
                new_yvals[row,col]=dataset.mean_feature[feature][imol,fnum]
            for col,label in enumerate(labels):
                pyplot.scatter(xvals,new_yvals[:,col],label=dataset.params[label_index]+' '+label)
            pyplot.xlabel(dataset.params[xval_index])
        else: #just plot feature values vs param or param combo            
            pyplot.scatter(xvals,dataset.mean_feature[feature][imol,:], label=mol)
            pyplot.xlabel('-'.join(dataset.params))
        pyplot.ylabel(mol+' '+feature)
        pyplot.legend()

def spatial_plot(data,dataset):
    numtrials=len(data.trials)
    for (fname,param) in dataset.ftuples:
        fig,axes=pyplot.subplots(len(data.molecules),numtrials, sharex=True)
        fig.suptitle(fname)
        for imol,mol in enumerate(data.molecules):
            for trial in range(numtrials):
               axes[imol,trial].imshow(dataset.spatial_means[param][mol][trial].T,extent=[0, np.max(dataset.time_set[param][mol]), float(list(data.spatial_dict.keys())[0]), float(list(data.spatial_dict.keys())[-1])],aspect='auto',origin='lower')
               #axes[imol,trial].colorbar()
            axes[imol,0].set_ylabel (mol +', location (um)')
        for trial in range(numtrials):
            axes[imol,trial].set_xlabel('time (ms)')
             
def plot_setup(plot_molecules,data,num_spines,plottype):
     pyplot.ion()
     if len(plot_molecules)>8:
          rows=int(np.round(np.sqrt(len(plot_molecules))))
          cols=int(np.ceil(len(plot_molecules)/float(rows)))
     else:
          cols=1
          rows=len(plot_molecules)
     fig,axes=pyplot.subplots(rows, cols,sharex=True) #n rows,  1 column #,figsize=(4*cols,rows)
     col_inc=[0.0,0.0]
     scale=['lin','lin']
     for i,paramset in enumerate(data.parlist):
          if len(paramset)>1:
               col_inc[i]=(len(colors.colors)-1)/(len(paramset)-1)
               if plottype==2 and num_spines>1:
                    col_inc[i]=(len(colors.colors)-1)/(len(paramset)*num_spines)
          else:
               col_inc[i]=0.0
     return fig,col_inc,scale

def get_color_label(parlist,params,colinc):
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
    plotlabel='-'.join([str(k) for k in params])
    return mycolor,plotlabel,par_index,map_index
 
def plottrace(plotmol,dataset,fig,colinc,scale,stimspines,plottype,textsize=12):
     num_spines=len(stimspines)
     print("plottrace: plotmol,parlist,parval:", plotmol,dataset.parlist,[p[1] for p in dataset.ftuples])
     axis=fig.axes  #increments col index 1st
     for (fname,param) in dataset.ftuples:
        #First, determine the color scaling
        if len(dataset.parlist)==0:
             mycolor=[0,0,0]
             plotlabel=''
        else:
            mycolor,plotlabel,par_index,map_index=get_color_label(dataset.parlist,param,colinc)
       #Second, plot each molecule
        for imol,mol in enumerate(plotmol):
            #axis[imol].autoscale(enable=True,tight=False)
            maxpoint=min(len(dataset.time_set[param][mol]),np.shape(dataset.file_set_conc[param][mol])[1])
            if num_spines>1 and plottype==2:
                for spnum,sp in enumerate(stimspines):
                     new_index=int((dataset.parlist[par_index].index(param[par_index])*num_spines+spnum)*colinc[par_index]*partial_scale)
                     new_col=colors2D[map_index].__call__(new_index+offset[map_index]) #colors.colors[new_index] #
                     axis[imol].plot(dataset.time_set[param][mol][0:maxpoint],np.mean(dataset.spine_means[param][mol][0:maxpoint],axis=0).T[spnum],
                                     label=plotlabel+sp.split('[')[-1][0:-1],color=new_col)
            else:
                 axis[imol].plot(dataset.time_set[param][mol][0:maxpoint],np.mean(dataset.file_set_conc[param][mol][0:maxpoint],axis=0),label=plotlabel,color=mycolor)
            axis[imol].set_ylabel(mol+' (nM)',fontsize=textsize, fontweight='bold')
            axis[imol].tick_params(labelsize=textsize)
            axis[imol].set_xlabel('Time (sec)',fontsize=textsize, fontweight='bold')
     axis[imol].legend(fontsize=legtextsize, loc='best')
     fig.canvas.draw()
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

def plot_signature(tot_species,dataset,figtitle,colinc,textsize,thresholds=[]):
    numcols=len(tot_species)
    #Will need to specify whether plotting spine (and non-spine) totals, or perhaps region totals
    #add additional total arrays in nrd_group
    numrows= 1 #one row for each region  For now, there is only an overall sum
    row=0 #once multiple rows, will loop over rows, and region as y label, mol as column heading?
    fig,axes=pyplot.subplots(numrows,numcols,sharex=True)
    fig.canvas.set_window_title(figtitle+'Totals')
    axis=fig.axes
    for i,(param,ss_tot) in enumerate(dataset.file_set_tot.items()):
        mycolor,plotlabel,par_index,map_index=get_color_label(dataset.parlist,param,colinc)
        for col,(mol,trace) in enumerate(ss_tot.items()): 
            #print('$$$$$$$$$$ pu.ps',param,mol,np.shape(trace),mycolor,plotlabel)
            ax=col+row*numrows
            newtime = np.linspace(0,dataset.endtime[param][mol], np.shape(trace)[1])/ms_to_sec #convert from ms to sec
            axis[ax].plot(newtime,np.mean(trace,axis=0),label=plotlabel,color=mycolor)
            axis[ax].set_title(mol+' TOTAL',fontsize=textsize)
            axis[ax].set_xlabel('Time (sec)',fontsize=textsize)
            axis[ax].tick_params(labelsize=textsize)
            axis[row*numrows].set_ylabel('Conc (nM)',fontsize=textsize)
    axis[0].legend(fontsize=legtextsize, loc='best')
    if len(thresholds): #this needs to be fixed.  Determine how to match thresholds to mol/sig
        r=(1,0)[row==0]
        axis[ax].plot([0,newtime[-1]],[thresholds[r],thresholds[r]],color='k',linestyle= 'dashed')
        axis[ax+1].plot([0,newtime[-1]],[thresholds[r+numcols],thresholds[r+numcols]],color='k',linestyle= 'dashed')
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
     fig.canvas.set_window_title(figtitle)
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
