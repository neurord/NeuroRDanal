from __future__ import print_function
from __future__ import division
import numpy as np
from matplotlib import pyplot

legtextsize=8
     
colors=pyplot.get_cmap('viridis')
#colors=pyplot.get_cmap('plasma')
colors2D=[pyplot.get_cmap('gist_heat'),pyplot.get_cmap('summer'),pyplot.get_cmap('Blues')]
offset=[0,0,63]  #avoid the light colors in low indices for the 'Blues' map
partial_scale=0.75 #avoid the very light colors.  Note that non-zero offset must be <= (1-partial_scale)*255

def space_avg(plot_molecules,whole_space_array,whole_time_array,trials,spatial_dict):
     for imol in range(len(plot_molecules)):
          fig,axes=pyplot.subplots(np.shape(whole_space_array[imol])[0], 1,sharex=True)
          fig.suptitle(plot_molecules[imol])
          for trial,tname in enumerate(trials):
               pyplot.subplot(np.shape(whole_space_array[imol])[0],1,trial+1)
               pyplot.imshow(whole_space_array[imol][trial].T,extent=[0, np.max(whole_time_array[imol][trial]), float(spatial_dict.keys()[0]), float(spatial_dict.keys()[-1])],aspect='auto',origin='lower')
               pyplot.colorbar()
               pyplot.ylabel (str(tname)+' (um)')
          pyplot.xlabel('time (ms)')
             
def plot_setup(plot_molecules,param_list,param_name,num_spines,plottype):
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
     for i,paramset in enumerate(param_list):
          if len(paramset)>1:
               col_inc[i]=(len(colors.colors)-1)/(len(paramset)-1)
               if plottype==2 and num_spines>1:
                    col_inc[i]=(len(colors.colors)-1)/(len(paramset)*num_spines-1)
          else:
               col_inc[i]=0.0
     return fig,col_inc,scale

def plottrace(plotmol,timearray,plotarray,parval,fig,colinc,scale,parlist,textsize,stimspines,plottype):
     num_spines=len(stimspines)
     print("plottrace: plotmol,parval,parlist:", plotmol,parval, parlist)
     axis=fig.axes  #increments col index 1st
     for pnum in range(len(parval)):
          #First, determine the color scaling
          if len(parlist)==0:
               mycolor=[0,0,0]
               plotlabel=''
          else:
               if np.shape(parlist[1])[0]==0 or np.shape(parlist[0])[0]==0:
                    if np.shape(parlist[1])[0]==0:
                         par_index=0
                    else:
                         par_index=1
                    color_index=int(parlist[par_index].index(parval[pnum])*colinc[par_index]*partial_scale)
                    mycolor=colors.colors[color_index]
                    plotlabel=parval[pnum]
               else:
                    if len(parlist[1])<len(parlist[0]):
                         par_index=0
                         map_index=parlist[1].index(parval[pnum][1])
                    else:
                         par_index=1
                         map_index=parlist[0].index(parval[pnum][0])
                    color_index=int(parlist[par_index].index(parval[pnum][par_index])*colinc[par_index]*partial_scale)
                    mycolor=colors2D[map_index].__call__(color_index+offset[map_index])
                    plotlabel=parval[pnum][0]+'-'+parval[pnum][1]
          #Second, plot each molecule
          for imol in range(len(plotmol)):
               #axis[imol].autoscale(enable=True,tight=False)
               #change label back to parval[pnum] after figures created
               maxpoint=min(len(timearray[imol][pnum][:]),len(plotarray[imol][pnum][:]))
               if num_spines>1 and plottype==2:
                    for spnum,sp in enumerate(stimspines):
                         new_index=int((parlist[par_index].index(parval[pnum])*num_spines+spnum)*colinc[par_index]*partial_scale)
                         #Note: this will not give expected color if 2 dimensions of parameters
                         new_col=colors.colors[new_index]
                         axis[imol].plot(timearray[imol][pnum][0:maxpoint],plotarray[imol][pnum][0:maxpoint].T[spnum],
                                         label=plotlabel+sp.split('[')[-1][0:-1],color=new_col)
               else:
                    axis[imol].plot(timearray[imol][pnum][0:maxpoint],plotarray[imol][pnum][0:maxpoint],label=plotlabel,color=mycolor)
               axis[imol].set_ylabel(plotmol[imol]+' (nM)',fontsize=textsize)
               axis[imol].tick_params(labelsize=textsize)
          axis[imol].set_xlabel('Time (sec)',fontsize=textsize)
     axis[0].legend(fontsize=legtextsize, loc='best')
     fig.canvas.draw()
     return

def plotss(plot_mol,xparval,ss):
    fig,axes=pyplot.subplots(len(plot_mol), 1,sharex=True)
    for imol,mol in enumerate(plot_mol):
        axes[imol].plot(xparval,ss[:,imol],'o',label=mol)
        axes[imol].set_ylabel('nM')
        if max(xparval)/min(xparval)>100:
            axes[imol].set_xscale("log")
        axes[imol].legend(fontsize=legtextsize, loc='best')
    fig.canvas.draw()
    return

def plot_signature(condition,traces,dt,figtitle,sign_title,textsize,thresholds,moretraces=[]):
     if len(moretraces):
          plot_ltd=1
     else:
          plot_ltd=0
     if len(np.shape(traces))==3:
          numrows=np.shape(traces)[2]
     else:
          numrows=1
     if len(moretraces):
          numcols=2
     else:
          numcols=1
     fig,axes=pyplot.subplots(numrows,numcols,sharex=True,figsize=(4,3))
     axis=fig.axes
     if numrows==1:
          if len(condition)>1:
               colinc=len(colors)/(len(condition)-1)
          else:
               colinc=0
          for i,cond in enumerate(condition):
               numpoints=np.shape(traces[i])[0]
               newtime = np.linspace(0,dt*(numpoints-1), numpoints+1)
               axis[0].plot(newtime,traces[i],label=cond,color=colors[int(i*colinc)])
               axis[0].legend(fontsize=legtextsize, loc='best')
               axis[0].set_ylabel('LTP sig (nM) ',fontsize=textsize)
               axis[0].set_xlabel('Time (sec)',fontsize=textsize)
               axis[0].tick_params(labelsize=textsize)
               if plot_ltd:
                    axis[1].plot(newtime,moretraces[i],color=colors[int(i*colinc)])
                    axis[1].set_ylabel('LTD sig (nM) ',fontsize=textsize)
                    axis[1].set_xlabel('Time (sec)',fontsize=textsize)
     else:
          domain=[]
          num_par=len(condition)/2
          for row,j in enumerate(range(0,numrows*numcols,numcols)):
               #thresh_index=np.ceil(row/numrows)
               domain.append(condition[0][row].split()[-1])
               for i,cond in enumerate(condition):
                    #the following assumes that first parameter has only two values
                    map_index=int( i/num_par )
                    color_index=int( i%num_par *(255/num_par) )
                    numpoints=np.shape(traces[i])[0]
                    newtime = np.linspace(0,dt*(numpoints-1), numpoints)
                    labl=cond[0][0:cond[0].rfind(' ')]
                    if j==0:
                         axis[j].plot(newtime,traces[i,:,row],label=labl,color=colors2D[map_index].__call__(color_index))
                    else:
                         axis[j].plot(newtime,traces[i,:,row],color=colors2D[map_index].__call__(color_index))
                    if plot_ltd:
                         axis[j+1].plot(newtime,moretraces[i,:,row],color=colors2D[map_index].__call__(color_index))
               axis[0].legend(fontsize=legtextsize, loc='lower right')
               axis[j].set_ylabel(domain[row],fontsize=textsize)
               axis[j].tick_params(labelsize=textsize)
               r=(1,0)[row==0]
               axis[j].plot([0,newtime[-1]],[thresholds[r],thresholds[r]],color='k',linestyle= 'dashed')
               if plot_ltd:
                    axis[j+1].plot([0,newtime[-1]],[thresholds[r+numcols],thresholds[r+numcols]],color='k',linestyle= 'dashed')
          if plot_ltd:
               axis[len(axis)-1].set_xlabel('Time (sec) LTD',fontsize=textsize)
               axis[len(axis)-2].set_xlabel('Time (sec) LTP',fontsize=textsize)
          else:
               axis[len(axis)-1].set_xlabel('Time (sec)',fontsize=textsize)
     fig.canvas.set_window_title(figtitle)
     fig.suptitle(sign_title)
     fig.canvas.draw()
     return

#The files are sorted on the parameters, assumes no hyphens in root filename
def file_tuple(fnames,params):
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
          print('pu: fname, par0:',fname,par0list)
          if len(params)>1:
               parval1=part_fname.split('-'+params[1])[-1] 
               ftuple.append((fname,(parval0,parval1)))
               if (parval1 not in par1list):
                    par1list.append(parval1)
               print('pu: par1:',par1list)
          else:
               ftuple.append((fname,parval0))
     return ftuple,[par0list,par1list]

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
