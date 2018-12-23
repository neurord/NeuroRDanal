from __future__ import print_function
from __future__ import division

import os
import numpy as np
from matplotlib import pyplot
import sys  
import glob
import h5utils
import plot_h5 as pu5
import h5py as h5
from scipy import optimize

Avogadro=6.023e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15
msec_per_sec=1000

window_size=1  #number of seconds on either side of peak value to average for maximum
showplot=1    #2 indicates plot the head conc, 0 means no plots
textsize=10

try:
    args = ARGS.split(",")
    print("ARGS =", ARGS, "commandline=", args)
    do_exit = False
except NameError: #NameError refers to an undefined variable (in this case ARGS)
    args = sys.argv[1:]
    print("commandline =", args)
    do_exit = True

ftuples,parlist,params=h5utils.argparse(args)
figtitle=args[0].split('/')[-1]+args[1]
input_molecules=args[2].split()
output_molecules=args[3].split()
plot_molecules=input_molecules+output_molecules

try:
    data.close()
except Exception:
    pass

parval=[]
numfiles=len(ftuples)
whole_plot_array=[]
whole_time_array=[]
for fnum,ftuple in enumerate(ftuples):
    data,maxvols,TotVol,trials,seeds,arraysize,p=h5utils.initialize(ftuple,numfiles,parval)
    print('done')
    plot_array=[]
    time_array=[]
    if fnum==0:
        molecules=h5utils.decode(data['model']['species'][:])
        #initialize plot stuff, arrays for static measurements, and find molecules among the output sets
        if numfiles>1:
            whole_plot_array=[[] for mol in plot_molecules]
            whole_time_array=[[] for mol in plot_molecules]
        num_mols=len(plot_molecules)
        baseline=np.zeros((arraysize,num_mols))
        auc=np.zeros((arraysize,num_mols))
        peakval=np.zeros((arraysize,num_mols))
    #which columns are needed for each molecule
    out_location,dt,rows=h5utils.get_mol_info(data,plot_molecules,maxvols)
    #
    #Which "rows" should be used for baseline value, specifed in args[3].  If different for each file then will have problems later
    sstart,ssend=h5utils.sstart_end(plot_molecules,args,4,out_location,dt,rows)
    voxel=0
    for mol in plot_molecules:
      if out_location[mol]!=-1:
        outset = list(out_location[mol]['location'].keys())[0]
        imol=out_location[mol]['location'][outset]['mol_index']
        tempConc=np.zeros((len(trials),out_location[mol]['samples']))
        time_array.append(data[trials[0]]['output'][outset]['times'][:]/msec_per_sec)
        #generate output files for these cases
        for trialnum,trial in enumerate(trials):
            tempConc[trialnum]=data[trial]['output'][outset]['population'][:,voxel,imol]/TotVol/mol_per_nM_u3
        plot_array.append(np.mean(tempConc,axis=0))
             #plot_array dimensions=numfiles x number of molecules x sample times
      else:
          if fnum==0 and molecule_name_issue==0:
              print("Choose molecules from:", molecules)
              molecule_name_issue=1
          time_array.append([])
          plot_array.append([])
    if numfiles>1:
        #plot_array dimensions=num molecules x sample times
        #whole_plot_array dimension  = num molecules*num files*sample time
        for mol in range(num_mols):
            whole_plot_array[mol].append(plot_array[mol])
            whole_time_array[mol].append(time_array[mol])

def peakval(pnum,plot_molecules,mol,sstart,ssend):
    imol=plot_molecules.index(mol)
    if out_location[mol]!=-1:
        dt=out_location[mol]['dt']
        window=int(window_size/dt)
        baseline=whole_plot_array[imol][pnum][sstart[imol]:ssend[imol]].mean()
        peakpt=whole_plot_array[imol][pnum][ssend[imol]:].argmax()+ssend[imol]
        #auc[pnum,imol]=np.sum(whole_plot_array[imol][pnum][ssend[imol]:endpt]-baseline[pnum,imol])*dt[imol]
        peakval=whole_plot_array[imol][pnum][peakpt-window:peakpt+window].mean()
    return peakval

def f1(conc,Kd,N):
    return conc**N/(conc**N+Kd**N)

def f2(conc,Kd,N):
    return (conc/(conc+Kd))**N

################# This part is new  compare to sig.py and nrdh5_anal.py
#reshape and print
inputvals=np.zeros((len(parlist[0]),len(parlist[1])))
outputvals=np.zeros((len(output_molecules),len(parlist[0]),len(parlist[1])))
print(params[0],' ',params[1],"    peakvals   ")
in_mol=input_molecules[0]
for p1,par1 in enumerate(parlist[0]):
    for p2,par2 in enumerate(parlist[1]):
        par_index=parval.index((par1,par2))
        print (' ',par1,' ',par2, ' ', in_mol,end=': ')
        inputvals[p1,p2]=peakval(par_index,plot_molecules,in_mol,sstart,ssend)
        print("%8.2f"%(inputvals[p1,p2]),end=', ')
        for imol,mol in enumerate(output_molecules):
            outputvals[imol,p1,p2]=peakval(par_index,plot_molecules,mol,sstart,ssend)
        print(outputvals[:,p1,p2])
                             
#  fit outputvalues vs inputvalues to sigmoid
pyplot.ion()
pyplot.figure()
fit1={}
markers=['ro','b*','kv','cs']
#fit2={}
for imol,mol in enumerate(output_molecules):
    fit1[mol]={}
    for p1,par1 in enumerate(parlist[0]):
        popt,pcov=optimize.curve_fit(f1,inputvals[p1,:],outputvals[imol,p1,:])
        fit1[mol][par1]={'Kd':popt[0],'N':popt[1]}
        popt,pcov=optimize.curve_fit(f2,inputvals[0,:],outputvals[0,0,:])
        fit2[par1]={'Kd':popt[0],'N':popt[1]}
        pyplot.plot(inputvals[p1,:],outputvals[imol,p1,:], markers[p1+imol*len(parlist[0])],label=mol+'-'+par1)
pyplot.legend()
print(fit1)
'''
if showplot:
    fig,col_inc,scale=pu5.plot_setup(plot_molecules,parlist,params,0,showplot)
    #need fnames
    fig.canvas.set_window_title(figtitle)
    pu5.plottrace(plot_molecules,whole_time_array,whole_plot_array,parval,fig,col_inc,scale,parlist,textsize,[],showplot)
    #
'''
    
