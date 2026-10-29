######### Read in json files from sscalculation, and  convert to region concentrations
import numpy as np
from nrd_output import nrdh5_output
from nrd_group import nrdh5_group
import plot_h5V2 as pu5
from matplotlib import pyplot as plt
plt.ion()
colors=plt.get_cmap('plasma')

def conc_difference(ctrl,results):
	####### Calculate difference between  control and param variation
	delta={mol:{} for mol in ctrl.keys()}
	for mol_reg in ctrl.keys():
		delta[mol_reg]=results[mol_reg]-ctrl[mol_reg]
	return delta

def hist_plot(df,variables,numbins=30):
	from matplotlib import pyplot as plt
	plt.ion()
	fig,ax=plt.subplots(len(variables),1)
	ax=fig.axes
	for i,var in enumerate(variables):
		x=np.abs(df[df[var]!=0][var]) #eliminate negative numbers
		zeros=df[df[var]==0]
		logbins = np.logspace(np.log10(x.min()),np.log10(x.max()),numbins)
		# histogram on linear scale
		ax[i].hist(x, bins=logbins)
		ax[i].set_xlabel(var+', num zeros='+str(len(zeros)))
		ax[i].set_xscale('log')
		ax[i].set_ylabel('frequency')
	fig.tight_layout()
	return fig

def flatten(input):
	output=[]
	for ls in input:
		if isinstance(ls,list):
			for x in ls:
				output.append(x)
		else:
			output.append(ls)
	return output

def key_regions_dict(key_regions,data):
	means_loc={}
	for reg in key_regions:
		if reg in data.region_list:
			means_loc[reg]=('region',list(data.region_list).index(reg)) 
			if data.region_list[0]=='default':
				means_loc[reg]=('region',list(data.region_list).index(reg)-1)#-1 is adhoc - the region_list has 'default', which doesn't appear in data.means
		elif reg in data.spinelist:
			means_loc[reg]=('spines',list(data.spinelist).index(reg))
	return means_loc

def calc_ss_ratio(outmol,data,key_regions,results):
	for numer in outmol.keys():
		if len(outmol[numer]): #if there are denominator molecules
			denom=0
			for dmol in outmol[numer]:
				denom+=np.mean(data.OverallMean[dmol][:,data.sstart[dmol]:])
		else:
			denom=1
		results[numer+'_'+'Overall']=np.mean(data.OverallMean[numer][:,data.sstart[dmol]:])/denom #ratio won't work if loop over other mols
		#kluge: specifying two out of multiple structures
		for reg in key_regions:
			if len(outmol[numer]): #if there are denominator molecules
				denom=0
				for dmol in outmol[numer]:
					denom+=np.mean(data.means[means_loc[reg][0]][dmol][:,data.sstart[dmol]:data.ssend[dmol],means_loc[reg][1]])
			else:
				denom=1
			results[numer+'_'+reg]=np.mean(data.means[means_loc[reg][0]][numer][:,data.sstart[dmol]:,means_loc[reg][1]])/denom
	return results

def fig_setup(data,colors,showplot):
	colinc=(len(colors.colors)-1)/(len(data.trials))
	if showplot[0]=='1': # one graph, only plotting overall conc - plot individual trials
		fig,ax=plt.subplots(len(data.molecules),1, sharex=True)
	elif showplot[0]=='3':   #each spine & non-spine on separate figure, color indexed by trial
		fig,ax=plt.subplots(len(data.molecules),len(data.spinelist), sharex=True)
	else:
		fig,ax=None
		colinc=1
		print('unrecognized showplot value')
	return fig,ax,colinc

def plot_ctrl(data,colors,showplot,figtitle):
	fig,ax,colinc=fig_setup(data,colors,showplot)
	if showplot[0]=='1': # one graph, only plotting overall conc - plot individual trials
		for imol,mol in enumerate(data.molecules):
			for t in range(len(data.trials)):
				mycolor=colors.colors[int(colinc*t)]
				ax[imol].plot(data.time[mol],data.OverallMean[mol][t,:],label='trial'+str(t),color=mycolor)
			ax[imol].plot([0,np.max(data.time[mol])],[np.mean(data.OverallMean[mol])]*2,color='gray')
			ax[imol].set_ylabel(mol+' (nM)')
		ax[imol].set_xlabel('Time (sec)')
		fig.canvas.manager.set_window_title(figtitle)
	elif showplot[0]=='3':   #each spine & non-spine on separate figure, color indexed by trial
		for spnum,sp in enumerate(data.spinelist):
			for imol,mol in enumerate(data.molecules):
				for t in range(len(data.trials)):
					mycolor=colors.colors[int(colinc*t)]
					ax[imol][spnum].plot(data.time[mol],data.means['spines'][sp][mol][t,:].T,
											label='trial'+str(t),color=mycolor)
				ax[imol][0].set_ylabel(mol+' (nM)')
				ax[imol][spnum].tick_params()
			ax[-1][spnum].set_xlabel('Time (sec)')
			fig[spnum].suptitle(figtitle+' '+sp)
	return fig

#replace with arg_parser?
outmol={'pCof':['Cof','Cofactin'], 'pLIMK':[]} #molecule(s) to save as representative of change in basal.  Each key is numerator, each value is denominator, which can be empty list
all_mols=flatten(list(outmol.keys())+list(outmol.values()))
dirname='C:/Users/klblackwell/Documents/Research/cofilin/robust/' #directory with random param variation files 
ctrlfile='Model_Cof-ctrl-4trbasal' 
fileroot='Model_Cof-basalctrl' #base name of h5 file  
par=['random']
showplot='10'
#obtain steady state from last 100 sec of simulation
start=3500
end=3600

#name of regions and structures
dendname="dend" #from morph file
spinehead="head" #from morph file
submembname=dendname+'sub'
key_regions=['dend','sa1[0]']

############################################################
######### Read in and process control file, create means_loc - specifies index in array to find molecules
ftuple=(dirname+ctrlfile+'.h5','ctrl')
data=nrdh5_output(ftuple)
data.rows_and_columns(all_mols,[start,end])
data.molecule_population()
if data.maxvols>1:
	data.region_structures(dendname,submembname,spinehead)#,stimspine) #stimspine is optionald
	data.average_over_voxels(['ctrl'])
	means_loc=key_regions_dict(key_regions,data)

ctrl_results={} #dictionary to store results from all parameter variations
ctrl_results=calc_ss_ratio(outmol,data,key_regions,ctrl_results)
#plot data to ensure correct
if int(showplot[0])>0:
	fig=plot_ctrl(data,colors,showplot,ctrlfile)

############## Analyze each random file in the specified directory
og=nrdh5_group(dirname+fileroot,par)
results={f[1]:{} for f in og.ftuples} #dictionary to store results from all parameter variations
delta_dict={f[1]:{} for f in og.ftuples}

#for fname in allfiles: #loop over all files
for fnum,ftuple in enumerate(og.ftuples):
	#if 'random' in fname:
	#	index=fname.split('.h5')[0].split('random')[-1]
	#else:
	#	index='ctrl'
	#ftuple=(dirname+h5_file,(index,))
	data=nrdh5_output(ftuple)
	data.rows_and_columns(all_mols,[start,end])
	data.molecule_population()
	if data.maxvols>1:
		data.region_structures(dendname,submembname,spinehead)#,stimspine) #stimspine is optionald
		data.average_over_voxels(og.params)

	#data.total_subspecies(tot_species,sub_species,params.start,weights=weight)
	results[ftuple[1]]=calc_ss_ratio(outmol,data,key_regions,results[ftuple[1]])

	########### Calculate difference to control ###############
	delta_dict[ftuple[1]]=conc_difference(ctrl_results,results[ftuple[1]])
	#create dictionary for output which contains both the difference

	if int(showplot[0])>0:
		fig=plot_ctrl(data,colors,showplot,par[0]+str(int(ftuple[-1][0])))

#########Output to file
import pandas as pd
#create dataframe for output and figures
df=pd.DataFrame(delta_dict).transpose()	
df_fig=hist_plot(df,['pCof_Overall','pCof_sa1[0]','pLIMK_Overall','pLIMK_sa1[0]'])
df_fig.suptitle('absolute differences')
outfile=dirname+fileroot+'_compare.csv'
df.to_csv(outfile, index_label='molecule', float_format='%.4f')
df = pd.DataFrame(results).transpose()
outfile=dirname+'robust_output.csv'
df.to_csv(outfile, index_label='index', float_format='%.4f')

