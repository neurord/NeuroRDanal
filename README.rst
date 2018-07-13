===========
NeuroRDanal
===========

**1. nrdh5_anal.py**
---------------------

Processes hdf5 output from NeuroRDv3 to produce graphs of molecules for one or more files, which use same morphology but can have different other parameters. To process multiple files, the first part of file name must be the same for all files, and parameter variations are specified as -par1-par2.  In addition, the set of files MUST use the same morphology file and MUST use the same discretization (specified in the top level model file).  If the set of files differ in morphology, then each file must be processed separately.
Graphs are generated for either a set of specified molecules, or all molecules if none are specified.  Basal value, peak and minimum are also printed, where basal value is calculated between two specified timepoints (basal_start basal_end).

First, download NeuroRDanal form the GitHub website.  Make sure they are in a subdirectory called "NeuroRDanal".  Set your python path to the directory containing NeuroRDanal.  For example, if NeuroRDanal is in the home directory of the user, then use the bash command: export $PYTHONPATH=$HOME.

To run the program from within python, type ARGS="subdir/fileroot,par1 par2,mol1 mol2,basal_start basal_end" then execfile('nrdh5_anal.py')
from outside python, type python nrdh5_anal "subdir/fileroot [par1 par2] [mol1 mol2] [basal_start basal_send]"
DO NOT PUT ANY SPACES NEXT TO THE COMMAS, DO NOT USE TABS
mol1 mol2, etc are the names of molecles to process
par1 and optionally par2 are used to construct filenames as "subdir/fileroot"+"-"+par1+"*"-"+par2+"*"
DO NOT use hyphens in filenames except for preceding parameter name
if no parameters specified, then fileroot needs to be full filename (excluding the .h5 extension)
If only a single file is specified, will plot multiple trials; if multiple files, plots the mean over trials for each file

Other parameters to adjust in program
1. outputavg - set to 1 to create region average output files to read into your favorite graphin software
2. showplot - set to 2 to plot the spine head concentration
3. stimspine - this should be the name of the spine head you want to plot with showplot=2
4. spinehead - this should be the name of your spinehead region

**2. sig.py**
---------------------
Calculate LTP/LTD signature from two sets of molecules, separately for spines and dendrites, by adding together the specified molecules, and then calculating area under the curve (and above the specified thresholds)

To run the program from within python, type ARGS="subdir/fileroot,par1 par2,LTPmol1 LTPmol2,LTDmol1 LTdmol2,basal_start basal_end, T_LTPd T_LTPsp T_LTDd T_LTDsp", then execfile('sig.py')
from outside python, type python sig.py "subdir/fileroot [par1 par2] [LTPmol1 LTPmol2] [LTDmol1 LTdmol2] [basal_start basal_end] [T_LTPd T_LTPsp T_LTDd T_LTDsp]"
LTPmol1 LTPmol2, etc are the names of molecles which produce LTP is sufficiently high (and hinder LTD)
LTDmol1 LTDmol2, etc are the names of molecles which produce LTD is sufficiently high (and hinder LTP)
T_LTPd T_LTPsp T_LTDd T_LTDsp are thresholds - defining "sufficiently high"

**3. plot_h5.py**
---------------------

Utilities used by nrdh5_anal.py and sig.py for creating graphs

**4. h5utils.py**
---------------------

Utilities used by nrdh5_anal.py and sig.py for creating region averages

**5. neurord_analysis.py**
---------------------------
Processes text file output from NeuroRDv3 to produce graphs of molecules for one or more files, which use same morphology but can have different other parameters. In other words, the set of files MUST use the same morphology file and MUST use the same discretization (specified in the top level model file).  If the set of files differ in morphology, then each file must be processed separately. Don't use this unless you can't get the hdf5 output to work. 

First, download NeuroRDanal form the GitHub website.  Make sure they are in a subdirectory called "NeuroRDanal".  Set your python path to the directory containing NeuroRDanal.  For example, if NeuroRDanal is in the home directory of the user, then use the bash command: export $PYTHONPATH=$HOME.

To run the program from within python, type ARGS="subdir/fileroot,par1 par2,mol1 mol2,basal_start basal_end" then execfile('neurord_analysis.py')
from outside python, type python neurord_analysis "subdir/fileroot [par1 par2] [mol1 mol2] [basal_start basal_end]"
DO NOT PUT ANY SPACES NEXT TO THE COMMAS, DO NOT USE TABS
mol1 mol2, etc are the names of molecles to process
par1 and optionally par2 are used to construct filenames as "subdir/fileroot"+"-"+par1+"*"-"+par2+"*"
basal_start and basal_end are the time, in seconds, prior to stimulation to use for calculating basal values
DO NOT use hyphens in filenames except for preceding parameter name
if no parameters specified, then fileroot needs to be full filename (excluding the .txt extension)

**6.header_parse.py**
---------------------
Utilities used by neurord_analysis for reading the first header line and determining which columns of data belong to which molecule, which voxel, and which region of the morphology.

**7. plot_utils.py**
--------------------
Utilities used by neurord_analysis for plotting the NeuroRD output

**8. sig2.py**
---------------
Program to read in the text file outputs of sig.py and generate a file of molecule-space-time samples - one line per file - for statistical analysis.  Alternatively, generate signature traces (normalized sum of a subset of the molecules) and plot them.
