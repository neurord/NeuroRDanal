from __future__ import print_function
from __future__ import division
import numpy as np
Avogadro=6.023e14 #to convert to nanoMoles
mol_per_nM_u3=Avogadro*1e-15

def header_parse(header,maxvols,prninfo):
    #Define some substrings used in reading header
    underscore='_'
    dot='.'
    #
    #number of voxels in the file
    volnum=np.zeros(maxvols,dtype=int)
    regionID=list()
    structType=list()
    molecules=list()
    #
    volnames=header.split()
    for i in range(maxvols):
        #print "volname",volnames[i+1]
        voxel_text=volnames[i+1].split(underscore)
        #voxel number is 2nd "phrase"
        volnum[i]=int(voxel_text[1])
        if i==0:
            minvolnum=int(volnum[0])
            print("min vox num", minvolnum,"max vox num", maxvols)
        #if output of only part of morphology, need to reduce maxvols
        if (i>0 and volnum[i]==volnum[0]):
            maxvols=int(volnum[i-1]-minvolnum+1)
            print("*******new max", maxvols)
            #need to exit loop
        else:
            dotyn=voxel_text[2].find(dot)
            if dotyn==-1:
                #region ID, such as dendrite, is 3d "phrase" if not spine
                #for shorter heading, only use first 3 char of region ID
                regionID.append(voxel_text[2][0:4])
            else:
                regionID.append(voxel_text[2].split(dot)[0])
            #sm vs cyt is 4th "phrase"
            structType.append(voxel_text[3][0:3])
    volnums=np.resize(volnum,maxvols)
    #These two lines not needed once loop exited
    regionID=np.resize(regionID,maxvols)
    structType=np.resize(structType,maxvols)
    if prninfo:
        print("in hparse: voxels", volnums)
        print("   regionID", np.size(regionID),regionID)
        print("   structures", np.size(structType),structType)
    #
    #Now determine the molecules in the file - last characters of volume string
    #skip all[0] which has time, increment by maxvols to find vol0 of each molecule
    for i,j in enumerate(range(1,len(volnames),maxvols)):
        molecules.append(volnames[j].split(underscore)[-1])
        if prninfo:
            print(j, volnames[j],molecules[i])
    return regionID,structType,molecules,volnums,maxvols

def read_mesh(meshname,prninfo):
    volcol=13
    #the x and y data can be used to create spatial averages over the dendrite
    #sum of volumes within a region will be used for calculating mean concentration
    meshdata=np.loadtxt(meshname,skiprows=1)
    if np.shape(meshdata)[0]==np.size(meshdata):
        print('mesh file:', meshdata)
        volume=meshdata[volcol]
        xloc=meshdata[1]
        yloc=meshdata[2]
        depth=meshdata[volcol+1]
        SurfaceLayer=np.abs(yloc-meshdata[8])
        maxvols=1
        TotVol=volume
    else:
        if prninfo:
            print('1st mesh file row:', meshdata[0,:])
        volume=meshdata[:,volcol]
        depth=meshdata[:,volcol+1]
        maxvols=len(volume)
        xloc=meshdata[:,1]
        yloc=meshdata[:,2]
        SurfaceLayer=np.abs(yloc-meshdata[:,8])
        TotVol=0
        for k in range(maxvols):
            TotVol+=volume[k]
    if prninfo:
        print("TotVol:", TotVol, "\ndepth", depth, "\ndeltaY", SurfaceLayer)
    return maxvols,volume,xloc,yloc,TotVol,SurfaceLayer

def region_volume(List,Vox,volume,prnvox):
    #This volume is in units of cubic microns, multiply by 1e-15 to convert to Liters
    print("\nFOR region avg: j,regionvox,vol:", end=' ')
    region_volume=np.zeros(len(List))
    for j in range(len(List)):
        for k in Vox[j]:
            region_volume[j]+=volume[k] 
        if prnvox:
            print(j, List[j],Vox[j],region_volume[j])
        else:
            print("not printed")
    return region_volume

def subvol_list(structType,regionID,volnum,fake):
    yes=1
    no=0
    #all voxels should be read in now with labels
    #extract number of unique regions (e.g. dendrite, or sa1[0]), 
    #and create list of subvolumes which contribute to that region
    regionList=list()
    regionVox=list()
    regionCol=list()
    regionStructList=list()
    regionStructVox=list()
    regionStructCol=list()

    for i in range(len(volnum)):
        if (structType[i] != fake):
            unique=yes
            regionStruct=regionID[i]+structType[i]
            uniqueStruct=yes
            #print "vol,regionStruct",i,regionStruct,"num regions:", len(regionList)
            #first construct list of region specific voxels
            j=0
            while ((unique==yes) and (j<len(regionList))):
                if regionID[i]==regionList[j]:
                    regionVox[j].append(volnum[i])
                    regionCol[j].append(i)
                    unique=no
                #endif
                j+=1
            #endwhile
            if (unique==yes):
                regionList.append(regionID[i])
                regionVox.append([])
                regionCol.append([])
                regionVox[len(regionList)-1].append(volnum[i])
                regionCol[len(regionList)-1].append(i)
            #endif
            #second construct list of region and structure specific voxels.
            j=0
            while ((uniqueStruct==yes) and (j<len(regionStructList))):
                if (regionStruct==regionStructList[j]):
                    regionStructVox[j].append(volnum[i])
                    regionStructCol[j].append(i)
                    uniqueStruct=no
                #endif
                j+=1
            #endwhile
            if (uniqueStruct==yes):
                regionStructList.append(regionStruct)
                regionStructVox.append([])
                regionStructCol.append([])
                regionStructVox[len(regionStructList)-1].append(volnum[i])
                regionStructCol[len(regionStructList)-1].append(i)
            #endif
    #end for i
    return regionList,regionVox,regionCol,regionStructList,regionStructVox,regionStructCol

#Not yet debugged
def spatial_average(xloc,yloc,bins,regionID,structType,volnum):
    binmin=np.zeros(bins+1)
    binvoxels=[]
    bincolumns=[]
    xdim=max(xloc)-min(xloc)
    ydim=max(yloc)-min(yloc)

    #First, create the bins
    #xdim is larger:
    if (xdim >= ydim):
        bininc=xdim/bins
        minloc=min(xloc)
        maxloc=max(xloc)
        loc=xloc
        spaceheader='#time, x='
    #ydim is larger:
    else:
        bininc=ydim/bins
        minloc=min(yloc)
        maxloc=max(yloc)
        loc=yloc
        spaceheader='#time, y='
#now that assignmenst are made, determine bins the same for x or y direciton
    for j in range(bins):
        binmin[j]=minloc+j*bininc
        binvoxels.append([])
        bincolumns.append([])
    binmin[bins]=maxloc+bininc/bins
     
    #Now assign voxels to bins, either by x or y location
    print("binmin: ",binmin[:])
    print("regions=", len(regionID), "structures=", len(structType), "meshfile=", len(volume), "volnums=", len(volnum))
    for k in range(len(volnum)):
        if (regionID[k]==dend):
            j=0
            assigned=no
            while ((j < bins) and (assigned==no)):
                if ((loc[k]>=binmin[j]) and (loc[k] < binmin[j+1])):
                    binvoxels[j].append(volnum[k])
                    bincolumns[j].append(k)
                    assigned=yes
                j=j+1
            
    #now calculate volume for these spatial bins
    SpatialVol=np.zeros(bins)
    print("\nFOR spatial average: j, vol, binvoxels, bincolums:")
    for j in range(bins):
        for k in binvoxels[j]:
            SpatialVol[j]+=volume[k]
        if (prnvox==1):
            print(j, SpatialVol[j], binvoxels[j], bincolumns[j])
        else:
            print("not printed")

    for j in range(bins):
        spaceheader=spaceheader+' '+str(format(binmin[j],'.1f'))+'-'+str(format(binmin[j]+bininc,'.1f'))
    print('Space header:', spaceheader)
    return

def region_means(data,regionList,regionCol,regionVol,time,molecule):
    RegionMeans=np.zeros((len(time),len(regionList)))
    header=''       #Header for output file
    for itime in range(len(time)):
        for j in range(len(regionList)):
            for k in regionCol[j]:
                RegionMeans[itime,j]+=data[itime,k]

    for j in range(len(regionList)):
        RegionMeans[:,j] /= regionVol[j]*mol_per_nM_u3
        #print "head",header,"mol",molecule,"reg",regionList[j]
        header=header+molecule+regionList[j]+' '       #Header for output file
    return header,RegionMeans

def readdata(fname,maxvols,molecules,finalrow=0):
    arrays=len(molecules)
    lines=[]
    timelist=[]
    convertarray=0
    if finalrow:
        f=open(fname, 'r')
        head=f.readline()
        convertarray=1
        for i in range(finalrow):
            templine=f.readline()
            timelist.append(float(templine.split()[0]))
            lines.append([int(y) for y in templine.split()[1:]])
    else:
        try:
            alldata=np.loadtxt(fname,skiprows=1)
            time=alldata[:,0]/1000
            data=alldata[:,1:alldata.shape[1]]
        except ValueError:
            f=open(fname, 'r')
            head=f.readline()
            convertarray=1
            for templine in f:
                timelist.append(float(templine.split()[0]))
                lines.append([int(y) for y in templine.split()[1:]])
    if convertarray:
        time=np.array(timelist[0:-1])/1000
        data=np.array(lines[0:-1])
    # eliminate the time column from the data, so that e.g., column 0 = voxel 0
    #
    #reshape the data to create a separate dimension for each molecule
    rows=np.shape(data)[0]
    if maxvols*arrays == np.shape(data)[1]:
            molecule_array=np.reshape(data, (rows,arrays,maxvols))
    else:
            print("UH OH! voxels:", maxvols, "molecules:", len(molecules), "columns:", np.shape(data)[1],"arrays",arrays)
    return time,molecule_array,rows

 
