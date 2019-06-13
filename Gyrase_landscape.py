###############################################
##Dmitry Sutormin, 2019##
##Gyrase Topo-Seq analysis##

####
#The only purpose - to compute by-position average of a set of wig files.
####

###############################################

#######
#Packages to be imported.
#######

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import cm as cm
from mpl_toolkits.mplot3d import Axes3D
import seaborn as sns
import mayavi
from mayavi import mlab
from mayavi.mlab import *
from tvtk.api import tvtk

#Dictionary of time points 
#'Replica name' : 'Path to wig file'
Dict_of_replicas={1 : "F:\Gyrase_landscape\Data\FE\Gyrase_tc_Topo-Seq_0min_average_FE_Early_Stat.wig",
                  2 : "F:\Gyrase_landscape\Data\FE\Gyrase_tc_Topo-Seq_5min_average_FE_Early_Stat.wig",
                  3 : "F:\Gyrase_landscape\Data\FE\Gyrase_tc_Topo-Seq_10min_average_FE_Early_Stat.wig",
                  4 : "F:\Gyrase_landscape\Data\FE\Gyrase_tc_Topo-Seq_15min_average_FE_Early_Stat.wig",
                  5 : "F:\Gyrase_landscape\Data\FE\Gyrase_tc_Topo-Seq_20min_average_FE_Early_Stat.wig",
                  6 : "F:\Gyrase_landscape\Data\FE\Gyrase_tc_Topo-Seq_25min_average_FE_Early_Stat.wig",
                  7 : "F:\Gyrase_landscape\Data\FE\Gyrase_tc_Topo-Seq_30min_average_FE_Early_Stat.wig",
                  8 : "F:\Gyrase_landscape\Data\FE\Gyrase_Cfx_10mkM_average_FE_Mid_exp.wig",
                  9 : "F:\Gyrase_landscape\Data\FE\Gyrase_tc_Topo-Seq_-3min_average_FE_Early_Stat.wig",
                  10 : "F:\Gyrase_landscape\Data\FE\Gyrase_Topo-Seq_Stationary_average_FE_Stat.wig",
                  }
#Deletions info (regions of a genome to mask).
Deletions=np.array([[274500, 372148], [793800, 807500], [1199000, 1214000]])
#Path to the csv with wig files merged to form a dataframe.
Pathin="F:\Gyrase_landscape\Data\Gyrase_landscape_binned_1000_eq_scaled.csv"
#Path for output correlation matrix.
Outpath="F:\Gyrase_landscape\\"


#######
#Parses WIG file.
#######

def wig_parsing(wigfile):
    print('Now is processing: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    NE_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            NE_values.append(float(line[0]))
    wigin.close()
    return NE_values


#########
##Compute correlation matrix and draw heatmaps.
#########

#Plot diagonal correlation matrix.
def correlation_matrix(df, cor_method, title, outpath):
    fig=plt.figure(figsize=(8,8), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    cax=ax1.imshow(df.corr(method=cor_method), interpolation="nearest", cmap=cmap, norm=None, vmin=-1, vmax=1)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    #Add colorbar, make sure to specify tick locations to match desired ticklabels.
    #Full scale:[-1.00, -0.95, -0.90, -0.85, -0.80, -0.75, -0.70, -0.65, -0.60, -0.55, -0.50, -0.45, -0.40, -0.35, -0.30, -0.25, -0.20, -0.15, -0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    fig.colorbar(cax, ticks=[-1.00, -0.90, -0.80, -0.70, -0.60, -0.50, -0.40, -0.30, -0.20, -0.10, 0.00, 0.10, 0.20, 0.30, 0.40, 0.50, 0.60, 0.70, 0.80, 0.90, 1.00], shrink=0.7)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()
    plt.close()
    return


#########
##Read input wig files and convert them.
#########

def read_wig_convert(Dict_of_replicas, binning_option):
    #Contains data of all replicas in separate arrays.
    dict_of_replicas={}
    dict_of_scalings={}
    min_scaling=10000000
    for replica_name, replica_path in Dict_of_replicas.items():
        replic_data=wig_parsing(replica_path)
        replic_data_ar=np.array(replic_data)
        data_mean=np.mean(replic_data_ar)
        dict_of_scalings[replica_name]=data_mean
        if min_scaling>data_mean:
            min_scaling=data_mean
        if binning_option>0:
            replic_data_ar_binned=replic_data_ar[:(replic_data_ar.size // binning_option) * binning_option].reshape(-1, binning_option).mean(axis=1)  
            dict_of_replicas[replica_name]=replic_data_ar_binned
        elif binning_option==0:
            dict_of_replicas[replica_name]=replic_data_ar
        print(replica_name, len(replic_data_ar), min_scaling)
    
    print(dict_of_scalings)
    #Scale datasets.
    dict_of_replicas_scaled={}
    for replica_name, replica_data in dict_of_replicas.items():
        replica_data_scaled=(replica_data/dict_of_scalings[replica_name])*1
        dict_of_replicas_scaled[replica_name]=replica_data_scaled
        print(f'{replica_name} mean, sum: {np.mean(replica_data_scaled)}, {np.sum(replica_data_scaled)}')
        
    Gyrase_dataframe=pd.DataFrame(dict_of_replicas_scaled)
    Gyrase_dataframe=Gyrase_dataframe[[1, 2, 3, 4, 5, 6, 7, 8, 9, 10]]
    return Gyrase_dataframe


##Implement smoothing!

#Gyrase_dataframe=read_wig_convert(Dict_of_replicas, 1000) 
#Gyrase_dataframe.to_csv(Pathin, sep='\t', index=False)

Gyrase_dataframe=pd.read_csv(Pathin, sep='\t', header=0)
print('Dataframe was read')

#correlation_matrix(Gyrase_dataframe, 'pearson', 'Correlation of Gyrase Topo-Seq time-points', Outpath+"Figures\Gyrase_Topo-Seq_time_points_correlation_matrix_binned_1000_eq_scaled.png")
#print('Correlation matrix was constructed!')


#########
##Subsamle dataframe and scale y and z axis.
#########

def scale_subsample_mask(Gyrase_dataframe, top, scale_genome, scale_points, scale_signal, Deletions, binning_option):
    #Subsample data
    Gyrase_dataframe=Gyrase_dataframe.head(top)
    #print(Gyrase_dataframe)
    Gyrase_dataframe_m=Gyrase_dataframe.as_matrix()
    Gyrase_dataframe_m=Gyrase_dataframe_m*scale_signal
    #print(Gyrase_dataframe_m)
    
    #Prepare data with special format for mayavi surface.
    df_dims=Gyrase_dataframe.shape
    rows_num=df_dims[0]
    cols_num=df_dims[1]

    x_list=[]
    y_list=[]
    for i in range(rows_num):
        x=[]
        y=[]
        for j in range(cols_num):
            x.append(i)
            y.append(j)
        x_list.append(x)
        y_list.append(y)
    
    #Scale data.
    x_array=np.array(x_list, np.float64)
    x_array=x_array/scale_genome
    y_array=np.array(y_list, np.float64)
    y_array=y_array/scale_points
    
    #Mask data.
    #print(Deletions)
    Deletions=Deletions//binning_option
    #print(Deletions)
    Deletions=Deletions/scale_genome
    #print(Deletions)
    x_array_mask_1=((Deletions[0][1]>=x_array) & (x_array>=Deletions[0][0]))
    x_array_mask_2=((Deletions[1][1]>=x_array) & (x_array>=Deletions[1][0]))
    x_array_mask_3=((Deletions[2][1]>=x_array) & (x_array>=Deletions[2][0]))
    x_array_mask_12=np.logical_or(x_array_mask_1, x_array_mask_2)
    x_array_mask_123=np.logical_or(x_array_mask_12, x_array_mask_3)
    return x_array, y_array, Gyrase_dataframe_m, x_array_mask_123

x_array, y_array, z_array, masking_array=scale_subsample_mask(Gyrase_dataframe, len(Gyrase_dataframe), 100, 0.2, 1, Deletions, 1000)


m=mayavi.mlab.surf(x_array, y_array, z_array, extent=[np.min(x_array), np.max(x_array), np.min(y_array), np.max(y_array), z_array[z_array>0].min(), 25], mask=masking_array)
mlab.axes(xlabel='Positions, nt', ylabel='Time points', zlabel='Fold enrichment, scaled')

#Make a grid.
xx=np.arange(np.min(x_array)-2, np.max(x_array), 1)
yy=np.arange(np.min(y_array)-1, np.max(y_array), 1)
zz=np.arange(z_array[z_array>0].min(), 25, 1)
xy=xz=np.zeros_like(xx) 
yx=yz=np.zeros_like(yy)
zx=zy=np.zeros_like(zz) 
#Decorate y ticks.
ylabels=['0 min', '5 min', '10 min', '15 min', '20 min', '25 min', '30 min', 'Mid exp', 'Tran stat', 'Stat']
for i in range(10):
    lensoffset=((np.max(y_array)-np.min(y_array))/9)*i
    mlab.plot3d(zx+np.max(x_array),zy+lensoffset,zz,line_width=1,tube_radius=0.1)
    mlab.plot3d(xx,xy+lensoffset,xz,line_width=1,tube_radius=0.1)
    mlab.text3d(-8, lensoffset, 0, ylabels[i], color=(1,1,1), opacity=1)   
#Decorate x ticks.
xlabels=['rRNA D', 'rRNA C', 'rRNA G', 'rRNA E', 'rRNA B', 'rRNA H', 'rRNA A', 'Mu SGS']
xticks=[34.24400, 36.90771, 27.25847, 42.12776, 34.66047, 2.23771, 35.97167, 6.57709]
for i in range(len(xticks)):
    lensoffset=xticks[i]
    mlab.plot3d(zx+lensoffset,zy+np.max(y_array),zz,line_width=1,tube_radius=0.1)
    mlab.plot3d(yx+lensoffset,yy,yz,line_width=1,tube_radius=0.1)
    mlab.text3d(lensoffset, -1, 0, xlabels[i], color=(1,1,1), opacity=1)   
mayavi.mlab.show()


#mlab.plot3d(yx,yy+lensoffset,yz,line_width=1,tube_radius=1)
#mlab.plot3d(zx+np.max(x_array),zy+lensoffset,zz,line_width=1,tube_radius=0.1)
#mlab.plot3d(xx,xy+lensoffset,xz,line_width=1,tube_radius=0.1)


#mayavi.mlab.test_surf()
#mlab.axes(xlabel='y', ylabel='x', zlabel='z')
#mayavi.mlab.show()


def test_quiver3d():
    x, y, z = np.mgrid[-2:3, -2:3, -2:3]
    r = np.sqrt(x ** 2 + y ** 2 + z ** 4)
    u = y * np.sin(r) / (r + 0.001)
    v = -x * np.sin(r) / (r + 0.001)
    w = np.zeros_like(z)
    obj = quiver3d(x, y, z, u, v, w, line_width=3, scale_factor=1)
    return obj

#test_quiver3d()
#mayavi.mlab.show()


"""
# Transform it to a long format
df=Gyrase_dataframe.unstack().reset_index()
df.columns=["X","Y","Z"]
df['X'] = df['X'].astype(int)
df['Y'] = df['Y'].astype(float)
dfx=df['X']
dfy=df['Y']
dfz=df['Z']
#dfxt=dfx.transpose()
#dfyt=dfy.transpose()
#dfzt=dfz.transpose()
dfxt=dfx.as_matrix()
dfyt=dfy.as_matrix()
dfyt+=1
dfyt=dfyt/1
dfzt=dfz.as_matrix()
dfzt=dfzt*10

print(dfxt)
print(dfyt)
print(dfzt)

#print(Gyrase_dataframe.as_matrix())
"""


""" 
# And transform the old column name in something numeric
df['X']=pd.Categorical(df['X'])
df['X']=df['X'].cat.codes
 
# Make the plot
fig = plt.figure()
ax = fig.gca(projection='3d')
ax.plot_trisurf(df['Y'], df['X'], df['Z'], cmap=plt.cm.viridis, linewidth=0.2)
plt.show()
 
# to Add a color bar which maps values to colors.
surf=ax.plot_trisurf(df['Y'], df['X'], df['Z'], cmap=plt.cm.viridis, linewidth=0.2)
fig.colorbar( surf, shrink=0.5, aspect=5)
plt.show()
 
# Rotate it
ax.view_init(30, 45)
plt.show()
 
# Other palette
ax.plot_trisurf(df['Y'], df['X'], df['Z'], cmap=plt.cm.jet, linewidth=0.01)
plt.show()

# vtk DataFile Version 2.0
Dataset test10k
ASCII
DATASET STRUCTURED_GRID
DIMENSIONS 10 10000 1
POINTS 100000 float
"""