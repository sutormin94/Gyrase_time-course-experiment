###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#The script takes raw GCSs data, returns only trusted GCSs, 
#computes GCSs shared between different conditions, 
#draws Venn diagrams of the sets overlappings, 
#writes GCSs sets.
###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn3, venn3_circles
from matplotlib import cm as cm
import collections
from collections import OrderedDict
import numpy as np
import pandas as pd
from pandas import DataFrame
from scipy.cluster import hierarchy as hc

#######
#Variables to be defined.
#######

print('Variables to be defined:')

#Path to the working directory
pwd="C:\Sutor\science\DNA-gyrase\Results\E_coli_synch_time-course_Topo-Seq\Seq_results\GCSs_analysis\GCSs_calling_0_01\\"

#Input data
path_to_tp_replicas={'-3_min' : {'-3_min_R1': pwd + "R1\\NN\\-3_min\\-3_min_raw_GCSs_called.txt",
                                      '-3_min_R2': pwd + "R2\\NN\\-3_min\\-3_min_raw_GCSs_called.txt"},
                          '0_min' : {'0_min_R1': pwd + "R1\\NN\\0_min\\0_min_raw_GCSs_called.txt",
                                     '0_min_R2': pwd + "R2\\NN\\0_min\\0_min_raw_GCSs_called.txt"},
                          '5_min' : {'5_min_R1': pwd + "R1\\NN\\5_min\\5_min_raw_GCSs_called.txt",
                                     '5_min_R2': pwd + "R2\\NN\\5_min\\5_min_raw_GCSs_called.txt"},
                          '10_min' : {'10_min_R1': pwd + "R1\\NN\\10_min\\10_min_raw_GCSs_called.txt",
                                      '10_min_R2': pwd + "R2\\NN\\10_min\\10_min_raw_GCSs_called.txt"},
                          '15_min' : {'15_min_R1': pwd + "R1\\NN\\15_min\\15_min_raw_GCSs_called.txt",
                                      '15_min_R2': pwd + "R2\\NN\\15_min\\15_min_raw_GCSs_called.txt"},
                          '20_min' : {'20_min_R1': pwd + "R1\\NN\\20_min\\20_min_raw_GCSs_called.txt",
                                      '20_min_R2': pwd + "R2\\NN\\20_min\\20_min_raw_GCSs_called.txt"},
                          '25_min' : {'25_min_R1': pwd + "R1\\NN\\25_min\\25_min_raw_GCSs_called.txt",
                                      '25_min_R2': pwd + "R2\\NN\\25_min\\25_min_raw_GCSs_called.txt"},
                          '30_min' : {'30_min_R1': pwd + "R1\\NN\\30_min\\30_min_raw_GCSs_called.txt",
                                      '30_min_R2': pwd + "R2\\NN\\30_min\\30_min_raw_GCSs_called.txt"}}

#Input data in one dict (for one output table contains all replices and raw GCSs)
path_to_replicas={'-3_min_R1': pwd + "R1\\NN\\-3_min\\-3_min_raw_GCSs_called.txt",
                  '-3_min_R2': pwd + "R2\\NN\\-3_min\\-3_min_raw_GCSs_called.txt",
                  '0_min_R1': pwd + "R1\\NN\\0_min\\0_min_raw_GCSs_called.txt",
                  '0_min_R2': pwd + "R2\\NN\\0_min\\0_min_raw_GCSs_called.txt",
                  '5_min_R1': pwd + "R1\\NN\\5_min\\5_min_raw_GCSs_called.txt",
                  '5_min_R2': pwd + "R2\\NN\\5_min\\5_min_raw_GCSs_called.txt",
                  '10_min_R1': pwd + "R1\\NN\\10_min\\10_min_raw_GCSs_called.txt",
                  '10_min_R2': pwd + "R2\\NN\\10_min\\10_min_raw_GCSs_called.txt",
                  '15_min_R1': pwd + "R1\\NN\\15_min\\15_min_raw_GCSs_called.txt",
                  '15_min_R2': pwd + "R2\\NN\\15_min\\15_min_raw_GCSs_called.txt",
                  '20_min_R1': pwd + "R1\\NN\\20_min\\20_min_raw_GCSs_called.txt",
                  '20_min_R2': pwd + "R2\\NN\\20_min\\20_min_raw_GCSs_called.txt",
                  '25_min_R1': pwd + "R1\\NN\\25_min\\25_min_raw_GCSs_called.txt",
                  '25_min_R2': pwd + "R2\\NN\\25_min\\25_min_raw_GCSs_called.txt",
                  '30_min_R1': pwd + "R1\\NN\\30_min\\30_min_raw_GCSs_called.txt",
                  '30_min_R2': pwd + "R2\\NN\\30_min\\30_min_raw_GCSs_called.txt"}

#Configuration of the output for the GCSs data in replicas.
Replicas_path_out=pwd + "GCSs_sets\\"
if not os.path.exists(Replicas_path_out):
    os.makedirs(Replicas_path_out)

All_conditions_name="All_time-points_GCSs"
All_trusted_conditions_name="All_time-points_trusted_GCSs"

#Outpath for Venn diagrams.
plot_outpath=Replicas_path_out + "Replicas_venn\\"
if not os.path.exists(plot_outpath):
    os.makedirs(plot_outpath)
    
#Outpath for emerging GCSs distributions.
seq_plot_outpath=Replicas_path_out + "Seq_comp_GCSs\\"
if not os.path.exists(seq_plot_outpath):
    os.makedirs(seq_plot_outpath)



#########
##Pars GCSs files, combine into one dataframe and write table contains all GCSs of all replicas of all time-points.
#########
#######
#Parsing raw GCSs coordinates, returns dictionary - GCSs_coordinate:N3E.
#######

def read_GCSs_file(GCSs_file_path):
    GCSs_dict={}
    GCSs_in=open(GCSs_file_path, 'r')
    for line in GCSs_in:
        line=line.rstrip().split('\t')
        if line[0] not in ['GCSs_coordinate']:
            GCSs_dict[int(line[0])]=float(line[1])
    GCSs_in.close()
    return GCSs_dict

#######
#Combines replicates into one GCSs table.
#######

def combine_replicates(replicas_dict, path_out, name):
    #Merges a range of replicates
    GCSs_replicas_dict={}
    names_ar=[]
    for key, value in replicas_dict.items(): #Iterates replicas
        names_ar.append(key)
        #Read file with raw GCSs
        Raw_GCSs_dict=read_GCSs_file(value)
        for k, v in Raw_GCSs_dict.items(): #Iterates raw GCSs
            #Table filling process initiation
            if len(names_ar)==1:
                GCSs_replicas_dict[k]=[v]
            #Table filling process continuing (the table already contains at least one GCSs set)
            else:
                #If GCSs is already in the table
                if k in GCSs_replicas_dict:
                    GCSs_replicas_dict[k].append(v)
                #If this is the first occurrence of the element in a NON empty table.
                else:
                    add_el=[]
                    for j in range(len(names_ar)-1):
                        add_el.append(0)
                    add_el.append(v)
                    GCSs_replicas_dict[k]=add_el
        #If table body line contains less elements than header does, hence add zero.
        for k, v in GCSs_replicas_dict.items():
            if len(v)<len(names_ar):
                GCSs_replicas_dict[k].append(0)
    #Sorting the list of dictionary keys.
    GCSs_replicas_dict_sorted=collections.OrderedDict(sorted(GCSs_replicas_dict.items()))
    #Writes merged GCSs data
    #fileout=open(path_out + name + '_GCSs_replicates.txt', 'w')
    #Header
    #fileout.write('GCSs_coordinate\t')
    #for i in names_ar:
    #    fileout.write(str(i) + '_N3E\t')
    #fileout.write('\n')
    #Body of the table
    #for k, v in GCSs_replicas_dict_sorted.items():
    #    fileout.write(str(k) + '\t')
    #    for i in GCSs_replicas_dict_sorted[k]:
    #        fileout.write(str(i) + '\t')
    #    fileout.write('\n')
    #fileout.close()
    return GCSs_replicas_dict, names_ar
 
#Prepares GCSs table for all time-points and replicas.
All_replicas_dict_of_ars, names_ar=combine_replicates(path_to_replicas, Replicas_path_out, All_conditions_name)      


#########
##Returns only trusted GCSs - observed in both biological replicas.
##Data organization: 1. coordinate of GCSs, 2.-3. N3E values for biological replicates 1-2
#########

#Check if GCS is trusted.
def trusted(ar):
    av_height=0
    ind=0
    for i in range(len(ar)):
        if ar[i]>0:
            ind=ind+1
            av_height=av_height+ar[i]
    if ind>1:
        return av_height/ind
    else:
        return "No signal"

#Iterate over GSCs check if GCS is trusted.
def trusted_GCSs_calling(GCSs_dictionary):
    ar=[]
    for k, v in GCSs_dictionary.items():
        if trusted(v)!="No signal":
            ar.append([k, trusted(v)])
    return ar

#Read files, assemble into one dataframe, check if GCSs are trusted.
def replicas_comb_trust_wrapper(dict_of_replicas_dict, path_out):
    dict_of_trusted_ars={}
    for name, replicas_dict in dict_of_replicas_dict.items():
        print('Now working with: ' + str(name))
        cur_GCSs_dict=combine_replicates(replicas_dict, path_out, name)[0]
        cur_GCSs_trusted=trusted_GCSs_calling(cur_GCSs_dict)
        print('Number of trusted GCSs for ' + str(name) + ' : ' + str(len(cur_GCSs_trusted)))
        print('Fraction of trusted GCSs for ' + str(name) + ' : ' + str(float(len(cur_GCSs_trusted))/len(cur_GCSs_dict)))
        dict_of_trusted_ars[name]=cur_GCSs_trusted
    return dict_of_trusted_ars

all_tp_trusted=replicas_comb_trust_wrapper(path_to_tp_replicas, Replicas_path_out)

#Write trusted GSCs data.
def write_GCSs_file(dict_of_ars, path_out):
    for k, v in dict_of_ars.items(): #Iterates lists to be written
        v.sort(key=lambda tup: tup[0])  #Sorting lists by the zero elements of the sublists they consist of 
        fileout=open(path_out + k + '_GCSs_trusted.txt', 'w')
        fileout.write('GCSs_coordinate\tN3E\n')
        for i in range(len(v)):
            fileout.write(str(v[i][0]) + '\t' + str(v[i][1]) + '\n')
        fileout.close()
    return

#write_GCSs_file(all_tp_trusted, Replicas_path_out)

#########
##Assemble ars of trusted GCSs into one dataframe.
#########
#######
#Convert dictionary of arrays into dictionary of dictionaries.
#######

#print(all_tp_trusted)

def ar_to_dict(dict_of_ars):
    dict_of_dicts={}
    for key, ar in dict_of_ars.items():
        dict_of_GCSs={}
        for gcs in ar:
            dict_of_GCSs[gcs[0]]=gcs[1]
        dict_of_dicts[key]=dict_of_GCSs
    return dict_of_dicts

all_tp_trusted_dicts=ar_to_dict(all_tp_trusted)
#print(all_tp_trusted_dicts)

#######
#Combine dictionaries of trusted GCSs into one dataframe and write it down.
#######

def combine_trusted(dict_of_dicts, path_out, name):
    #Merges a range of replicates
    GCSs_replicas_dict={}
    names_ar=[]
    for key, value in dict_of_dicts.items(): #Iterates dicts of trusted GCSs
        names_ar.append(key)
        Trusted_GCSs_dict=value
        for k, v in Trusted_GCSs_dict.items(): #Iterates trusted GCSs
            #Table filling process initiation
            if len(names_ar)==1:
                GCSs_replicas_dict[k]=[v]
            #Table filling process continuing (the table already contains at least one GCSs set)
            else:
                #If GCSs is already in the table
                if k in GCSs_replicas_dict:
                    GCSs_replicas_dict[k].append(v)
                #If this is the first occurrence of the element in a NON empty table.
                else:
                    add_el=[]
                    for j in range(len(names_ar)-1):
                        add_el.append(0)
                    add_el.append(v)
                    GCSs_replicas_dict[k]=add_el
        #If table body line contains less elements than header does, hence add zero.
        for k, v in GCSs_replicas_dict.items():
            if len(v)<len(names_ar):
                GCSs_replicas_dict[k].append(0)
    #Sorting the list of dictionary keys.
    GCSs_replicas_dict_sorted=collections.OrderedDict(sorted(GCSs_replicas_dict.items()))
    #Writes merged trusted GCSs data
    #fileout=open(path_out + name + '_trusted_GCSs.txt', 'w')
    #Header
    #fileout.write('GCSs_coordinate\t')
    #for i in names_ar:
    #    fileout.write(str(i) + '_N3E\t')
    #fileout.write('\n')
    #Body of the table
    #for k, v in GCSs_replicas_dict_sorted.items():
    #    fileout.write(str(k) + '\t')
    #    for i in GCSs_replicas_dict_sorted[k]:
    #        fileout.write(str(i) + '\t')
    #    fileout.write('\n')
    #fileout.close()    
    return

#combine_trusted(all_tp_trusted_dicts, Replicas_path_out, All_trusted_conditions_name)


#########
##Make Venn diagrams for biological replicas.
#########
#######
#GCSs shared between two arrays (biological replicas).
#######

def pairs_construction(ar1, ar2):
    double=[]
    for i in range(len(ar1)):
        for j in range(len(ar2)):
            if ar1[i][0]==ar2[j][0]:
                double.append([ar1[i][0], ar1[i][1], ar2[j][1]]) #GCSs coordinate, N3E_1, N3E_2 
    return double

#######
#Parses replicas, overlaps lists of GCSs, output data for Venn diagram construction.
#######

def replicates_parsing_to_list_and_overlapping(replicas_dict, name):
    #Parsing
    GCSs_dict={}
    for k, v in replicas_dict.items(): #Iterate replicas.
        GCSs_dict[k]=[]
        for c, h in read_GCSs_file(v).items(): #Iterate GCSs.
            GCSs_dict[k].append([c, h])
    #Overlapping
    one_two=pairs_construction(GCSs_dict[name+'_R1'], GCSs_dict[name+'_R2'])
    #Venn input description (for 2 sets): one, two, one_two
    venn_input=[len(GCSs_dict[name+'_R1'])-len(one_two), 
                len(GCSs_dict[name+'_R2'])-len(one_two), 
                len(one_two)]
    return venn_input

def replicas_venn_wrapper(dict_of_replicas_dict, plot_outpath):
    for name, replicas_dict in dict_of_replicas_dict.items():
        print('Now working with: ' + str(name))
        venn2data=replicates_parsing_to_list_and_overlapping(replicas_dict, name)
        venn2(subsets = (venn2data), set_labels = (name+"_R1", name+"_R2"))
        #plt.savefig(plot_outpath+name+'_R1_R2_venn.png', dpi=320)
        plt.close()    
    return

#replicas_venn_wrapper(path_to_tp_replicas, plot_outpath)


#########
##Compute correlation matrix and draw heatmaps.
#########
#######
#Construct correlation matrix between replicas.
#######

#Plot diagonal correlation matrix.
def correlation_matrix(df, cor_method, title, outpath):
    fig=plt.figure(figsize=(8,8), dpi=100)
    ax1=fig.add_subplot(111)
    cmap=cm.get_cmap('rainbow', 30)
    cax=ax1.imshow(df.corr(method=cor_method), interpolation="nearest", cmap=cmap)
    ax1.grid(True, which='minor', linestyle="--", linewidth=0.5, color="black")
    plt.title(title)
    labels=list(df)
    ax1.set_xticks(np.arange(len(labels)))
    ax1.set_yticks(np.arange(len(labels)))    
    ax1.set_xticklabels(labels, fontsize=12, rotation=90)
    ax1.set_yticklabels(labels, fontsize=12)
    #Add colorbar, make sure to specify tick locations to match desired ticklabels
    fig.colorbar(cax, ticks=[-0.10, -0.05, 0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00])
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()
    return

#Plot dendrogram on the basis of correlation matrix.
def correlation_dendrogram(df, cor_method, title, outpath):
    corr_inv=1-df.corr(method=cor_method) #compute correlation and inverse to distance
    corr_inv_condensed=hc.distance.squareform(corr_inv) #convert to condensed
    z=hc.linkage(corr_inv_condensed, method='average')
    dendrogram=hc.dendrogram(z, labels=corr_inv.columns, leaf_rotation=90)
    plt.title(title)
    plt.tight_layout()
    plt.savefig(outpath, dpi=400, figsize=(8, 8))
    plt.show()  
    plt.close()
    return

#Plotting the heatmap for all replicas.
All_replicas_df=pd.read_csv(Replicas_path_out+All_conditions_name+'_GCSs_replicates.txt', sep='\t', header=(0)).drop(labels=['GCSs_coordinate', 'Unnamed: 17'], axis='columns')
#correlation_matrix(All_replicas_df, 'spearman', "Correlation between all samples of\ntime-course gyrase Topo-Seq experiment", Replicas_path_out+"Samples_correlation\\All_tc_samples_cm_spearman.png")
#correlation_matrix(All_replicas_df, 'pearson', "Correlation between all samples of\ntime-course gyrase Topo-Seq experiment", Replicas_path_out+"Samples_correlation\\All_tc_samples_cm_pearson.png")
#correlation_matrix(All_replicas_df, 'kendall', "Correlation between all samples of\ntime-course gyrase Topo-Seq experiment", Replicas_path_out+"Samples_correlation\\All_tc_samples_cm_kendall.png")
#correlation_dendrogram(All_replicas_df, 'spearman', "Dendrogram of distances between all samples of time-course\ngyrase Topo-Seq experiment", Replicas_path_out+"Samples_correlation\\All_tc_samples_cm_spearman_dendrogram.png")
#correlation_dendrogram(All_replicas_df, 'pearson', "Dendrogram of distances between all samples of time-course\ngyrase Topo-Seq experiment", Replicas_path_out+"Samples_correlation\\All_tc_samples_cm_pearson_dendrogram.png")
#correlation_dendrogram(All_replicas_df, 'kendall', "Dendrogram of distances between all samples of time-course\ngyrase Topo-Seq experiment", Replicas_path_out+"Samples_correlation\\All_tc_samples_cm_kendall_dendrogram.png")


#Plotting the heatmap for all trusted GCSs lists.
All_trusted_GSCs_df=pd.read_csv(Replicas_path_out+All_trusted_conditions_name+'_trusted_GCSs.txt', sep='\t', header=(0)).drop(labels=['GCSs_coordinate', 'Unnamed: 9'], axis='columns')
#correlation_matrix(All_trusted_GSCs_df, 'spearman', "Correlation between all samples of time-course\ngyrase Topo-Seq experiment (trusted GCSs)", Replicas_path_out+"Samples_correlation\\All_tc_trusted_samples_cm_spearman.png")
#correlation_matrix(All_trusted_GSCs_df, 'pearson', "Correlation between all samples of time-course\ngyrase Topo-Seq experiment (trusted GCSs)", Replicas_path_out+"Samples_correlation\\All_tc_trusted_samples_cm_pearson.png")
#correlation_matrix(All_trusted_GSCs_df, 'kendall', "Correlation between all samples of time-course\ngyrase Topo-Seq experiment (trusted GCSs)", Replicas_path_out+"Samples_correlation\\All_tc_trusted_samples_cm_kendall.png")
#correlation_dendrogram(All_trusted_GSCs_df, 'spearman', "Dendrogram of distances between all samples of time-course\ngyrase Topo-Seq experiment (trusted GCSs)", Replicas_path_out+"Samples_correlation\\All_tc_trusted_samples_cm_spearman_dendrogram.png")
#correlation_dendrogram(All_trusted_GSCs_df, 'pearson', "Dendrogram of distances between all samples of time-course\ngyrase Topo-Seq experiment (trusted GCSs)", Replicas_path_out+"Samples_correlation\\All_tc_trusted_samples_cm_pearson_dendrogram.png")
#correlation_dendrogram(All_trusted_GSCs_df, 'kendall', "Dendrogram of distances between all samples of time-course\ngyrase Topo-Seq experiment (trusted GCSs)", Replicas_path_out+"Samples_correlation\\All_tc_trusted_samples_cm_kendall_dendrogram.png")

 
#########
##Sequentially compare trusted GCSs of time-points.
######### 

#Compare two dects of GCSs.
def compare_two_dicts(GCSs_dict_1, GSCs_dict_2):
    common_GCSs={}
    new_GCSs={}
    for key, GCS in GSCs_dict_2.items():
        if key in GCSs_dict_1:
            common_GCSs[key]=GCS
        else:
            new_GCSs[key]=GCS
    return common_GCSs, new_GCSs

#Convert one dict to two ars.
def dict_to_ars(dictionary):
    ar_keys=[]
    ar_values=[]
    for key, value in dictionary.items():
        ar_keys.append(key)
        ar_values.append(value)
    return ar_keys, ar_values

#Sequentially compare time-points.
def compare_trusted_GSCs_dicts(dict_of_dicts, path_out):
    names_ar=[]
    for key, gcss_dict in dict_of_dicts.items():
        names_ar.append(key)
    print(names_ar)
    
    New_GCSs_ar_dict={}
    for i in range(len(names_ar)-1):
        common_GCSs, new_GCSs=compare_two_dicts(dict_of_dicts[names_ar[i]], dict_of_dicts[names_ar[i+1]])
        coordinates_ar, N3E_ar=dict_to_ars(new_GCSs)
        New_GCSs_ar_dict[names_ar[i+1]+'_'+names_ar[i]]=coordinates_ar
        print(names_ar[i+1], names_ar[i], len(new_GCSs))
    
    #Parameters
    genome_len=4647454
    bins=np.linspace(0,genome_len,1001)
    ticks1=[0, 500000, 1000000, 1500000, 2000000, 2500000, 3000000, 3500000, 4000000, 4500000]
    xticknames1=['', '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500']
    colors=['#7FCE79', '#BAE85C', '#ff878b', '#8991ff', '#ac5eff', '#50b3ff', '#ffd75e']
    plot_names=['plot1', 'plot2', 'plot3', 'plot4', 'plot5', 'plot6', 'plot7']
    Y_labels=['0min\\\n-3min', '5min\\\n0min', '10min\\\n5min', '15min\\\n10min', '20min\\\n15min', '25min\\\n20min', '30min\\\n25min']
    yter=1592477
    yori=3711828
    #GCSs data plotting.
    fig, plot_names=plt.subplots(7,1,figsize=(11,15), dpi=100)
    i=0
    Histo_comp_dict={} #Will contain computed histogramm data (bins and values)
    for key, value in New_GCSs_ar_dict.items():
        plot_names[i].set_xlim(0, genome_len)
        plot_names[i].set_xticks(ticks1, minor=False)
        plot_names[i].set_xticks([yter, yori], minor=True)
        plot_names[i].set_xticklabels(xticknames1)
        plt.setp(plot_names[i].set_xticklabels(xticknames1), rotation=0, fontsize=14)
        plot_names[i].locator_params(axis='y', nbins=6)
        plot_names[i].tick_params(axis='x', which='major', labelsize=19)
        Histo_comp_dict[key]=plot_names[i].hist(value, bins, facecolor=colors[i], alpha=0.7, linewidth=0.1, edgecolor='black') #Plot histo and save computed histogramm data (bins and values)
        plot_names[i].tick_params(axis='y', which='major', pad=7, labelsize=15)
        plot_names[i].set_ylabel(Y_labels[i], size=22, labelpad=8, rotation=90)
        i+=1
    plt.tight_layout()
    fig.savefig(path_out+"New_trusted_GCSs_time-course_1000_bins.png", figsize=(11,15), dpi=400)
    plt.close()
    return

compare_trusted_GSCs_dicts(all_tp_trusted_dicts, seq_plot_outpath)
 
print('Script ended its work succesfully!') 
 