###############################################
##Dmitry Sutormin, 2019##
##Topo-Seq analysis##

#The script takes results of scanning procedure (WIG file), returns scores for
#GCSs (writes TAB files contain coordinate\tN3E\tScore info), 
#computes Pearson correlation between N3E and score, plots (Score, N3E) scatter plots.

###############################################

#######
#Packages to be imported.
#######

import os
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.offsetbox import AnchoredText
import numpy as np
import scipy
from scipy.stats import pearsonr
import colorsys

#######
#Variables to be defined.
#######

#Path to the input WIG score file.
path_to_res_score_files="C:\Sutor\science\DNA-gyrase\Results\GCSs_sets_and_motifs\Score_tracks\E_coli_w3110_G_Mu_score.wig"
#Path to the working directory.
pwd="F:\E_coli_Topo-Seqs\Gyrase_Topo-Seq_data_coverage_depth\GSCs_calling\\"

#Input data - GCSs, TAB.
path_to_GCSs_files={'EP_Cfx': pwd + "Trusted_GCSs_calling_0.01\EP_Cfx_10mkM_trusted_GCSs.txt",
                    'EP_RifCfx': pwd + "Trusted_GCSs_calling_0.01\EP_RifCfx_trusted_GCSs.txt",
                    'EP_Micro': pwd + "Trusted_GCSs_calling_0.01\EP_Micro_trusted_GCSs.txt",
                    'EP_Oxo': pwd + "Trusted_GCSs_calling_0.01\EP_Oxo_trusted_GCSs.txt"}

#Output data - GCSs, TAB (with score info added).
path_to_GCSs_sets_with_score=pwd+"Trusted_GCSs_h_s\\"
if not os.path.exists(path_to_GCSs_sets_with_score):
    os.makedirs(path_to_GCSs_sets_with_score)
Outpath_to_GCSs_files={'EP_Cfx': path_to_GCSs_sets_with_score + "EP_Cfx_10mkM_GCSs_trusted_h_s_0.01.txt",
                       'EP_RifCfx': path_to_GCSs_sets_with_score + "EP_RifCfx_GCSs_trusted_h_s_0.01.txt",
                       'EP_Micro': path_to_GCSs_sets_with_score + "EP_Micro_GCSs_trusted_h_s_0.01.txt",
                       'EP_Oxo': path_to_GCSs_sets_with_score + "EP_Oxo_GCSs_trusted_h_s_0.01.txt"}
#Output directory for N3E-score plots.
N3E_score_plot_path=path_to_GCSs_sets_with_score

#######
#Parsing and preparing score files.
#######

def Parsing_score(path_s):
    #Parse input WIG score file.
    filein=open(path_s, "r")
    score_ar=[]
    for line in filein:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            score_ar.append(float(line[0]))
    filein.close()
    return score_ar

#######
#Trusted GCSs data parsing.
#Calculate and return score for GCSs.
#Write GCSs files with scores info obtained.
#######

def GCSs_parsing_score_returning_info_writing(input_dict, score_ar, output_dict):
    GCSs_sets_dict={}
    for k, v in input_dict.items():
        GCSs_dict={}
        filein=open(v, 'r')
        fileout=open(output_dict[k], 'w')
        fileout.write('GCSs_coordinate\tN3E\tScore\n')
        for line in filein:
            line=line.rstrip().split('\t')
            if line[0] not in ['GCSs_coordinate']:
                GCSs_dict[int(line[0])]=[float(line[1]), (score_ar[int(line[0])-1+1]+score_ar[int(line[0])-1+4])/2]
                fileout.write(str(line[0]) + '\t' + str(line[1]) + '\t' + str((score_ar[int(line[0])-1+1]+score_ar[int(line[0])-1+4])/2) + '\n')
            else:
                continue
        filein.close()
        fileout.close()
        GCSs_sets_dict[k]=GCSs_dict
        print('Number of trusted GCSs for ' + str(k) + ' : ' + str(len(GCSs_dict)))
    return GCSs_sets_dict

#######
#Calculate N3E-score Pearson correlation.
#######

def correlation_h_s(GCSs_sets_dict):
    HS_dict={}
    for ab, GCSs_set in GCSs_sets_dict.items():
        #[[N3E data], [Score data]]
        HS_dict[ab]=[] 
        N3E=[]
        Score=[]
        for k, v in GCSs_set.items():
            N3E.append(v[0])
            Score.append(v[1])
        #Sorting.
        Score, N3E = (list(t) for t in zip(*sorted(zip(Score, N3E))))
        print('Paerson correlation (N3E, score) for ' + str(ab) + ' : ' + str(scipy.stats.pearsonr(N3E, Score)))
        N3E_log=np.log((np.array(N3E)))
        print('Paerson correlation (ln(N3E), score) for ' + str(ab) + ' : ' + str(scipy.stats.pearsonr(N3E_log, Score)))
        HS_dict[ab].append(N3E)
        HS_dict[ab].append(Score)
    return HS_dict

#######
#Visualize (Score, N3E) dependencies.
#######

def Plot_N3E_score(HS_dict, plot_path):
    #List dictionary keys to iterate over further.
    list_of_dict_keys=list(HS_dict.keys())
    print(list_of_dict_keys)
    
    #Generate array of colours.
    #Stolen from https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
    N=len(list_of_dict_keys)
    HSV_tuples=[(x*1.0/N, 0.5, 0.5) for x in range(N)]
    #print(HSV_tuples)
    RGB_tuples=list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
    #print(RGB_tuples)
    
    #Compute grid dimentions.  
    grid_dim=int(np.round(np.sqrt(len(HS_dict)), 0))
    if grid_dim*grid_dim>=len(HS_dict):
        grid_width=grid_dim
        grid_length=grid_dim
    else:
        grid_width=grid_dim
        grid_length=grid_dim+1  
    
    #Plot the data.   
    fig=plt.figure(figsize=(7.5*grid_width,7.5*grid_length), dpi=100)
    for j in range(grid_length):
        for i in range(grid_width):
            plot0=plt.subplot2grid((grid_length, grid_width),(j, i)) 
            index_var=((j*(grid_width)) + i)
            #print(index_var)
            if index_var<=len(list_of_dict_keys)-1:
                now_working_with=list_of_dict_keys[index_var]
            else:
                break
            fit=np.polyfit(HS_dict[now_working_with][1], HS_dict[now_working_with][0], 1)
            print(fit)
            fit_fn=np.poly1d(fit)  
            plot0.plot(HS_dict[now_working_with][1], HS_dict[now_working_with][0], 'o', fillstyle='none', color=RGB_tuples[index_var], markeredgecolor=RGB_tuples[index_var], markersize=2, alpha=0.8)
            plot0.plot(HS_dict[now_working_with][1], fit_fn(HS_dict[now_working_with][1]), '--k', label='y='+str(round(fit[0], 3))+'x+'+str(round(fit[1], 3)))
            pearson_corr=pearsonr(HS_dict[now_working_with][1], HS_dict[now_working_with][0])
            plot0.annotate(f'Pearson corr={np.round(pearson_corr[0],3)}\np-value={np.round(pearson_corr[1],5)}', xy=(0.5, 0.8), xycoords='axes fraction', size=20)
            plot0.set_xlabel('GCSs score', size=17) 
            plot0.set_ylabel('GCSs N3E', size=17)
            plot0.legend(loc='upper right', fontsize=17)
            plot0.set_title(now_working_with, size=18)
    plt.tight_layout()
    plt.savefig(plot_path + "GCSs_N3E_score_correlations.png", dpi=400, figsize=(7.5*grid_width,7.5*grid_length)) 
    plt.close()
    return
    

#######
#Visualize Score, N3E distributions.
#######

def Plot_N3E_score_distribution(HS_dict, score_ar, plot_path):
    #List dictionary keys to iterate over further.
    list_of_dict_keys=list(HS_dict.keys())
    print(list_of_dict_keys)
    
    #Generate array of colours.
    #Stolen from https://stackoverflow.com/questions/876853/generating-color-ranges-in-python
    N=len(list_of_dict_keys)
    HSV_tuples=[(x*1.0/N, 0.5, 0.5) for x in range(N)]
    #print(HSV_tuples)
    RGB_tuples=list(map(lambda x: colorsys.hsv_to_rgb(*x), HSV_tuples))
    #print(RGB_tuples)
    
    #Compute grid dimentions.  
    grid_dim=int(np.round(np.sqrt(len(HS_dict)), 0))
    if grid_dim*grid_dim>len(HS_dict):
        grid_width=grid_dim
        grid_length=grid_dim
    else:
        grid_width=grid_dim
        grid_length=grid_dim+1  
    
    #Plot the data.   
    fig=plt.figure(figsize=(7.5*grid_width,7.5*grid_length), dpi=100)
    for j in range(grid_length):
        for i in range(grid_width):
            index_var=((j*(grid_width)) + i)
            #print(index_var)
            if index_var<=len(list_of_dict_keys)-1:
                now_working_with=list_of_dict_keys[index_var]
            else:
                break            

            #N3E  
            plot0=plt.subplot2grid((grid_length*2, grid_width),(j*2, i)) 
            plot0.hist(HS_dict[now_working_with][0], color=RGB_tuples[index_var], edgecolor='black', alpha=0.8)
            plot0.annotate('Mean N3E='+str(round(np.mean(HS_dict[now_working_with][0]),2)), xy=(0.45, 0.9), xycoords='axes fraction', size=15)
            plot0.set_xlabel('GCSs N3E', size=17)
            plot0.set_ylabel('Number of GCSs', size=17)
            plot0.set_title(now_working_with, size=18)
            #Score
            plot1=plt.subplot2grid((grid_length*2, grid_width),((j*2)+1, i)) 
            plot1.hist(HS_dict[now_working_with][1], color=RGB_tuples[index_var], edgecolor='black', alpha=0.8)
            plot1.annotate('Mean score=\n'+str(round(np.mean(HS_dict[now_working_with][1]),2)), xy=(0.60, 0.8), xycoords='axes fraction', size=15)
            plot1.set_xlabel('GCSs score', size=17)
            plot1.set_ylabel('Number of GCSs', size=17)            
    
    #Whole genome scores distribution
    plot8=plt.subplot2grid((grid_length*2, grid_width),(j*2,i))  
    plot8.hist(score_ar, color='#5762ff', edgecolor='black', alpha=0.8)
    plot8.annotate('Mean score='+str(round(np.mean(score_ar),2)), xy=(0.65, 0.9), xycoords='axes fraction',  weight="bold", size=18)
    plot8.set_xlabel('Score', size=17) 
    plot8.set_ylabel('Number of genome positions', size=17)
    plot8.set_title('Scores distribution for $\it{E. coli}$ DY330 MuSGS genome', size=18)
    plt.tight_layout()
    plt.savefig(plot_path + "GCSs_N3E_score_distributions.png", dpi=400, figsize=(7.5*grid_width,7.5*grid_length)) 
    plt.close()     
    return


#######
#Wrapps all the functions together.
#######

def Wrapper(path_s, input_dict, output_dict, plot_path):
    score_track=Parsing_score(path_s)
    GCSs_info_dict=GCSs_parsing_score_returning_info_writing(input_dict, score_track, output_dict)
    Sorted_for_plot_dict=correlation_h_s(GCSs_info_dict)
    Plot_N3E_score(Sorted_for_plot_dict, plot_path)
    Plot_N3E_score_distribution(Sorted_for_plot_dict, score_track, plot_path)
    return
    
Wrapper(path_to_res_score_files, path_to_GCSs_files, Outpath_to_GCSs_files, N3E_score_plot_path)   

print('Script ended its work succesfully!')