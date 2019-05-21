###############################################
##Dmitry Sutormin, 2018##
##Topo-Seq analysis##

#Time-course gyrase experiment analysis.
#The script takes WIG-files with -IP-Cfx coverage, smoothes them and identifies the 
#replisome position.
###############################################

#######
#Packages to be imported.
#######

import os
from os import listdir
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numpy import diff
import scipy
from scipy import stats
from Bio import SeqIO
from Bio.SeqUtils import GC as GC_count
import matplotlib.patheffects as PathEffects
from fourier import *

#######
#Variables to be defined.
#######

#Deletions data.
Deletions="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Deletions_w3110_G_Mu_SGS.broadPeak"

#Folder with raw WIG files.
WIG_input_raw="F:\Gyrase_time-course_experiment\Reads_eq\WIG_files\Raw_data\\"
#Output folder for NSF data.
Output_directory="F:\Gyrase_time-course_experiment\Reads_eq\\"
#Folder with WIG files contains normalized and smoothed -IP+Cfx data.
WIG_control_NS_input="F:\Gyrase_time-course_experiment\Reads_eq\WIG_files\\NS\\"
#Folder with WIG files contains normalized and smoothed +IP+Cfx data.
WIG_exp_NS_input="F:\Gyrase_time-course_experiment\Reads_eq\WIG_files\\NS\\"

##Reviewed untill here


#Folder with WIG files (initial normalization, smoothing, filtering).
WIG_input_NSF="F:\Gyrase_time-course_experiment\WIG_files\Raw_data\-IP+Cfx_R1\\"
#Folder with WIG files (additional normalization on neutral region).
WIG_input_NN="F:\Gyrase_time-course_experiment\WIG_files\\"
#Folder with WIG files contains normalized and smoothed -IP+Cfx data.
WIG_control_NS_input="F:\Gyrase_time-course_experiment\WIG_files\\NS_NN\-IP+Cfx\R1\\"
#Folder with WIG files contains normalized, smoothed and filtered -IP+Cfx data.
WIG_control_NSF_input="F:\Gyrase_time-course_experiment\WIG_files\\NSF_NN\\"
#Folder with WIG files contains normalized and smoothed +IP+Cfx data.
WIG_exp_NS_input="F:\Gyrase_time-course_experiment\WIG_files\\NS_NN\+IP+Cfx\R1\\"
#Folder with WIG files contains normalized, smoothed and filtered +IP+Cfx data.
WIG_exp_NSF_input="F:\Gyrase_time-course_experiment\WIG_files\\NSF_NN\+IP+Cfx\R1\\"

#Folder with WIG files contains normalized, smoothed, filtered, and sequentially divided -IP+Cfx data.
WIG_cont_NSSDF_input="F:\Gyrase_time-course_experiment\WIG_files\Sequential_division\-IP+Cfx\R2\\NSF\\"
#Folder with WIG files contains normalized, smoothed, filtered (or not filtered), and sequentially divided +IP+Cfx data.
WIG_exp_NSSDF_input="F:\Gyrase_time-course_experiment\WIG_files\Sequential_division\\"

#Folder with WIG files contains normalized, smoothed, filtered, sequentially divided and differentiated -IP+Cfx data.
WIG_control_NSF_seq_div_dif_input="F:\Gyrase_time-course_experiment\WIG_files\Sequential_division\Dif\-IP+Cfx\\"

#Arrays plotting
WIG_input_1="F:\Gyrase_time-course_experiment\WIG_files\Sequential_division\-IP+Cfx\R1\\NSF\\"
WIG_input_2="F:\Gyrase_time-course_experiment\WIG_files\Sequential_division\-IP+Cfx\R2\\NSF\\"

#Output folder for Div data.
#Output_directory_div="C:\Sutor\science\DNA-gyrase\Results\E_coli_synch_time-course_Topo-Seq\Seq_results\Coverage_analysis\\"
Output_directory_div="F:\Gyrase_time-course_experiment\\"

#######
#Opens and reads BroadPeak file with deletions coordinates.
#Example:
#NC_007779.1_w3110_Mu\t274500\t372148\tDel_0\t10\t.\t1.0\t-1\t-1
#Which corresponds to GenomeID\tStart\tEnd\tFeature_name\tEnrichment\tStrand\tScore\tP-value\tE-value
#######

def deletions_info(del_path):
    del_ar=[]
    filein=open(del_path, 'r')
    for line in filein:
        line=line.rstrip().split('\t')
        del_ar.append([int(line[1]), int(line[2])])
    filein.close()
    return del_ar

#######
#Read directory with input files.
#######

def read_wig_files(wig_dir):
    input_files=listdir(wig_dir)
    wig_data_ar=[]
    for filename in input_files:
        wig_data=wig_parsing(wig_dir + filename)
        wig_data_ar.append(wig_data)
    return wig_data_ar


#######
#Parses WIG files.
#Computes a total number coverage.
#######

def wig_parsing(wigfile):
    print('Now is reading: ' + str(wigfile))
    wigin=open(wigfile, 'r')
    C_values=[]
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            C_values.append(float(line[0]))
    wigin.close()
    return C_values


#######
#Returns nearby NE value if current position falls into deleted region of the genome.
#######

def get_value(i, ends, deletions):
    if i<0: #coordinate is out of the left genome border (start)
        j=len(ends)+i
    elif i>=len(ends): #coordinate is out of the right genome border (end)
        j=i-len(ends)
    else: #coordinate is within the genome borders
        check_in_del=0
        for dl in deletions: #check if coordinate falls into deletion
            if dl[1]>=i>=dl[0]:
                j=dl[1]-dl[0]+i+1
                check_in_del=1
        if check_in_del==0:
            j=i
    return ends[j]


#######
#Smooth coverage tracks.
#######

def Smoothing(ends, deletions, window):
    smoothed=[]
    #Calculating the value for the first genome position
    mean=0.0
    window_float=float(window)
    for i in range(-window, window):
        mean=mean + get_value(i, ends, deletions)
    mean=mean/(2*window_float)
    smoothed.append(mean)
    #Calculating values for the part of the genome remains
    for i in range(1, len(ends)):
        mean=mean + (get_value(i+window, ends, deletions) - get_value(i-window, ends, deletions))/(2*window_float)
        smoothed.append(mean)
    return smoothed



#######
#Write the wig file.
#######

def write_wig(input_array, fileout_path, flag, strain_id, n_round):
    fileout=open(fileout_path, 'w+')
    fileout.write('track type=wiggle_0 name="'+str(flag)+'" autoScale=off viewLimits=0.0:25.0'+'\n'+'fixedStep chrom='+str(strain_id)+' start=1 step=1\n')   
    for i in range(len(input_array)):
        fileout.write(str(round(input_array[i], n_round))+'\n')
    fileout.close()    
    return


#######
#Wrapp functions for normalization, smoothing, filtering and plotting together.
#######

def NSF_wrapper(wig_dir, window, FF, n_round, del_path, strain_id, path_out, DS_name, ylim):
    #Read deletions info.
    deletions=deletions_info(del_path)
    #Read data.
    wig_data=read_wig_files(f'{wig_dir}\{DS_name}\\')
    
    #Identify minimal and maximal coverage.
    mean_coverage_ar=[]
    for sample in wig_data:
        mean_coverage_ar.append(np.mean(sample[0:2500000]))
    print(f'Mean coverages of neutral regions:\n {mean_coverage_ar}')
    Min_mean_coverage=min(mean_coverage_ar)
    print(f'Min mean coverage: {Min_mean_coverage}')  
    Max_mean_coverage=max(mean_coverage_ar)
    print(f'Max mean coverage: {Max_mean_coverage}')     
    #Plot distribution of mean coverages.
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle('Distribution of NR mean coverages', fontsize=20)
    plot1=plt.subplot()     
    plot1.hist(mean_coverage_ar, facecolor="green", alpha=0.7, linewidth=1, edgecolor='black')
    plot1.annotate(f'Min cov={round(Min_mean_coverage,2)}', xy=(0.60, 0.8), xycoords='axes fraction', size=15)
    plot1.annotate(f'Max cov={round(Max_mean_coverage,2)}', xy=(0.60, 0.7), xycoords='axes fraction', size=15)
    plot1.set_xlabel('NR mean coverage', size=17)
    plot1.set_ylabel('Number of samples', size=17)
    plt.show()  
      
    #Normalize samples coverage (add pseudocount +1).
    print("Samples are normalizing...")
    wig_data_norm=[]
    i=0
    for sample in wig_data:
        print(f'Now sample {i} out of {len(wig_data)-1} is normalizing... with mean_coverage: {mean_coverage_ar[i]}')
        sample_norm=[1.0 * (x + 1) * Min_mean_coverage/mean_coverage_ar[i] for x in sample]
        wig_data_norm.append(sample_norm)
        i+=1
        
    #Smooth normalized coverage.
    print("Samples are smoothing...")
    wig_data_norm_sm=[]
    i=0
    for sample_norm in wig_data_norm:
        print(f'Now sample {i} out of {len(wig_data)-1} is smoothing...')
        sample_norm_sm=Smoothing(sample_norm, deletions, window)
        wig_data_norm_sm.append(sample_norm_sm)
        i+=1
    
    #Data Fourier-filtration (optional).
    if FF[0]=='T':
        print('Samples are filtering...')
        #Rank of the last harmonic to be used for reconstruction.
        diap=FF[1]
        ff=Fourier()
        X=np.arange(0,len(wig_data_norm[0]))
        #Retained harmonics.
        retained_d={}
        for hnum in range(diap):
            retained_d[int(hnum)]=None 
        retained=list(retained_d)
        #Filtering.
        wig_data_norm_ff=[]
        i=0
        for sample_norm in wig_data_norm:
            print(f'Now sample {i} out of {len(wig_data)-1} is filtering...')
            ft=ff.fourier(np.array(sample_norm))
            ft_reduced=ff.reduce(ft, retained)
            sample_norm_filtered=ff.recover(X, ft_reduced)
            wig_data_norm_ff.append(sample_norm_filtered)
    elif FF[0]=='F':
        print('Fourier filtration is ommitting...')
        
    #Plot data.
    if not os.path.exists(f'{path_out}Figures\\NS'):
        os.makedirs(f'{path_out}Figures\\NS')
    if not os.path.exists(f'{path_out}Figures\\Raw'):
        os.makedirs(f'{path_out}Figures\\Raw')
    if FF[0]=='T':
        if not os.path.exists(f'{path_out}Figures\\NSF'):
            os.makedirs(f'{path_out}Figures\\NSF') 
    for i in range(len(wig_data)):
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print(f'Now sample {i} out of {len(wig_data)-1} is plotting...')
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        Xcor=np.arange(0,len(wig_data[i]))
        plot1=plt.subplot() 
        if FF[0]=='T': #If FF was performed.
            plt.suptitle(f'E. coli synchronization {DS_name} data smoothing {window}bp and \n F filtering {diap}h', fontsize=20)
            plot1.plot(Xcor, wig_data_norm[i], '-', label='Raw data N', color='black', linewidth=0.1)
            plot1.plot(Xcor, wig_data_norm_sm[i], '--', label=f'Data NSm {window}bp', color='#ee4e4e', linewidth=1)
            plot1.plot(Xcor, wig_data_norm_ff[i], '--', label=f'Data FF {diap}h', color='#D17E7E', linewidth=1)    
            plot1.set_ylim(ylim[0], ylim[1])
        elif FF[0]=='F': #If FF wasn't performed.
            plt.suptitle(f'E. coli synchronization {DS_name} data smoothing {window}bp', fontsize=20)
            plot1.plot(Xcor, wig_data_norm[i], '-', label='Raw data N', color='black', linewidth=0.1)
            plot1.plot(Xcor, wig_data_norm_sm[i], '--', label=f'Data NSm {window}bp', color='#ee4e4e', linewidth=1)  
            plot1.set_ylim(ylim[0], ylim[1])
        plot1.set_xlabel('Genome position, nt', size=17)
        plot1.set_ylabel('Coverage depth', size=17)
        plot1.legend(loc='upper left', fontsize=20)
        plt.show()
        if FF[0]=='T':
            plt.savefig(f'{path_out}Figures\\NSF\{DS_name}_{i}_N_{window}bp_S_{diap}h_FF.png', dpi=300, figsize=(16, 8))
        elif FF[0]=='F':
            plt.savefig(f'{path_out}Figures\\NS\{DS_name}_{i}_N_{window}bp_S.png', dpi=300, figsize=(16, 8))
        plt.close()
        
        #Raw data plotting (before N and S).
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        Xcor=np.arange(0,len(wig_data[i]))
        plot1=plt.subplot() 
        plt.suptitle(f'E. coli synchronization {DS_name} raw data', fontsize=20)
        plot1.plot(Xcor, wig_data[i], '-', label='Raw data', color='black', linewidth=0.1)
        plot1.annotate(f'Mean coverage NR={round(mean_coverage_ar[i],0)}', xy=(0.25, 0.95), xycoords='axes fraction', size=17)
        plot1.set_ylim(ylim[0], ylim[1])
        plot1.set_xlabel('Genome position, nt', size=17)
        plot1.set_ylabel('Coverage depth', size=17)
        plot1.legend(loc='upper left', fontsize=20)
        plt.show()
        plt.savefig(f'{path_out}Figures\\Raw\{DS_name}_{i}_raw.png', dpi=300, figsize=(16, 8))
        plt.close()
        
    #Write the WIG files.
    if not os.path.exists(f'{path_out}WIG_files\\N\{DS_name}'):
        os.makedirs(f'{path_out}WIG_files\\N\{DS_name}') 
    if not os.path.exists(f'{path_out}WIG_files\\NS\{DS_name}'):
        os.makedirs(f'{path_out}WIG_files\\NS\{DS_name}')
    if FF[0]=='T':
        if not os.path.exists(f'{path_out}WIG_files\\NF\{DS_name}'):
            os.makedirs(f'{path_out}WIG_files\\NF\{DS_name}')     
    for i in range(len(wig_data)):
        print(f'Now sample {i} out of {len(wig_data)-1} is writing...')
        write_wig(wig_data_norm[i], f'{path_out}WIG_files\\N\{DS_name}\\{DS_name}_{i}_N.wig', i, strain_id, n_round)
        write_wig(wig_data_norm_sm[i], f'{path_out}WIG_files\\NS\{DS_name}\\{DS_name}_{i}_N_{window}bp_S.wig', i, strain_id, n_round)
        if FF[0]=='T':
            write_wig(wig_data_norm_ff[i], f'{path_out}WIG_files\\NF\{DS_name}\\{DS_name}_{i}_N_{diap}h_FF.wig', i, strain_id, n_round)  
    return

NSF_wrapper(WIG_input_raw, 1000, ['F', 20], 15, Deletions, "NC_007779.1_w3110_Mu", Output_directory, '-IP+Cfx_R1', [-50, 1350])
NSF_wrapper(WIG_input_raw, 1000, ['F', 20], 15, Deletions, "NC_007779.1_w3110_Mu", Output_directory, '+IP+Cfx_R1', [-500, 8000])
NSF_wrapper(WIG_input_raw, 1000, ['F', 20], 15, Deletions, "NC_007779.1_w3110_Mu", Output_directory, '-IP+Cfx_R2', [-50, 1350])
NSF_wrapper(WIG_input_raw, 1000, ['F', 20], 15, Deletions, "NC_007779.1_w3110_Mu", Output_directory, '+IP+Cfx_R2', [-500, 8000])


#######
#Wrapp functions for normalization of +IP+Cfx samples by corresponding smoothed -IP+Cfx, 
#plotting and writing the data.
#######

def NPaired_wrapper(wig_dir_control, wig_dir_exp, n_round, strain_id, path_out, Replic):
    #Read data -Cfx+IP smoothed or filtered.
    wig_data_controls=read_wig_files(f'{wig_dir_control}\\-IP+Cfx_{Replic}') 
    #Read data +Cfx+IP smoothed or not.
    wig_data_exp=read_wig_files(f'{wig_dir_exp}\\+IP+Cfx_{Replic}')
    #Pairwise division +Cfx+IP/-Cfx+IP.
    wig_data_exp_div=[]
    for i in range(len(wig_data_controls)):
        exp_div=[]
        for j in range(len(wig_data_controls[i])):
            j_pos_div=wig_data_exp[i][j]/wig_data_controls[i][j]
            exp_div.append(j_pos_div)
        wig_data_exp_div.append(exp_div)
        
    #Data plotting and writing.
    if not os.path.exists(f'{path_out}Figures\Exp_NS_div_Cont_NS\\{Replic}'):
        os.makedirs(f'{path_out}Figures\Exp_NS_div_Cont_NS\\{Replic}')      
    if not os.path.exists(f'{path_out}WIG_files\Exp_NS_div_Cont_NS\\{Replic}'):
        os.makedirs(f'{path_out}WIG_files\Exp_NS_div_Cont_NS\\{Replic}')     
    X=np.arange(0,len(wig_data_exp_div[0]))
    ori_position=int((3712059+3711828)/2)
    dif_position=int((1593613+1593585)/2)     
    time_array=['-3min', '0min', '5min', '10min', '15min', '20min', '25min', '30min']
    for i in range(len(wig_data_exp_div)):
        print(f'Sample {i} out of {len(wig_data_exp_div)-1} is plotting...')
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit        
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        plt.suptitle("E. coli synchronization data +IP+Cfx/-IP+Cfx", fontsize=20)
        plot1=plt.subplot() 
        plot1.plot(X, wig_data_exp_div[i], '-', label=f'{time_array[i]} +IP+Cfx_NS/-IP+Cfx_NS', color='#550000', linewidth=1)
        plot1.plot(X, wig_data_controls[i], '-', label=f'{time_array[i]} -IP+Cfx NS', color='#D17E7E', linewidth=3) #NS: 600, 500
        plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
        plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
        plot1.set_xticks([dif_position, ori_position], minor=True)
        plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
        plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
        plot1.tick_params(axis='both', which='major', labelsize=17)  
        plot1.set_ylim(-0.5, 4) #NS: -0.5, 7;
        plot1.set_xlabel('Genome position, kb', size=20)
        plot1.set_ylabel('Coverage depth normalized', size=20)
        plot1.legend(loc='upper left', fontsize=25)
        plot1.annotate('', xytext=(ori_position, 0.0), xy=(ori_position, 0.350), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #+IP: 0, 200; -IP: 100, 200
        plt.show()
        plt.savefig(f'{path_out}Figures\Exp_NS_div_Cont_NS\\{Replic}\\{time_array[i]}_+IP+Cfx_{Replic}_NS_div_-IP+Cfx_{Replic}_NS_-IP+Cfx_{Replic}_NS.png', dpi=300, figsize=(16, 8))
        plt.close()
        write_wig(wig_data_exp_div[i], f'{path_out}WIG_files\Exp_NS_div_Cont_NS\\{Replic}\\{i}_+IP+Cfx_{Replic}_NS_div_-IP+Cfx_NS_{Replic}.wig', f'{i}_+IP+Cfx_{Replic}_NS/-IP+Cfx_{Replic}_NS', strain_id, n_round)
    return

#NPaired_wrapper(WIG_control_NS_input, WIG_exp_NS_input, 15, "NC_007779.1_w3110_Mu", Output_directory, 'R1')
#NPaired_wrapper(WIG_control_NS_input, WIG_exp_NS_input, 15, "NC_007779.1_w3110_Mu", Output_directory, 'R2')

##Reviewed untill here

#######
#Wrapp functions for plotting of sets of data of equial dimension (e.g. -IP+Cfx NS or NSF).
#######

def Coverage_plotting(wig_input_dir_NS, wig_input_dir_NSF, path_out):
    #Read input data NS.
    wig_data_NS=read_wig_files(wig_input_dir_NS) 
    #Read input data NSF.
    wig_data_NSF=read_wig_files(wig_input_dir_NSF)     
    #Data plotting.
    X=np.arange(0,len(wig_data_NS[0]))
    ori_position=int((3712059+3711828)/2)
    dif_position=int((1593613+1593585)/2)    
    time_array=['-3min', '0min', '5min', '10min', '15min', '20min', '25min', '30min']
    for i in range(len(wig_data_NS)):
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print("Sample is plotting...")
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        plt.suptitle("E. coli synchronization experiment data -IP+Cfx sd NS R1 and IP+Cfx sd NSF R2", fontsize=20)
        plot1=plt.subplot() 
        plot1.plot(X, wig_data_NS[i], '-', label=time_array[i] + ' -IP+Cfx_NSF R1', color='#550000', linewidth=3) #+IP NS: 0.2
        plot1.plot(X, wig_data_NSF[i], '-', label=time_array[i] + ' -IP+Cfx_NSF R2', color='#D17E7E', linewidth=3)
        plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
        plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
        plot1.set_xticks([dif_position, ori_position], minor=True)
        plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
        plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
        plot1.tick_params(axis='both', which='major', labelsize=17)  
        plot1.set_ylim(0.7, 1.3) #+IP: -500, 5500; -IP: 0, 1200
        plot1.set_xlabel('Genome position, kb', size=20)
        plot1.set_ylabel('Coverage depth normalized', size=20)
        plot1.legend(loc='upper left', fontsize=25)
        plot1.annotate('', xytext=(ori_position, 0.85), xy=(ori_position, 0.95), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #+IP: 0, 200; -IP: 100, 200
        plt.show()
        plt.savefig(path_out  + 'Figures\Sequential_division\-IP+Cfx\R1_R2_comp\\' + str(time_array[i]) + "_-IP+Cfx_R1_sd_NSF_vs_-IP+Cfx_R2_sd_NSF.png", dpi=300, figsize=(16, 8))
        plt.close()
    return

#Coverage_plotting(WIG_input_1, WIG_input_2, Output_directory)

#######
#Performes Fourier filtering of a list of arrays.
#######

def Fourier_filtering(wig_dir_in, diap, strain_id, name, path_out):
    #Read wig data.
    wig_data=read_wig_files(wig_dir_in) 
    #Data Fourier-filtration.
    print("Samples are filtering...")
    ff=Fourier()
    X=np.arange(0,len(wig_data[0]))
    #Harmonics to keep.
    retained_d={}
    for hnum in range(diap):
        retained_d[int(hnum)]=None 
    retained=list(retained_d)    
    #Data filtering.
    for i in range(len(wig_data)):
        #Fourier transformations.
        print("Now sample " + str(i+1) + " out of " + str(len(wig_data)) + " is filtering...")
        ft=ff.fourier(np.array(wig_data[i]))
        ft_reduced=ff.reduce(ft, retained)
        print(ft_reduced)
        wig_data_filtered=ff.recover(X, ft_reduced)
        print(X[:10])
        print(wig_data[i][:10])
        print(wig_data_filtered[:10])    
        write_wig(wig_data_filtered, path_out, name+"_filtered_"+str(diap)+"_fh", strain_id) 
    return

#######
#Wrapp functions for sequential division of time-points t(i+1)/t(i).
#######

def Seq_div_wrapper(wig_dir_control, wig_dir_exp, n_round, strain_id, path_out, prefix, suffix):
    #Read data -Cfx+IP.
    wig_data_controls=read_wig_files(wig_dir_control) 
    #Read data +Cfx+IP.
    wig_data_exp=read_wig_files(wig_dir_exp)
    #Sequental division t(i+1)/t(i).
    wig_data_seq_div=[]
    for i in range(len(wig_data_exp)-1): #For real sequential division t(i+1)/t(i).
    #for i in range(len(wig_data_exp)):  #For division on some fixed time-point (e.g. 10 min).
        seq_div=[]
        for j in range(len(wig_data_exp[i])):
            j_seq_div=wig_data_exp[i+1][j]/wig_data_exp[i][j] #For real sequential division t(i+1)/t(i).
            #j_seq_div=wig_data_exp[i][j]/wig_data_exp[4][j] #For division on some fixed time-point (e.g. 15 min).
            seq_div.append(j_seq_div)
        wig_data_seq_div.append(seq_div)
    #Data plotting.
    X=np.arange(0,len(wig_data_exp[0]))
    ori_position=int((3712059+3711828)/2)
    dif_position=int((1593613+1593585)/2)    
    time_array=['-3min', '0min', '5min', '10min', '15min', '20min', '25min', '30min']    
    for i in range(len(wig_data_seq_div)):       
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print("Sample is plotting...")
        #mpl.rcParams['agg.path.chunksize']=10000
        #plt.figure(figsize=(16, 8), dpi=100)
        #plt.suptitle("E. coli synchronization data sequential division t(i+1)/t(i)", fontsize=20)
        #plot1=plt.subplot() 
        #plot1.plot(X, wig_data_seq_div[i], '-', label=time_array[i+1] + '/' + time_array[i] + ' -IP+Cfx_NSF', color='#550000', linewidth=1) #+IP NS: 0.2
        #plot1.plot(X, np.array(wig_data_controls[i])/450, '--', label=time_array[i] + ' -IP+Cfx_NSF', color='#D17E7E', linewidth=2) #+IP NS R1: 500; +IP NS R1: 450;
        #plot1.plot(X, np.array(wig_data_controls[i+1])/450, '--', label=time_array[i+1] + ' -IP+Cfx_NSF', color='#882D60', linewidth=2) #+IP NS R1: 500; +IP NS R1: 450;
        #plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
        #plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
        #plot1.set_xticks([dif_position, ori_position], minor=True)
        #plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
        #plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
        #plot1.tick_params(axis='both', which='major', labelsize=17)  
        #plot1.set_ylim(0.0, 2.5) #+IP: ; -IP: 
        #plot1.set_xlabel('Genome position, kb', size=20)
        #plot1.set_ylabel('Coverage depth normalized', size=20)
        #plot1.legend(loc='upper left', fontsize=25)
        #plot1.annotate('', xytext=(ori_position, 0.3), xy=(ori_position, 0.550), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #+IP: ; -IP: 
        #plt.show()
        #plt.savefig(path_out  + "Figures\Sequential_division\-IP+Cfx\R2\\" + str(i+1) + "_div_" +str(i)+ "_-IP+Cfx_NSF_R2.png", dpi=300, figsize=(16, 8))
        #plt.close()
        #Write the WIG files.
        write_wig(wig_data_seq_div[i], path_out + prefix + str(i+1) + "_div_" + str(i) + suffix, str(i+1)+"/"+str(i), strain_id, n_round)  #NS or NSF in case of +IP+Cfx     
    return

#Seq_div_wrapper(WIG_control_NSF_input+"-IP+Cfx\R1\\", WIG_control_NSF_input+"-IP+Cfx\R1\\", 15, "NC_007779.1_w3110_Mu", Output_directory_div, "WIG_files\Sequential_division\-IP+Cfx\R1\\NSF\\", "_-IP+Cfx_NSF_R1.wig")

#######
#Wrapp functions for sequential deduction of time-points t(i+1)-t(i).
#######

def Seq_ded_wrapper(wig_dir_control, wig_dir_exp, diap, del_path, strain_id, path_out):
    #Read deletions info.
    deletions=deletions_info(del_path)
    #Read data -Cfx+IP smoothed and filtered.
    wig_data_controls=read_wig_files(wig_dir_control) 
    #Read data +Cfx+IP smoothed.
    wig_data_exp=read_wig_files(wig_dir_exp)
    #Pairwise division +Cfx+IP/-Cfx+IP.
    wig_data_exp_div=[]
    for i in range(len(wig_data_controls)):
        exp_div=[]
        for j in range(len(wig_data_controls[i][0])):
            j_pos_div=wig_data_exp[i][0][j]/wig_data_controls[i][0][j]
            exp_div.append(j_pos_div)
        wig_data_exp_div.append(exp_div) 
    #Sequental deduction t(i+1)-t(i)
    wig_data_seq_ded=[]
    for i in range(len(wig_data_exp_div)-1):
        seq_ded=[]
        for j in range(len(wig_data_exp_div[i])):
            j_seq_ded=wig_data_exp_div[i+1][j]-wig_data_exp_div[i][j]
            seq_ded.append(j_seq_ded)
        wig_data_seq_ded.append(seq_ded)
        
    #Data Fourier-filtration.
    print("Samples are filtering...")
    ff=Fourier()
    X=np.arange(0,len(wig_data_seq_ded[0]))
    
    retained_d={}
    for hnum in range(diap):
        retained_d[int(hnum)]=None 
    retained=list(retained_d)
    
    #Data plotting.
    for i in range(len(wig_data_seq_ded)):
        #Fourier transformations.
        print("Now sample " + str(i+1) + " out of " + str(len(wig_data_seq_ded)) + " is filtering...")
        ft=ff.fourier(np.array(wig_data_seq_ded[i]))
        ft_reduced=ff.reduce(ft, retained)
        print(ft_reduced)
        wig_data_seq_ded_filtered=ff.recover(X, ft_reduced)
        print(X[:10])
        print(wig_data_seq_ded[i][:10])
        print(wig_data_seq_ded_filtered[:10])        
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print("Sample is plotting...")
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        plt.suptitle("E. coli synchronization data sequential deduction t(i+1)-t(i)", fontsize=20)
        plot1=plt.subplot() 
        plot1.plot(X, wig_data_seq_ded[i], '-', label='(+IP+Cfx_NS/-IP+Cfx_NSF)i+1 - (+IP+Cfx_NS/-IP+Cfx_NSF)i Original', color='#550000', linewidth=0.2)
        plot1.plot(X, wig_data_seq_ded_filtered, '--', label='(+IP+Cfx_NS/-IP+Cfx_NSF)i+1 - (+IP+Cfx_NS/-IP+Cfx_NSF)i filtered', color='#D17E7E', linewidth=3)
        plot1.set_xlabel('Genome position, nt', size=17)
        plot1.set_ylabel('Coverage depth normalized', size=17)
        plot1.legend(loc='upper left')
        plt.show()
        plt.savefig(path_out  + "Figures\Sequential_ded\\" + str(i+1) + "_ded_" +str(i)+ "_plus_20fh_R1.png", dpi=300, figsize=(16, 8))
        plt.close()
        #Write the WIG files.
        write_wig(wig_data_seq_ded[i], path_out + "WIG_files\Seq_ded\\" + str(i+1) + "_ded_" + str(i) + "_R1.wig", str(i+1)+"/"+str(i), strain_id)
        write_wig(wig_data_seq_ded_filtered, path_out + "WIG_files\Seq_ded\\" + str(i+1) + "_ded_" + str(i) + "_20fh_R1.wig", str(i+1)+"/"+str(i), strain_id)        
    return

#Seq_ded_wrapper(WIG_control_NSF_input, WIG_exp_NS_input, 20, Deletions, "NC_007779.1_w3110_Mu", Output_directory_div)

#######
#Wrapp functions for sequential division of time-points t(i+1)/t(i) for controls NSF.
#######

def Seq_div_or_ded_of_controls_wrapper(wig_dir_control, wig_dir_exp, diap, del_path, strain_id, path_out):
    #Read deletions info.
    deletions=deletions_info(del_path)
    #Read data -Cfx+IP smoothed and filtered.
    wig_data_controls=read_wig_files(wig_dir_control) 
    #Sequental division t(i+1)/t(i) or deduction t(i+1)-t(i).
    wig_data_controls_seq_div=[]
    for i in range(len(wig_data_controls)-1): #For real sequential division t(i+1)/t(i) or deduction t(i+1)-t(i).
    #for i in range(len(wig_data_controls)):  #For division on some fixed time-point (e.g. 10 min).
        seq_div=[]
        for j in range(len(wig_data_controls[i][0])):
            #j_seq_div=wig_data_controls[i+1][0][j]/wig_data_controls[i][0][j] #For real sequential division t(i+1)/t(i).
            j_seq_div=wig_data_controls[i+1][0][j]-wig_data_controls[i][0][j] #For real sequential deduction t(i+1)-t(i).
            #j_seq_div=wig_data_controls[i][0][j]/wig_data_controls[4][0][j] #For division on some fixed time-point (e.g. 15 min).
            seq_div.append(j_seq_div)
        wig_data_controls_seq_div.append(seq_div)
        
    X=np.arange(0,len(wig_data_controls_seq_div[0]))
    
    #Data plotting.
    for i in range(len(wig_data_controls_seq_div)):     
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print("Sample is plotting...")
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        #plt.suptitle("E. coli synchronization -IP+Cfx data sequential division t(i+1)/t(i)", fontsize=20) #Div
        plt.suptitle("E. coli synchronization -IP+Cfx data sequential deduction t(i+1)-t(i)", fontsize=20) #Ded
        plot1=plt.subplot() 
        #plot1.plot(X, wig_data_controls_seq_div[i]/np.mean(wig_data_controls_seq_div[i][0:2500000]), '-', label='(-IP+Cfx_NSF)i+1/(-IP+Cfx_NSF)i', color='#D17E7E', linewidth=3) #Div
        plot1.plot(X, wig_data_controls_seq_div[i]+np.abs(np.mean(wig_data_controls_seq_div[i][0:2500000])), '-', label='(-IP+Cfx_NSF)i+1-(-IP+Cfx_NSF)i', color='#D17E7E', linewidth=3) #Ded
        plot1.set_xlabel('Genome position, nt', size=17)
        plot1.set_ylabel('Coverage depth normalized', size=17)
        #plot1.set_ylim(0.85, 1.35) #Div
        plot1.set_ylim(-20, 150) #Ded
        plot1.legend(loc='upper left')
        plt.show()
        #plt.savefig(path_out  + "Figures\Sequential_ded\-IP+Cfx\R1\\" + str(i+1) + "_div_" +str(i)+ "_-IP+Cfx_NSF_R1.png", dpi=300, figsize=(16, 8)) #Div
        plt.savefig(path_out  + "Figures\Sequential_ded\-IP+Cfx\R2\\" + str(i+1) + "_ded_" +str(i)+ "_-IP+Cfx_NSF_R2.png", dpi=300, figsize=(16, 8)) #Ded
        plt.close()
        #Write the WIG files.
        #write_wig(wig_data_controls_seq_div[i], path_out + "WIG_files\Seq_ded\-IP+Cfx\R1\\" + str(i+1) + "_div_" + str(i) + "_-IP+Cfx_NSF_R1.wig", str(i+1)+"/"+str(i), strain_id) #Div
        write_wig(wig_data_controls_seq_div[i], path_out + "WIG_files\Seq_ded\-IP+Cfx\R2\\" + str(i+1) + "_ded_" + str(i) + "_-IP+Cfx_NSF_R2.wig", str(i+1)+"/"+str(i), strain_id) #Ded 
    return

#Seq_div_or_ded_of_controls_wrapper(WIG_control_NSF_input, WIG_exp_NS_input, 20, Deletions, "NC_007779.1_w3110_Mu", Output_directory_div)

#######
#Wrapp functions for camparison of changes in gyrase enrichment (+IP+Cfx div or ded) with coverage growth (-IP+Cfx NSF).
#######

def Gyr_div_or_ded_vs_contr_NSF(wig_dir_control_NFS, wig_dir_exp_seq_div, diap, del_path, strain_id, path_out):
    #Read deletions info.
    deletions=deletions_info(del_path)
    #Read data -Cfx+IP smoothed and filtered.
    wig_data_controls_NFS=read_wig_files(wig_dir_control_NFS) 
    #Read data +Cfx+IP sequential divided or deduced.
    wig_data_exp_SD=read_wig_files(wig_dir_exp_seq_div)  
    
    #Plot NFS control and SD experimental data.
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    
    #plot1.plot(xcoord, np.array(wig_data_controls_NFS[0][0])/550, ':', label='-3 min NFS', color='#550000', linewidth=2)
    #plot1.plot(xcoord, wig_data_exp_SD[0][0], '--', label='0 min/-3 min SD', color='#801515', linewidth=2)
    #plot1.plot(xcoord, np.array(wig_data_controls_NFS[1][0])/550, ':', label='0 min NFS', color='#801515', linewidth=2)
    #plot1.plot(xcoord, wig_data_exp_SD[1][0], '--', label='5 min/0 min SD', color='#AA3939', linewidth=2)
    #plot1.plot(xcoord, np.array(wig_data_controls_NFS[2][0])/550, ':', label='5 min NFS', color='#AA3939', linewidth=2)
    
    #For derivative calculation.
    dx=1
    xcoord=np.arange(0,4647999)
    dxcoord=np.arange(0,4647998)
    
    #-3 min to 0 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot()     
    plot1.plot(xcoord, [x+1 for x in wig_data_exp_SD[0][0]], '--', label='-3 min/10 min SD', color='#AA3939', linewidth=2)
    NFSN=np.array(wig_data_controls_NFS[1][0])/np.mean(wig_data_controls_NFS[1][0][0:2500000])
    plot1.plot(xcoord, NFSN, ':', label='0 min NFS', color='#A8658B', linewidth=2) #550 R1
    dNFSN=diff(NFSN)/dx
    plot1.plot(dxcoord, dNFSN, ':', label='0 min dNFS', color='#69073E', linewidth=1) #550 R1
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([int((1593613+1593585)/2), int((3712059+3711828)/2)], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    plot1.set_ylim(0.6, 2)
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    plot1.annotate('', xytext=(int((3712059+3711828)/2), 0.7), xy=(int((3712059+3711828)/2), 0.8), arrowprops=dict(arrowstyle='->', facecolor='green'), color='black', weight="bold", size=15)
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_ded_vs_NSF_vs_dNSF\\" + "R1_dNSF_0_min.png", dpi=300, figsize=(16, 8))
    plt.close()      
    
    #0 min to 5 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot()     
    plot1.plot(xcoord, [x+1 for x in wig_data_exp_SD[1][0]], '--', label='0 min/10 min SD', color='#AA3939', linewidth=2)
    NFSN=np.array(wig_data_controls_NFS[2][0])/np.mean(wig_data_controls_NFS[2][0][0:2500000])
    #plot1.plot(xcoord, NFSN, ':', label='5 min NFS', color='#A8658B', linewidth=2) #550 R1
    dNFSN=diff(NFSN)/dx
    plot1.plot(dxcoord, dNFSN, ':', label='5 min dNFS', color='#69073E', linewidth=1) #550 R1
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([int((1593613+1593585)/2), int((3712059+3711828)/2)], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    #plot1.set_ylim(0.6, 2)
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #plot1.annotate('', xytext=(int((3712059+3711828)/2), 0.7), xy=(int((3712059+3711828)/2), 0.8), arrowprops=dict(arrowstyle='->', facecolor='green'), color='black', weight="bold", size=15)
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_ded_vs_NSF_vs_dNSF\\" + "R1_dNSF_5_min.png", dpi=300, figsize=(16, 8))
    plt.close()         
    
    #5 min to 10 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot()     
    #plot1.plot(xcoord, [x+1 for x in wig_data_exp_SD[2][0]], '--', label='10 min/5 min SDed', color='#AA3939', linewidth=2)
    NFSN=np.array(wig_data_controls_NFS[3][0])/np.mean(wig_data_controls_NFS[3][0][0:2500000])
    #plot1.plot(xcoord, NFSN, ':', label='10 min NFS', color='#A8658B', linewidth=2) #550 R1
    dNFSN=diff(NFSN)/dx
    plot1.plot(dxcoord, dNFSN, ':', label='10 min dNFS', color='#69073E', linewidth=1) #550 R1
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([int((1593613+1593585)/2), int((3712059+3711828)/2)], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    #plot1.set_ylim(0.6, 2)
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #plot1.annotate('', xytext=(int((3712059+3711828)/2), 0.7), xy=(int((3712059+3711828)/2), 0.8), arrowprops=dict(arrowstyle='->', facecolor='green'), color='black', weight="bold", size=15)
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_ded_vs_NSF_vs_dNSF\\" + "R1_dNSF_10_min.png", dpi=300, figsize=(16, 8))
    plt.close()     
    
    #10 min to 15 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot() 
    #plot1.plot(xcoord, [x+1 for x in wig_data_exp_SD[3][0]], '--', label='15 min/10 min SDed', color='#AA3939', linewidth=2)
    NFSN=np.array(wig_data_controls_NFS[4][0])/np.mean(wig_data_controls_NFS[4][0][0:2500000])
    #plot1.plot(xcoord, NFSN, ':', label='15 min NFS', color='#A8658B', linewidth=2) #550 R1; 
    dNFSN=diff(NFSN)/dx
    plot1.plot(dxcoord, dNFSN, ':', label='15 min dNFS', color='#69073E', linewidth=1) #550 R1
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([int((1593613+1593585)/2), int((3712059+3711828)/2)], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    #plot1.set_ylim(0.6, 2)
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    #plot1.annotate('', xytext=(int((3712059+3711828)/2), 0.7), xy=(int((3712059+3711828)/2), 0.8), arrowprops=dict(arrowstyle='->', facecolor='green'), color='black', weight="bold", size=15)
    #Mark local maxima.
    #max_index1=wig_data_exp_SD[3][0][3500000:4000000].index(max(wig_data_exp_SD[3][0][3500000:4000000]))+3500000
    #plot1.annotate('', xytext=(max_index1, 1.4), xy=(max_index1, 1.35), arrowprops=dict(arrowstyle='->', facecolor='black'), color='black', weight="bold", size=15)
    #Title axis.
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_ded_vs_NSF_vs_dNSF\\" + "R1_dNSF_15_min.png", dpi=300, figsize=(16, 8))
    plt.close()     
    
    #15 min to 20 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot() 
    #plot1.plot(xcoord, [x+1 for x in wig_data_exp_SD[4][0]], '--', label='20 min/15 min SDed', color='#AA3939', linewidth=2)
    NFSN=np.array(wig_data_controls_NFS[5][0])/np.mean(wig_data_controls_NFS[5][0][0:2500000])
    #plot1.plot(xcoord, NFSN, ':', label='20 min NFS', color='#A8658B', linewidth=2) #550 R1; 
    dNFSN=diff(NFSN)/dx
    plot1.plot(dxcoord, dNFSN, ':', label='20 min dNFS', color='#69073E', linewidth=1) #550 R1
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([int((1593613+1593585)/2), int((3712059+3711828)/2)], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    #plot1.set_ylim(0.6, 2)
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    #plot1.annotate('', xytext=(int((3712059+3711828)/2), 0.7), xy=(int((3712059+3711828)/2), 0.8), arrowprops=dict(arrowstyle='->', facecolor='green'), color='black', weight="bold", size=15)
    #Mark local maxima.
    #max_index1=wig_data_exp_SD[4][0][3300000:3800000].index(max(wig_data_exp_SD[4][0][3300000:3800000]))+3300000
    #plot1.annotate('', xytext=(max_index1, 1.4), xy=(max_index1, 1.35), arrowprops=dict(arrowstyle='->', facecolor='black'), color='black', weight="bold", size=15)
    #max_index2=wig_data_exp_SD[4][0][3800000:4200000].index(max(wig_data_exp_SD[4][0][3800000:4200000]))+3800000
    #plot1.annotate('', xytext=(max_index2, 1.4), xy=(max_index2, 1.35), arrowprops=dict(arrowstyle='->', facecolor='black'), color='black', weight="bold", size=15)    
    #Title axis.    
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_ded_vs_NSF_vs_dNSF\\" + "R1_dNSF_20_min.png", dpi=300, figsize=(16, 8))
    plt.close()     
    
    #20 min to 25 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot() 
    #plot1.plot(xcoord, [x+1 for x in wig_data_exp_SD[5][0]], '--', label='25 min/20 min SDed', color='#AA3939', linewidth=2)
    NFSN=np.array(wig_data_controls_NFS[6][0])/np.mean(wig_data_controls_NFS[6][0][0:2500000])
    #plot1.plot(xcoord, NFSN, ':', label='25 min NFS', color='#A8658B', linewidth=2) #500 R1; 
    dNFSN=diff(NFSN)/dx
    plot1.plot(dxcoord, dNFSN, ':', label='25 min dNFS', color='#69073E', linewidth=1) #550 R1
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([int((1593613+1593585)/2), int((3712059+3711828)/2)], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    #plot1.set_ylim(0.6, 2)
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    #plot1.annotate('', xytext=(int((3712059+3711828)/2), 0.7), xy=(int((3712059+3711828)/2), 0.8), arrowprops=dict(arrowstyle='->', facecolor='green'), color='black', weight="bold", size=15)
    #Mark local maxima.
    #max_index1=wig_data_exp_SD[5][0][3000000:3500000].index(max(wig_data_exp_SD[5][0][3000000:3500000]))+3000000
    #plot1.annotate('', xytext=(max_index1, 1.4), xy=(max_index1, 1.35), arrowprops=dict(arrowstyle='->', facecolor='black'), color='black', weight="bold", size=15)
    #max_index2=wig_data_exp_SD[5][0][4000000:4500000].index(max(wig_data_exp_SD[5][0][4000000:4500000]))+4000000
    #plot1.annotate('', xytext=(max_index2, 1.4), xy=(max_index2, 1.35), arrowprops=dict(arrowstyle='->', facecolor='black'), color='black', weight="bold", size=15)    
    #Title axis.    
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_ded_vs_NSF_vs_dNSF\\" + "R1_dNSF_25_min.png", dpi=300, figsize=(16, 8))
    plt.close() 
    
    #25 min to 30 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot() 
    #plot1.plot(xcoord, [x+1 for x in wig_data_exp_SD[6][0]], '--', label='30 min/25 min SDed', color='#AA3939', linewidth=2)  
    NFSN=np.array(wig_data_controls_NFS[7][0])/np.mean(wig_data_controls_NFS[7][0][0:2500000])
    #plot1.plot(xcoord, NFSN, ':', label='30 min NFS', color='#A8658B', linewidth=2) #480 R1
    dNFSN=diff(NFSN)/dx
    plot1.plot(dxcoord, dNFSN, ':', label='30 min dNFS', color='#69073E', linewidth=1) #550 R1
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([int((1593613+1593585)/2), int((3712059+3711828)/2)], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    #plot1.set_ylim(0.6, 2)
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    #plot1.annotate('', xytext=(int((3712059+3711828)/2), 0.7), xy=(int((3712059+3711828)/2), 0.8), arrowprops=dict(arrowstyle='->', facecolor='green'), color='black', weight="bold", size=15)
    #Mark local maxima.
    #max_index1=wig_data_exp_SD[6][0][2700000:3500000].index(max(wig_data_exp_SD[6][0][2700000:3500000]))+2700000
    #plot1.annotate('', xytext=(max_index1, 1.4), xy=(max_index1, 1.35), arrowprops=dict(arrowstyle='->', facecolor='black'), color='black', weight="bold", size=15)
    #max_index2=wig_data_exp_SD[6][0][4200000:4600000].index(max(wig_data_exp_SD[6][0][4200000:4600000]))+4200000
    #plot1.annotate('', xytext=(max_index2, 1.4), xy=(max_index2, 1.35), arrowprops=dict(arrowstyle='->', facecolor='black'), color='black', weight="bold", size=15)    
    #Title axis.      
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_ded_vs_NSF_vs_dNSF\\" + "R1_dNSF_30_min.png", dpi=300, figsize=(16, 8))
    plt.close()        
    return

#Gyr_div_or_ded_vs_contr_NSF(WIG_control_NSF_input, WIG_exp_NSSDF_input, 20, Deletions, "NC_007779.1_w3110_Mu", Output_directory_div)

#######
#Computes derivative of an array.
#######

def Der_comp(wig_dir_in, dx, n_round, strain_id, path_out, prefix, suffix):
    #Read wig data.
    wig_data=read_wig_files(wig_dir_in) 
    #Data differentiating.
    print("Samples are differentiating...")  
    for i in range(len(wig_data)):
        #Derivative calculation.
        print("Now sample " + str(i+1) + " out of " + str(len(wig_data)) + " is differentiating...")        
        dif_array=diff(np.array(wig_data[i]))/dx
        write_wig(dif_array, path_out + prefix + str(i+1) + "_div_" + str(i) + suffix, str(i+1)+"/"+str(i), strain_id, n_round)       
    return 

#Der_comp(WIG_input_1, 0.000001, 15, "NC_007779.1_w3110_Mu", Output_directory_div, "WIG_files\Sequential_division\Dif\-IP+Cfx\R1\\NSF\\", "_-IP+Cfx_NSF_R1_dif.wig")
#Der_comp(WIG_input_2, 0.000001, 15, "NC_007779.1_w3110_Mu", Output_directory_div, "WIG_files\Sequential_division\Dif\-IP+Cfx\R2\\NSF\\", "_-IP+Cfx_NSF_R2_dif.wig")

#######
#Wrapp functions for plotting of sets of data of equial dimension (e.g. -IP+Cfx NSF, d((-IP+Cfx NSF)i+1/(-IP+Cfx NSF)i), (+IP+Cfx)i+1/(+IP+Cfx)i).
#######

def Coverage_plotting_locate_gyr_and_rep(wig_input_dir_NSF, wig_input_seq_div_dif_NSF, wig_exp_seq_div_NSF, path_out, prefix, suffix):
    #Read input data NSF.
    wig_data_NSF=read_wig_files(wig_input_dir_NSF) 
    #Read input data NSF sequentially divided and differentiated.
    wig_data_NSF_sdd=read_wig_files(wig_input_seq_div_dif_NSF)   
    #Read experimental data NSF sequentially divided.
    wig_data_NSF_exp_sd=read_wig_files(wig_exp_seq_div_NSF)     
    #Data plotting.
    X=np.arange(0,len(wig_data_NSF[0]))
    dX=np.arange(0,len(wig_data_NSF_sdd[0]))
    ori_position=int((3712059+3711828)/2)
    dif_position=int((1593613+1593585)/2)    
    time_array=['-3min', '0min', '5min', '10min', '15min', '20min', '25min', '30min']
    for i in range(len(wig_data_NSF_exp_sd)):
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print("Sample is plotting...")
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        plt.suptitle("E. coli synchronization experiment data +IP+Cfx sd NS, -IP+Cfx NSF, -IP+Cfx sdd NSF", fontsize=20)
        plot1=plt.subplot() 
        plot1.plot(X, np.array(wig_data_NSF[i])/400, '-', label=time_array[i] + ' -IP+Cfx_NSF', color='#D17E7E', linewidth=2) #-IP+Cfx coverage t(i); 500 R1
        plot1.plot(X, np.array(wig_data_NSF[i+1])/400, '-', label=time_array[i+1] + ' -IP+Cfx_NSF', color='#882D60', linewidth=2) #-IP+Cfx coverage t(i+1); 500 R1
        plot1.plot(X, wig_data_NSF_exp_sd[i], '-', label=time_array[i+1] + '/' + time_array[i] + ' +IP+Cfx_NS', color='#550000', linewidth=0.2) #+IP+Cfx coverage t(i+1)/t(i)
        plot1.plot(0.0, 1.0, ':', label='d(' + time_array[i+1] + '/' + time_array[i] + ' -IP+Cfx_NSF)', color='black', linewidth=1) #Legend for -IP+Cfx dif of t(i+1)/t(i)
        plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
        plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
        plot1.set_xticks([dif_position, ori_position], minor=True)
        plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
        plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
        plot1.tick_params(axis='both', which='major', labelsize=17)  
        plot1.set_ylim(0.4, 2.2) #+IP: -500, 5500; -IP: 0, 1200
        plot1.set_xlabel('Genome position, kb', size=20)
        plot1.set_ylabel('Coverage depth normalized', size=20)
        plot1.legend(loc='upper left', fontsize=25)
        plot1.annotate('', xytext=(ori_position, 0.55), xy=(ori_position, 0.65), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #+IP: 0, 200; -IP: 100, 200
        plot2=plot1.twinx()
        plot2.plot(dX, np.array(wig_data_NSF_sdd[i]), ':', label='d(' + time_array[i+1] + '/' + time_array[i] + ' -IP+Cfx_NSF)', color='black', linewidth=1) #-IP+Cfx dif of t(i+1)/t(i)
        plot2.set_ylabel('Derivative of (-IP+Cfx NSF)i+1/(-IP+Cfx NSF)i', size=20)
        plot2.tick_params(axis='y', which='major', labelsize=17) 
        plot2.set_ylim(-1.5, 2.5)
        plot2.set_yticks([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5], minor=False)
        plot2.set_yticklabels(np.array([-1.5,-1.0,-0.5,0.0,0.5,1.0,1.5,2.0,2.5])/1000000, minor=False)        
        plt.show()
        plt.savefig(path_out  + prefix + str(i+1) + '_' + str(i) + suffix, dpi=300, figsize=(16, 8)) #, format='eps'
        plt.close()
        print('Position of the left replisome at moment ' + time_array[i+1] + ': ' + str(wig_data_NSF_sdd[i].index(max(wig_data_NSF_sdd[i]))))
        print('Position of the right replisome at moment ' + time_array[i+1] + ': ' + str(wig_data_NSF_sdd[i].index(min(wig_data_NSF_sdd[i]))))
    return

#Coverage_plotting_locate_gyr_and_rep(WIG_control_NSF_input+"-IP+Cfx\R1\\", WIG_control_NSF_seq_div_dif_input+"R1\\NSF\\", WIG_exp_NSSDF_input+"+IP+Cfx\R1\\NSF\\", Output_directory, 'Figures\Dif_of_cont_vs_cont_vs_exp\R1\\NSF\\', "_-IP+Cfx_NSF_vs_-IP+Cfx_sdd_NSF_vs_+IP+Cfx_sd_NSF_R1.eps")
#Coverage_plotting_locate_gyr_and_rep(WIG_control_NSF_input+"-IP+Cfx\R2\\", WIG_control_NSF_seq_div_dif_input+"R2\\NSF\\", WIG_exp_NSSDF_input+"+IP+Cfx\R2\\NS\\", Output_directory, 'Figures\Dif_of_cont_vs_cont_vs_exp\R2\\NS\\', "_-IP+Cfx_NSF_vs_-IP+Cfx_sdd_NSF_vs_+IP+Cfx_sd_NS_R2.png")

#######
#Wrapp functions for replisome (Seq div or ded -IP+Cfx) and gyrase (Seq div or ded +IP+Cfx) enrichment positions identification.
#######

def Repl_Gyr_ident(wig_dir_control_seq_ch, wig_dir_exp_seq_ch, diap, del_path, strain_id, path_out):
    #Read deletions info.
    deletions=deletions_info(del_path)
    #Read data -Cfx+IP seuentially divided ar deduced.
    wig_data_controls_SD=read_wig_files(wig_dir_control_seq_ch) 
    #Read data +Cfx+IP sequentially divided or deduced.
    wig_data_exp_SD=read_wig_files(wig_dir_exp_seq_ch)  
    
    #Plot NFS control and SD experimental data.
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    
    #plot1.plot(xcoord, np.array(wig_data_controls_SD[0][0])/550, ':', label='-3 min NFS', color='#550000', linewidth=2)
    #plot1.plot(xcoord, wig_data_exp_SD[0][0], '--', label='0 min/-3 min SD', color='#801515', linewidth=2)
    #plot1.plot(xcoord, np.array(wig_data_controls_SD[1][0])/550, ':', label='0 min NFS', color='#801515', linewidth=2)
    #plot1.plot(xcoord, wig_data_exp_SD[1][0], '--', label='5 min/0 min SD', color='#AA3939', linewidth=2)
    #plot1.plot(xcoord, np.array(wig_data_controls_SD[2][0])/550, ':', label='5 min NFS', color='#AA3939', linewidth=2)
    
    dx=0.00001
    xcoord=np.arange(0,4647999)
    xcoord_dif=np.arange(0,4647998)
    ori_position=int((3712059+3711828)/2)
    dif_position=int((1593613+1593585)/2)
        
    
    #-3 min to 0 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot()     
    plot1.plot(xcoord, wig_data_exp_SD[0][0], '--', label='0 min/-3 min +IP+Cfx', color='#AA3939', linewidth=2) #Div
    plot1.plot(xcoord, np.array(wig_data_controls_SD[0][0])/np.mean(wig_data_controls_SD[0][0][0:2500000]), ':', label='0 min/-3 min -IP+Cfx', color='#A8658B', linewidth=2) #Div
    plot1.plot(xcoord_dif, der_comp(wig_data_controls_SD[0][0], dx, 'Div'), ':', label='d(0 min/-3 min)', color='#69073E', linewidth=1) #Div    
    #plot1.plot(xcoord, wig_data_exp_SD[0][0], '--', label='0 min--3 min SD', color='#AA3939', linewidth=2) #Ded
    #plot1.plot(xcoord, (np.array(wig_data_controls_SD[0][0])+np.abs(np.mean(wig_data_controls_SD[0][0][0:2500000])))/280, ':', label='0 min--3 min -IP+Cfx', color='#A8658B', linewidth=2) #Ded
    #plot1.plot(xcoord_dif, der_comp(wig_data_controls_SD[0][0], dx, 'Ded'), ':', label='d(0 min--3 min)', color='#69073E', linewidth=1) #Ded
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([dif_position, ori_position], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    plot1.set_ylim(0.6, 1.5) #Div: 0.6, 1.5, Ded: -0.5, 1.0
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    plot1.annotate('', xytext=(ori_position, 0.8), xy=(ori_position, 0.85), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #Div: 0.8, 0.85, Ded: -0.35, -0.25
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_div_vs_seq_div\R1\\" + "R1_+IP+Cfx_seq_div_vs_-IP+Cfx_seq_div_vs_d(-IP+Cfx_seq_div)_-3_min_to_0_min_transition.png", dpi=300, figsize=(16, 8))
    plt.close()      
    
    #0 min to 5 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot()     
    plot1.plot(xcoord, wig_data_exp_SD[1][0], '--', label='5 min/0 min +IP+Cfx', color='#AA3939', linewidth=2) #Div
    plot1.plot(xcoord, np.array(wig_data_controls_SD[1][0])/np.mean(wig_data_controls_SD[1][0][0:2500000]), ':', label='5 min/0 min -IP+Cfx', color='#A8658B', linewidth=2) #Div
    plot1.plot(xcoord_dif, der_comp(wig_data_controls_SD[1][0], dx, 'Div'), ':', label='d(5 min/0 min)', color='#69073E', linewidth=1) #Div
    #plot1.plot(xcoord, wig_data_exp_SD[1][0], '--', label='5 min-0 min SD', color='#AA3939', linewidth=2) #Ded
    #plot1.plot(xcoord, (np.array(wig_data_controls_SD[1][0])+np.abs(np.mean(wig_data_controls_SD[1][0][0:2500000])))/280, ':', label='5 min-0 min -IP+Cfx', color='#A8658B', linewidth=2) #Ded
    #plot1.plot(xcoord_dif, der_comp(wig_data_controls_SD[1][0], dx, 'Ded'), ':', label='d(5 min-0 min)', color='#69073E', linewidth=1) #Ded
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([dif_position, ori_position], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    plot1.set_ylim(0.6, 1.5) #Div: 0.6, 1.5, Ded: -0.5, 1.0
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    plot1.annotate('', xytext=(ori_position, 0.8), xy=(ori_position, 0.85), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #Div: 0.8, 0.85, Ded: -0.35, -0.25
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_div_vs_seq_div\R1\\" + "R1_+IP+Cfx_seq_div_vs_-IP+Cfx_seq_div_vs_d(-IP+Cfx_seq_div)_0_min_to_5_min_transition.png", dpi=300, figsize=(16, 8))
    plt.close()         
    
    #5 min to 10 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot()     
    plot1.plot(xcoord, wig_data_exp_SD[2][0], '--', label='10 min/5 min +IP+Cfx', color='#AA3939', linewidth=2) #Div
    plot1.plot(xcoord, np.array(wig_data_controls_SD[2][0])/np.mean(wig_data_controls_SD[2][0][0:2500000]), ':', label='10 min/5 min -IP+Cfx', color='#A8658B', linewidth=2) #Div
    dcoverage=der_comp(wig_data_controls_SD[2][0], dx, 'Div') #Div
    plot1.plot(xcoord_dif, dcoverage, ':', label='d(10 min/5 min)', color='#69073E', linewidth=1) #Div    
    #plot1.plot(xcoord, wig_data_exp_SD[2][0], '--', label='10 min-5 min SD', color='#AA3939', linewidth=2) #Ded
    #plot1.plot(xcoord, (np.array(wig_data_controls_SD[2][0])+np.abs(np.mean(wig_data_controls_SD[2][0][0:2500000])))/280, ':', label='10 min-5 min -IP+Cfx', color='#A8658B', linewidth=2) #Ded
    #dcoverage=der_comp(wig_data_controls_SD[2][0], dx, 'Ded') #Ded
    #plot1.plot(xcoord_dif, dcoverage, ':', label='d(10 min-5 min)', color='#69073E', linewidth=1) #Ded
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([dif_position, ori_position], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    plot1.set_ylim(0.6, 1.5) #Div: 0.6, 1.5, Ded: -0.5, 1.0
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    plot1.annotate('', xytext=(ori_position, 0.8), xy=(ori_position, 0.85), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #Div: 0.8, 0.85, Ded: -0.35, -0.25
    #Mark local maxima in gyrase enrichment.
    max_index_gyr1=wig_data_exp_SD[2][0][3500000:4000000].index(max(wig_data_exp_SD[2][0][3500000:4000000]))+3500000
    plot1.annotate('', xytext=(max_index_gyr1, 1.4), xy=(max_index_gyr1, 1.35), arrowprops=dict(arrowstyle='->', color='#AA3939', linewidth=2), color='#AA3939', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in coverage enrichment.
    max_index_rep1=wig_data_controls_SD[2][0][3500000:4000000].index(max(wig_data_controls_SD[2][0][3500000:4000000]))+3500000
    plot1.annotate('', xytext=(max_index_rep1, 1.4), xy=(max_index_rep1, 1.35), arrowprops=dict(arrowstyle='->', color='#A8658B', linewidth=2), color='#A8658B', weight="bold", size=15)   #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in replisome enrichment.
    max_index_replisome=dcoverage[2500000:4600000].tolist().index(max(dcoverage[2500000:4600000]))+2500000
    plot1.annotate('', xytext=(max_index_replisome, 1.32), xy=(max_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35   
    #Mark local minima in replisome enrichment.
    min_index_replisome=dcoverage[2500000:4600000].tolist().index(min(dcoverage[2500000:4600000]))+2500000    
    plot1.annotate('', xytext=(min_index_replisome, 1.32), xy=(min_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35  
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_div_vs_seq_div\R1\\" + "R1_+IP+Cfx_seq_div_vs_-IP+Cfx_seq_div_vs_d(-IP+Cfx_seq_div)_10_min_to_5_min_transition.png", dpi=300, figsize=(16, 8))
    plt.close()   
    print('Distance between gyr and cov maxima during 5 to 10 min transition, bp: ' + str(np.abs(max_index_gyr1-max_index_rep1)))
    
    #10 min to 15 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot() 
    plot1.plot(xcoord, wig_data_exp_SD[3][0], '--', label='15 min/10 min +IP+Cfx', color='#AA3939', linewidth=2) #Div
    plot1.plot(xcoord, np.array(wig_data_controls_SD[3][0])/np.mean(wig_data_controls_SD[3][0][0:2500000]), ':', label='15 min/10 min -IP+Cfx', color='#A8658B', linewidth=2) #Div
    dcoverage=der_comp(wig_data_controls_SD[3][0], dx, 'Div') #Div
    plot1.plot(xcoord_dif, dcoverage, ':', label='d(15 min/10 min)', color='#69073E', linewidth=1) #Div   
    #plot1.plot(xcoord, wig_data_exp_SD[3][0], '--', label='15 min-10 min SD', color='#AA3939', linewidth=2) #Ded
    #plot1.plot(xcoord, (np.array(wig_data_controls_SD[3][0])+np.abs(np.mean(wig_data_controls_SD[3][0][0:2500000])))/280, ':', label='15 min-10 min -IP+Cfx', color='#A8658B', linewidth=2) #Ded
    #dcoverage=der_comp(wig_data_controls_SD[3][0], dx, 'Ded') #Ded
    #plot1.plot(xcoord_dif, dcoverage, ':', label='d(15 min-10 min)', color='#69073E', linewidth=1) #Ded
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([dif_position, ori_position], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    plot1.set_ylim(0.6, 1.5) #Div: 0.6, 1.5, Ded: -0.5, 1.0
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    plot1.annotate('', xytext=(ori_position, 0.8), xy=(ori_position, 0.85), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #Div: 0.8, 0.85, Ded: -0.35, -0.25
    #Mark local maxima in gyrase enrichment.
    max_index_gyr1=wig_data_exp_SD[3][0][3500000:4000000].index(max(wig_data_exp_SD[3][0][3500000:4000000]))+3500000
    plot1.annotate('', xytext=(max_index_gyr1, 1.4), xy=(max_index_gyr1, 1.35), arrowprops=dict(arrowstyle='->', color='#AA3939', linewidth=2), color='#AA3939', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in coverage enrichment.
    max_index_rep1=wig_data_controls_SD[3][0][3500000:4000000].index(max(wig_data_controls_SD[3][0][3500000:4000000]))+3500000
    plot1.annotate('', xytext=(max_index_rep1, 1.4), xy=(max_index_rep1, 1.35), arrowprops=dict(arrowstyle='->', color='#A8658B', linewidth=2), color='#A8658B', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in replisome enrichment.
    max_index_replisome=dcoverage[2500000:4600000].tolist().index(max(dcoverage[2500000:4600000]))+2500000
    plot1.annotate('', xytext=(max_index_replisome, 1.32), xy=(max_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35   
    #Mark local minima in replisome enrichment.
    min_index_replisome=dcoverage[2500000:4600000].tolist().index(min(dcoverage[2500000:4600000]))+2500000    
    plot1.annotate('', xytext=(min_index_replisome, 1.32), xy=(min_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35       
    #Title axis.
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_div_vs_seq_div\R1\\" + "R1_+IP+Cfx_seq_div_vs_-IP+Cfx_seq_div_vs_d(-IP+Cfx_seq_div)_15_min_to_10_min_transition.png", dpi=300, figsize=(16, 8))
    plt.close()  
    print('Distance between gyr and cov maxima during 10 to 15 min transition, bp: ' + str(np.abs(max_index_gyr1-max_index_rep1)))
    
    #15 min to 20 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot() 
    plot1.plot(xcoord, wig_data_exp_SD[4][0], '--', label='20 min/15 min +IP+Cfx', color='#AA3939', linewidth=2) #Div
    plot1.plot(xcoord, np.array(wig_data_controls_SD[4][0])/np.mean(wig_data_controls_SD[4][0][0:2500000]), ':', label='20 min/15 min -IP+Cfx', color='#A8658B', linewidth=2) #Div
    dcoverage=der_comp(wig_data_controls_SD[4][0], dx, 'Div') #Div
    plot1.plot(xcoord_dif, dcoverage, ':', label='d(20 min/15 min)', color='#69073E', linewidth=1)    #Div
    #plot1.plot(xcoord, wig_data_exp_SD[4][0], '--', label='20 min-15 min SD', color='#AA3939', linewidth=2) #Ded
    #plot1.plot(xcoord, (np.array(wig_data_controls_SD[4][0])+np.abs(np.mean(wig_data_controls_SD[4][0][0:2500000])))/280, ':', label='20 min-15 min -IP+Cfx', color='#A8658B', linewidth=2) #Ded
    #dcoverage=der_comp(wig_data_controls_SD[4][0], dx, 'Ded')
    #plot1.plot(xcoord_dif, dcoverage, ':', label='d(20 min-15 min)', color='#69073E', linewidth=1)
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([dif_position, ori_position], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    plot1.set_ylim(0.6, 1.5) #Div: 0.6, 1.5, Ded: -0.5, 1.0
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    plot1.annotate('', xytext=(ori_position, 0.8), xy=(ori_position, 0.85), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15)   #Div: 0.8, 0.85, Ded: -0.35, -0.25
    #Mark local maxima in gyrase enrichment.
    max_index_gyr1=wig_data_exp_SD[4][0][3200000:ori_position].index(max(wig_data_exp_SD[4][0][3200000:ori_position]))+3200000
    plot1.annotate('', xytext=(max_index_gyr1, 1.4), xy=(max_index_gyr1, 1.35), arrowprops=dict(arrowstyle='->', color='#AA3939', linewidth=2), color='#AA3939', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    max_index_gyr2=wig_data_exp_SD[4][0][ori_position:4500000].index(max(wig_data_exp_SD[4][0][ori_position:4500000]))+ori_position
    plot1.annotate('', xytext=(max_index_gyr2, 1.4), xy=(max_index_gyr2, 1.35), arrowprops=dict(arrowstyle='->', color='#AA3939', linewidth=2), color='#AA3939', weight="bold", size=15)  #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in coverage enrichment.
    max_index_rep1=wig_data_controls_SD[4][0][3200000:ori_position].index(max(wig_data_controls_SD[4][0][3200000:ori_position]))+3200000
    plot1.annotate('', xytext=(max_index_rep1, 1.4), xy=(max_index_rep1, 1.35), arrowprops=dict(arrowstyle='->', color='#A8658B', linewidth=2), color='#A8658B', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    max_index_rep2=wig_data_controls_SD[4][0][ori_position:4500000].index(max(wig_data_controls_SD[4][0][ori_position:4500000]))+ori_position
    plot1.annotate('', xytext=(max_index_rep2, 1.4), xy=(max_index_rep2, 1.35), arrowprops=dict(arrowstyle='->', color='#A8658B', linewidth=2), color='#A8658B', weight="bold", size=15)    #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in replisome enrichment.
    max_index_replisome=dcoverage[2500000:4600000].tolist().index(max(dcoverage[2500000:4600000]))+2500000
    plot1.annotate('', xytext=(max_index_replisome, 1.32), xy=(max_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35  
    #Mark local minima in replisome enrichment.
    min_index_replisome=dcoverage[2500000:4600000].tolist().index(min(dcoverage[2500000:4600000]))+2500000    
    plot1.annotate('', xytext=(min_index_replisome, 1.32), xy=(min_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35       
    #Title axis.    
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_div_vs_seq_div\R1\\" + "R1_+IP+Cfx_seq_div_vs_-IP+Cfx_seq_div_vs_d(-IP+Cfx_seq_div)_20_min_to_15_min_transition.png", dpi=300, figsize=(16, 8))
    plt.close()     
    print('Distance between gyr and cov maxima on the left replichore during 15 to 20 min transition, bp: ' + str(np.abs(max_index_gyr1-max_index_rep1)))
    print('Distance between gyr and cov maxima on the right replichore during 15 to 20 min transition, bp: ' + str(np.abs(max_index_gyr2-max_index_rep2)))
    print('Distance between gyr and rep maxima on the left replichore during 15 to 20 min transition, bp: ' + str(max_index_replisome-max_index_gyr1))
    print('Distance between gyr and rep maxima on the right replichore during 15 to 20 min transition, bp: ' + str(max_index_gyr2-min_index_replisome))
    
    #20 min to 25 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot() 
    plot1.plot(xcoord, wig_data_exp_SD[5][0], '--', label='25 min/20 min +IP+Cfx', color='#AA3939', linewidth=2) #Div
    plot1.plot(xcoord, np.array(wig_data_controls_SD[5][0])/np.mean(wig_data_controls_SD[5][0][0:2500000]), ':', label='25 min/20 min -IP+Cfx', color='#A8658B', linewidth=2) #Div
    dcoverage=der_comp(wig_data_controls_SD[5][0], dx, 'Div') #Div
    plot1.plot(xcoord_dif, dcoverage, ':', label='d(25 min/20 min)', color='#69073E', linewidth=1)    #Div
    #plot1.plot(xcoord, wig_data_exp_SD[5][0], '--', label='25 min-20 min SD', color='#AA3939', linewidth=2) #Ded
    #plot1.plot(xcoord, (np.array(wig_data_controls_SD[5][0])+np.abs(np.mean(wig_data_controls_SD[5][0][0:2500000])))/280, ':', label='25 min-20 min -IP+Cfx', color='#A8658B', linewidth=2) #Ded
    #dcoverage=der_comp(wig_data_controls_SD[5][0], dx, 'Ded')
    #plot1.plot(xcoord_dif, dcoverage, ':', label='d(25 min-20 min)', color='#69073E', linewidth=1)
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([dif_position, ori_position], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    plot1.set_ylim(0.6, 1.5) #Div: 0.6, 1.5, Ded: -0.5, 1.0
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    plot1.annotate('', xytext=(ori_position, 0.8), xy=(ori_position, 0.85), arrowprops=dict(arrowstyle='->', color='green'), color='black', weight="bold", size=15) #Div: 0.8, 0.85, Ded: -0.35, -0.25
    #Mark local maxima in gyrase enrichment.
    max_index_gyr1=wig_data_exp_SD[5][0][3000000:ori_position].index(max(wig_data_exp_SD[5][0][3000000:ori_position]))+3000000
    plot1.annotate('', xytext=(max_index_gyr1, 1.4), xy=(max_index_gyr1, 1.35), arrowprops=dict(arrowstyle='->', color='#AA3939', linewidth=2), color='#AA3939', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    max_index_gyr2=wig_data_exp_SD[5][0][ori_position:4500000].index(max(wig_data_exp_SD[5][0][ori_position:4500000]))+ori_position
    plot1.annotate('', xytext=(max_index_gyr2, 1.4), xy=(max_index_gyr2, 1.35), arrowprops=dict(arrowstyle='->', color='#AA3939', linewidth=2), color='#AA3939', weight="bold", size=15)  #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in coverage enrichment.
    max_index_rep1=wig_data_controls_SD[5][0][3200000:ori_position].index(max(wig_data_controls_SD[5][0][3200000:ori_position]))+3200000
    plot1.annotate('', xytext=(max_index_rep1, 1.4), xy=(max_index_rep1, 1.35), arrowprops=dict(arrowstyle='->', color='#A8658B', linewidth=2), color='#A8658B', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    max_index_rep2=wig_data_controls_SD[5][0][ori_position:4500000].index(max(wig_data_controls_SD[5][0][ori_position:4500000]))+ori_position
    plot1.annotate('', xytext=(max_index_rep2, 1.4), xy=(max_index_rep2, 1.35), arrowprops=dict(arrowstyle='->', color='#A8658B', linewidth=2), color='#A8658B', weight="bold", size=15)      #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in replisome enrichment.
    max_index_replisome=dcoverage[2500000:4600000].tolist().index(max(dcoverage[2500000:4600000]))+2500000
    plot1.annotate('', xytext=(max_index_replisome, 1.32), xy=(max_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35  
    #Mark local minima in replisome enrichment.
    min_index_replisome=dcoverage[2500000:4600000].tolist().index(min(dcoverage[2500000:4600000]))+2500000    
    plot1.annotate('', xytext=(min_index_replisome, 1.32), xy=(min_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35        
    #Title axis.    
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_div_vs_seq_div\R1\\" + "R1_+IP+Cfx_seq_div_vs_-IP+Cfx_seq_div_vs_d(-IP+Cfx_seq_div)_25_min_to_20_min_transition.png", dpi=300, figsize=(16, 8))
    plt.close() 
    print('Distance between gyr and cov maxima on the left replichore during 20 to 25 min transition, bp: ' + str(np.abs(max_index_gyr1-max_index_rep1)))
    print('Distance between gyr and cov maxima on the right replichore during 20 to 25 min transition, bp: ' + str(np.abs(max_index_gyr2-max_index_rep2)))  
    print('Distance between gyr and rep maxima on the left replichore during 20 to 25 min transition, bp: ' + str(max_index_replisome-max_index_gyr1))
    print('Distance between gyr and rep maxima on the right replichore during 20 to 25 min transition, bp: ' + str(max_index_gyr2-min_index_replisome))
    
    #25 min to 30 min transition.
    print("Samples are plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("E. coli gyrase enrichment and replisome positioning", fontsize=20)
    plot1=plt.subplot() 
    plot1.plot(xcoord, wig_data_exp_SD[6][0], '--', label='30 min/25 min +IP+Cfx', color='#AA3939', linewidth=2) #Div
    plot1.plot(xcoord, np.array(wig_data_controls_SD[6][0])/np.mean(wig_data_controls_SD[6][0][0:2500000]), ':', label='30 min/25 min -IP+Cfx', color='#A8658B', linewidth=2) #Div
    dcoverage=der_comp(wig_data_controls_SD[6][0], dx, 'Div') #Div
    plot1.plot(xcoord_dif, dcoverage, ':', label='d(30 min/25 min)', color='#69073E', linewidth=1) #Div   
    #plot1.plot(xcoord, wig_data_exp_SD[6][0], '--', label='30 min-25 min SD', color='#AA3939', linewidth=2) #Ded
    #plot1.plot(xcoord, (np.array(wig_data_controls_SD[6][0])+np.abs(np.mean(wig_data_controls_SD[6][0][0:2500000])))/280, ':', label='30 min-25 min -IP+Cfx', color='#A8658B', linewidth=2) #Ded
    #dcoverage=der_comp(wig_data_controls_SD[6][0], dx, 'Ded') #Ded
    #plot1.plot(xcoord_dif, dcoverage, ':', label='d(30 min-25 min)', color='#69073E', linewidth=1) #Ded
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([dif_position, ori_position], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    plot1.set_ylim(0.6, 1.5) #Div: 0.6, 1.5, Ded: -0.5, 1.0
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    plot1.annotate('', xytext=(ori_position, 0.8), xy=(ori_position, 0.85), arrowprops=dict(arrowstyle='->', color='green'), color='green', weight="bold", size=15) #Div: 0.8, 0.85, Ded: -0.35, -0.25
    #Mark local maxima in gyrase enrichment.
    max_index_gyr1=wig_data_exp_SD[6][0][2700000:ori_position].index(max(wig_data_exp_SD[6][0][2700000:ori_position]))+2700000
    plot1.annotate('', xytext=(max_index_gyr1, 1.4), xy=(max_index_gyr1, 1.35), arrowprops=dict(arrowstyle='->', color='#AA3939', linewidth=2), color='#AA3939', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    max_index_gyr2=wig_data_exp_SD[6][0][ori_position:4600000].index(max(wig_data_exp_SD[6][0][ori_position:4600000]))+ori_position
    plot1.annotate('', xytext=(max_index_gyr2, 1.4), xy=(max_index_gyr2, 1.35), arrowprops=dict(arrowstyle='->', color='#AA3939',linewidth=2), color='#AA3939', weight="bold", size=15)  #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    #Mark local maxima in coverage enrichment.
    max_index_rep1=wig_data_controls_SD[6][0][3000000:ori_position].index(max(wig_data_controls_SD[6][0][3000000:ori_position]))+3000000
    plot1.annotate('', xytext=(max_index_rep1, 1.4), xy=(max_index_rep1, 1.35), arrowprops=dict(arrowstyle='->', color='#A8658B', linewidth=2), color='#A8658B', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5 
    max_index_rep2=wig_data_controls_SD[6][0][ori_position:4600000].index(max(wig_data_controls_SD[6][0][ori_position:4600000]))+ori_position
    plot1.annotate('', xytext=(max_index_rep2, 1.4), xy=(max_index_rep2, 1.35), arrowprops=dict(arrowstyle='->', color='#A8658B', linewidth=2), color='#A8658B', weight="bold", size=15) #Div: 1.4, 1.35; Ded: 0.6, 0.5    
    #Mark local maxima in replisome enrichment.
    max_index_replisome=dcoverage[2500000:4600000].tolist().index(max(dcoverage[2500000:4600000]))+2500000
    plot1.annotate('', xytext=(max_index_replisome, 1.32), xy=(max_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35  
    #Mark local minima in replisome enrichment.
    min_index_replisome=dcoverage[2500000:4600000].tolist().index(min(dcoverage[2500000:4600000]))+2500000    
    plot1.annotate('', xytext=(min_index_replisome, 1.32), xy=(min_index_replisome, 1.27), arrowprops=dict(arrowstyle='->', color='#69073E', linewidth=2), color='#69073E', weight="bold", size=15)   #Div: 1.32, 1.27; Ded: 0.45, 0.35    
    #Title axis.      
    plot1.set_xlabel('Genome position, kb', size=20)
    plot1.set_ylabel('Coverage depth normalized', size=20)
    plot1.legend(loc='upper left', fontsize=25)
    plt.show()
    plt.savefig(path_out + "Figures\Seq_div_vs_seq_div\R1\\" + "R1_+IP+Cfx_seq_div_vs_-IP+Cfx_seq_div_vs_d(-IP+Cfx_seq_div)_30_min_to_25_min_transition.png", dpi=300, figsize=(16, 8))
    plt.close()  
    print('Distance between gyr and cov maxima on the left replichore during 25 to 30 min transition, bp: ' + str(np.abs(max_index_gyr1-max_index_rep1)))
    print('Distance between gyr and cov maxima on the right replichore during 25 to 30 min transition, bp: ' + str(np.abs(max_index_gyr2-max_index_rep2)))  
    print('Distance between gyr and rep maxima on the left replichore during 25 to 30 min transition, bp: ' + str(max_index_replisome-max_index_gyr1))
    print('Distance between gyr and rep maxima on the right replichore during 25 to 30 min transition, bp: ' + str(max_index_gyr2-min_index_replisome))
    return

#Repl_Gyr_ident(WIG_cont_NSSDF_input, WIG_exp_NSSDF_input, 20, Deletions, "NC_007779.1_w3110_Mu", Output_directory_div)

#######
#Calculate and plot replication rate.
#######

def calc_repl_rate():
    time=np.array([0, 5, 10, 15])
    distance_left_R1=np.array([85274.5,225483.5,412856.5,690296.5])
    LR1_ap=np.polyfit(time,distance_left_R1,1)
    distance_right_R1=np.array([89843.5,227152.5,471353.5,669413.5])
    RR1_ap=np.polyfit(time,distance_right_R1,1)
    distance_left_R2=np.array([80333.5,233104.5,468141.5,699638.5])
    LR2_ap=np.polyfit(time,distance_left_R2,1)
    distance_right_R2=np.array([102829.5,253880.5,464516.5,775229.5])
    RR2_ap=np.polyfit(time,distance_right_R2,1)
    Rates_ar=[LR1_ap[0], RR1_ap[0], LR2_ap[0], RR2_ap[0]]
    Lag_ar=[LR1_ap[1], RR1_ap[1], LR2_ap[1], RR2_ap[1]]
    Mean_Rep_Rate=np.mean(Rates_ar)
    SD_Rep_Rate=np.std(Rates_ar)
    Mean_Lag=np.mean(Lag_ar)
    SD_Lag=np.std(Lag_ar)
    bottom_edge=(time*(Mean_Rep_Rate-SD_Rep_Rate))+Mean_Lag-SD_Lag
    upper_edge=(time*(Mean_Rep_Rate+SD_Rep_Rate))+Mean_Lag+SD_Lag
    print(Mean_Rep_Rate, SD_Rep_Rate)
    plt.figure(figsize=(10, 7), dpi=100)
    plt.suptitle("E. coli replisome rate", fontsize=20)
    plot1=plt.subplot() 
    plot1.fill_between(time, bottom_edge, upper_edge, facecolor='green', alpha=0.3, zorder=10)
    plot1.plot(time, distance_left_R1, '-', label='R1 left replichore', color='#D17E7E', linewidth=2)
    plot1.plot(time, bottom_edge, '-', color='black', linewidth=0.2)
    plot1.plot(time, upper_edge, '-', color='black', linewidth=0.2)
    plot1.plot(time, distance_right_R1, '-', label='R1 right replichore', color='#882D60', linewidth=2)
    plot1.plot(time, distance_left_R2, '--', label='R2 left replichore', color='#D17E7E', linewidth=2)
    plot1.plot(time, distance_right_R2, '--', label='R2 right replichore', color='#882D60', linewidth=2)
    
    plot1.annotate(r'$Replisome\ rate=$'+str(np.round(Mean_Rep_Rate/1000,1))+r'$kb\pm$'+str(np.round(SD_Rep_Rate/1000,1))+r'$kb$', xytext=(7.5, 200000), xy=(7.5, 200000), color='black', weight="bold", size=15) #Div: 0.8, 0.85, Ded: -0.35, -0.25
    
    plot1.set_yticks(np.arange(0,800001,200000), minor=False)
    plot1.set_yticklabels(np.arange(0,801,200), minor=False)
    
    plot1.set_xticks(np.arange(0,16,5), minor=False)
    plot1.set_xticklabels(np.arange(15,31,5), minor=False)    

    plot1.tick_params(axis='both', which='major', labelsize=17)  
    plot1.set_xlabel('Time from synchronization start, min', size=20)
    plot1.set_ylabel('Distance from Ori, kb', size=20)
    plot1.legend(loc='upper left', fontsize=20)
    plt.show()
    #plt.savefig(path_out  + prefix + str(i+1) + '_' + str(i) + suffix, dpi=300, figsize=(16, 8)) #, format='eps'
    #plt.close()    
    return

#calc_repl_rate()


#######
#Compare smoothing by sliding window of different width with Fourier filtration.
#######

def Window_vs_Fourier(wig_in_raw, wig_in_NS, wig_in_ff, del_path, path_out):
    #Read deletions info.
    deletions=deletions_info(del_path)
    #Read data raw.
    wig_data_RAW=wig_parsing(wig_in_raw) 
    #Read data NS.
    wig_data_NS=wig_parsing(wig_in_NS) 
    #Read data FF.
    wig_data_FF=wig_parsing(wig_in_ff)        
    
    #Smooth normalized coverage.
    print("Samples are smoothing...")
    windows=[5000, 10000, 50000, 100000]
    wig_data_NS_sm=[]
    i=0
    for win_width in windows:
        print("Now sample " + str(i+1) + " out of " + str(len(windows)) + " is smoothing...")
        sample_sm=Smoothing(wig_data_NS, deletions, win_width)
        wig_data_NS_sm.append(sample_sm)
        i+=1   
    wig_data_RAW_sm=[]
    i=0
    for win_width in windows:
        print("Now sample " + str(i+1) + " out of " + str(len(windows)) + " is smoothing...")
        sample_sm=Smoothing(wig_data_RAW, deletions, win_width)
        wig_data_RAW_sm.append(sample_sm)
        i+=1  
    
    #Plot raw data, NS data, FF data.
    X=np.arange(0,len(wig_data_NS))
    ori_position=int((3712059+3711828)/2)
    dif_position=int((1593613+1593585)/2)     
    #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
    #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
    print("Sample is plotting...")
    mpl.rcParams['agg.path.chunksize']=10000
    plt.figure(figsize=(16, 8), dpi=100)
    plt.suptitle("Methods for coverage smoothing comparision", fontsize=20)
    plot1=plt.subplot()    
    
    plot1.plot(X, wig_data_RAW, '--', label='Raw coverage', color='#550000', linewidth=0.2)
    plot1.plot(X, wig_data_NS, '--', label='Coverage NS W1kb', color='#801515', linewidth=0.5)
    plot1.plot(X, wig_data_FF, '--', label='Coverge NS W1kb FF20', color='#D17E7E', linewidth=1)
    plot1.plot(X, wig_data_NS_sm[0], '--', label='Coverage NS W5kb', color='#490029', linewidth=1)
    plot1.plot(X, wig_data_NS_sm[1], '--', label='Coverage NS W10kb', color='#69073E', linewidth=1)
    plot1.plot(X, wig_data_NS_sm[2], '--', label='Coverage NS W50kb', color='#882D60', linewidth=1)
    plot1.plot(X, wig_data_NS_sm[3], '--', label='Coverage NS W100kb', color='#A8658B', linewidth=1)    
    plot1.plot(X, wig_data_RAW_sm[0], '--', label='Coverage RAW W5kb', color='#490029', linewidth=1)
    plot1.plot(X, wig_data_RAW_sm[1], '--', label='Coverage RAW W10kb', color='#69073E', linewidth=1)
    plot1.plot(X, wig_data_RAW_sm[2], '--', label='Coverage RAW W50kb', color='#882D60', linewidth=1)
    plot1.plot(X, wig_data_RAW_sm[3], '--', label='Coverage RAW W100kb', color='#A8658B', linewidth=1)    
    
    plot1.set_xticks([0,500000,1000000,1500000,2000000,2500000,3000000,3500000, 4000000,4500000], minor=False)
    plot1.set_xticklabels([0, '500', '1000', '1500', '2000', '2500', '3000', '3500', '4000', '4500'], minor=False)
    plot1.set_xticks([dif_position, ori_position], minor=True)
    plot1.set_xticklabels(['dif', 'Ori'], va='top', minor=True)
    plot1.tick_params(axis='x', which='minor', direction='in', pad=-22, labelsize=17)
    #plot1.set_ylim(0.6, 2)
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    #Mark origin
    #plot1.annotate('', xytext=(ori_position, 0.7), xy=(ori_position, 0.8), arrowprops=dict(arrowstyle='->', facecolor='green'), color='black', weight="bold", size=15)
    plot1.set_xlabel('Genome position, nt', size=17)
    plot1.set_ylabel('Coverage depth smoothed', size=17)
    plot1.legend(loc='upper left')
    plt.show()
    plt.savefig(path_out  + "Figures\SW_vs_FF\\" + "-IP+Cfx_R1_30min_RAW_RAWW_NS_NSW_FF_copm.png", dpi=300, figsize=(16, 8))
    plt.close()     
    return

#Window_vs_Fourier(WIG_input_NN+'Raw_data\-IP+Cfx_R1\\7_DSu_16_S96_edt_for_rev_depth.wig', WIG_input_NN+'NS\-IP+Cfx_R1\\7_-IP+Cfx_R1_normalized_smoothed_1000bp.wig',
#                  WIG_input_NN+'NSF\-IP+Cfx_R1\\7_-IP+Cfx_R1_normalized_smoothed_1000bp_filtered_20fh.wig', Deletions, Output_directory)

#######
#Calculates ratio between neutral region coverage and coverage at Ori.
#######

def calc_ratio(wig_ar):
    ori_position=int((3712059+3711828)/2)
    cov_neutral=np.mean(wig_ar[0:2000000])
    cov_ori=wig_ar[ori_position]
    rep_ratio=cov_ori/cov_neutral
    return rep_ratio


#######
#Calculate fraction of cells undergo replication.
#######

def Frac_of_rep(wig_in_raw_dir_R1, wig_in_raw_dir_R2, del_path, path_out):
    '''
    #Read deletions info.
    deletions=deletions_info(del_path)
    #Read data raw R1.
    wig_data_RAW_R1=read_wig_files(wig_in_raw_dir_R1)   
    #Read data raw R2.
    wig_data_RAW_R2=read_wig_files(wig_in_raw_dir_R2)     
    
    #Smooth normalized coverage.
    print("Samples are smoothing...")
    windows=[1000, 5000, 10000]
    wig_data_RAW_sm_R1=[]
    wig_data_RAW_sm_R2=[]
    for win_width in windows:
        some_window_width_R1=[]
        some_window_width_R2=[]
        for i in range(len(wig_data_RAW_R1)):
            print("Now sample " + str(i+1) + " out of " + str(len(wig_data_RAW_R1)) + " is smoothing...")
            sample_sm_R1=Smoothing(wig_data_RAW_R1[i], deletions, win_width)
            some_window_width_R1.append(sample_sm_R1)  
            sample_sm_R2=Smoothing(wig_data_RAW_R2[i], deletions, win_width)
            some_window_width_R2.append(sample_sm_R2)              
        wig_data_RAW_sm_R1.append(some_window_width_R1)
        wig_data_RAW_sm_R2.append(some_window_width_R2)
    
    #Calculate Ori to Neutral region coverage ratio.
    #NS smoothed additionally ratios.
    NS_sm_ratios_R1=[]
    NS_sm_ratios_R2=[]
    for i in range(len(wig_data_RAW_sm_R1)):
        win_width_ratios_R1=[]
        win_width_ratios_R2=[]
        for j in range(len(wig_data_RAW_sm_R1[i])):
            time_point_ratio_R1=calc_ratio(wig_data_RAW_sm_R1[i][j])
            win_width_ratios_R1.append(time_point_ratio_R1)
            time_point_ratio_R2=calc_ratio(wig_data_RAW_sm_R2[i][j])
            win_width_ratios_R2.append(time_point_ratio_R2)            
        NS_sm_ratios_R1.append(win_width_ratios_R1)
        NS_sm_ratios_R2.append(win_width_ratios_R2)
    '''
    
    #Sample is plotting...
    NS_sm_ratios_R1=[[1.0810971233078916, 1.0116723339859484, 1.0746432048545937, 1.1593970526147852, 1.313672815046678, 1.6453463391390668, 1.87642899450788, 2.0094741195423764],
                     [1.1168755139816495, 1.0856453615856183, 1.0908837338378339, 1.1771383965652389, 1.3285875952833743, 1.6271142527630722, 1.8680633313596464, 2.007372358057889],
                     [1.1128374238875716, 1.0899834883812058, 1.0796765427141313, 1.1565044968741585, 1.3322937163670774, 1.630821398281644, 1.8510050702451664, 1.9882109446575125]]
    NS_sm_ratios_R2=[[1.1583278715481657, 1.104359974303709, 1.0956745067856233, 1.2628976677706059, 1.3858047454370084, 1.5857462985250037, 1.84994207114619, 2.022672577752544],
                     [1.1628091237712161, 1.115230195789552, 1.1013557631635666, 1.2218642720776254, 1.3906649948976513, 1.6345465113318387, 1.8752097232479301, 2.038586178173467],
                     [1.1584063555918824, 1.1220269490727266, 1.1025915160933062, 1.2018507136460288, 1.3774946096489158, 1.6221849894704616, 1.841048088632014, 2.006613169697742]]
    #Zero point position: 1.0908042112415093
    #[0.9, 1.0908042112415093, 1.1908042112415094, 1.2908042112415092, 1.3908042112415093, 1.4908042112415094, 1.5908042112415093, 1.6908042112415094, 1.7908042112415092, 1.8908042112415093, 1.9908042112415094, 2.090804211241509]
    #[-10   0  10  20  30  40  50  60  70  80  90 100]
    
    #Plot fraction of cells undergo replication.
    X=np.arange(0,8)  
    print("Sample is plotting...")
    plt.figure(figsize=(10, 7), dpi=100)
    plt.suptitle("Fraction of cells undergo replication", fontsize=20)
    plot1=plt.subplot()
    plot1.plot(X, NS_sm_ratios_R1[0], ':', label='Raw W1kb R1', color='black', linewidth=1)
    plot1.plot(X, NS_sm_ratios_R1[1], ':', label='RAW W5kb R1', color='#69073E', linewidth=1)
    plot1.plot(X, NS_sm_ratios_R1[2], ':', label='RAW W10kb R1', color='#D17E7E', linewidth=1)
    plot1.plot(X, NS_sm_ratios_R2[0], '-', label='RAW W1kb R2', color='black', linewidth=1)
    plot1.plot(X, NS_sm_ratios_R2[1], '-', label='RAW W5kb R2', color='#69073E', linewidth=1)
    plot1.plot(X, NS_sm_ratios_R2[2], '-', label='RAW W10kb R2', color='#D17E7E', linewidth=1)  
    
    plot1.set_xticks(X, minor=False)
    plot1.set_xticklabels(['-3', '0', '5', '10', '15', '20', '25', '30'], minor=False)
    zero_ypoint=np.mean([NS_sm_ratios_R1[0][2], NS_sm_ratios_R1[1][2], NS_sm_ratios_R1[2][2], NS_sm_ratios_R2[0][2], NS_sm_ratios_R2[1][2], NS_sm_ratios_R2[2][2]])
    print('Zero point position: ' + str(zero_ypoint))
    yticks_coords=[]
    yticks_coords.append(zero_ypoint -0.1)
    for i in range(11):
        yticks_coords.append(zero_ypoint + i*0.1)
    plot1.set_yticks(yticks_coords, minor=False)
    plot1.set_yticklabels(np.arange(-10,101,10), minor=False)
    print(yticks_coords)
    print(np.arange(-10,101,10))
    plot1.set_yticks([zero_ypoint, zero_ypoint+0.95], minor=True)
    plot1.yaxis.grid(True, which='minor', linewidth=0.4, linestyle='--', color='black')    
    plot1.tick_params(axis='both', which='major', direction='out', labelsize=17)
    plot1.set_xlabel('Time-points, min', size=17)
    plot1.set_ylabel('Fraction of replicating cells, %', size=17)
    plot1.legend(loc='upper left', fontsize=20)
    plt.show()
    plt.savefig(path_out  + "Figures\Fraction_of_replicated\\" + "Frac_of_rep_R1_R2_adj.png", dpi=300, figsize=(16, 8))
    plt.close()     
    return

#Frac_of_rep(WIG_input_NN+'Raw_data\-IP+Cfx_R1\\', WIG_input_NN+'Raw_data\-IP+Cfx_R2\\', Deletions, Output_directory)



    
#Primarily for colour codes.
#Plot normalized and smoothed coverage.
#Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
#See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
#print("Samples are plotting...")
#mpl.rcParams['agg.path.chunksize']=10000
#xcoord=np.arange(0,4647999)
#plt.figure(figsize=(16, 8), dpi=100)
#plt.suptitle("E. coli synchronization", fontsize=20)
#plot1=plt.subplot() 
#plot1.plot(xcoord, wig_data_norm_sm[0], '--', label='-3 min', color='#550000', linewidth=1)
#plot1.plot(xcoord, wig_data_norm_sm[1], '--', label='0 min', color='#801515', linewidth=1)
#plot1.plot(xcoord, wig_data_norm_sm[2], '--', label='5 min', color='#AA3939', linewidth=1)
#plot1.plot(xcoord, wig_data_norm_sm[3], '--', label='10 min', color='#D17E7E', linewidth=1)
#plot1.plot(xcoord, wig_data_norm_sm[4], '--', label='15 min', color='#490029', linewidth=1)
#plot1.plot(xcoord, wig_data_norm_sm[5], '--', label='20 min', color='#69073E', linewidth=1)
#plot1.plot(xcoord, wig_data_norm_sm[6], '--', label='25 min', color='#882D60', linewidth=1)
#plot1.plot(xcoord, wig_data_norm_sm[7], '--', label='30 min', color='#A8658B', linewidth=1)  
#plot1.set_xlabel('Genome position, nt', size=17)
#plot1.set_ylabel('Coverage depth', size=17)
#plot1.legend(loc='upper left')
#plt.show()
#plt.savefig(path_out + 'Replic_2_+IP+Cfx_coverage_dynamic_1000.png', dpi=300, figsize=(16, 8))
#plt.close()  

'''
#######
#Normalize with neutral regioin coverage (isn't influenced by replication).
#######

def norm_neutral(wig_dir, n_round, strain_id, path_out, prefix, suffix):
    #Read data.
    wig_data=read_wig_files(wig_dir)
    #Identify minimal coverage.
    mean_coverage_ar=[]
    for sample in wig_data:
        mean_coverage_ar.append(np.mean(sample[0:2500000]))
    Min_mean_coverage=min(mean_coverage_ar)
    print('Min_mean_coverage: ' + str(Min_mean_coverage))    
    #Normalize samples coverage (add pseudocount +1).
    print("Samples are normalizing...")
    wig_data_norm=[]
    i=0
    for sample in wig_data:
        print("Now sample " + str(i+1) + " out of " + str(len(wig_data)) + " is normalizing...")
        sample_norm=[x * Min_mean_coverage/mean_coverage_ar[i] for x in sample]
        wig_data_norm.append(sample_norm)
        write_wig(sample_norm, path_out + prefix + str(i) + suffix, i, strain_id, n_round) #filtered_20fh_
        i+=1    
    return

#norm_neutral(WIG_input_NN+"\\NS\-IP+Cfx_R1\\", 15, "NC_007779.1_w3110_Mu", Output_directory, 'WIG_files\\NS_NN\-IP+Cfx\R1\\', "_-IP+Cfx_R1_normalized_smoothed_1000bp_nn.wig")
#norm_neutral(WIG_input_NN+"\\NSF\-IP+Cfx_R1\\", 15, "NC_007779.1_w3110_Mu", Output_directory, 'WIG_files\\NSF_NN\-IP+Cfx\R1\\', "_-IP+Cfx_R1_normalized_smoothed_1000bp_filtered_20fh_nn.wig")
'''