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

#Folder with WIG files.
WIG_input="C:\Sutor\science\DNA-gyrase\Results\E_coli_synch_time-course_Topo-Seq\Seq_results\WIG.tar\WIG\+IP+Cfx_R2\\"
#Deletions data.
Deletions="C:\Sutor\science\DNA-gyrase\scripts\Gyrase_Topo-seq\Additional_genome_features\Deletions_w3110_G_Mu_SGS.broadPeak"
#Output folder for NSF data.
Output_directory="C:\Sutor\science\DNA-gyrase\Results\E_coli_synch_time-course_Topo-Seq\Seq_results\Coverage_analysis\\"

#Folder with WIG files contains normalized, smoothed and filtered -IP+Cfx data.
WIG_control_NSF_input="C:\Sutor\science\DNA-gyrase\Results\E_coli_synch_time-course_Topo-Seq\Seq_results\Coverage_analysis\WIG_files\\NSF\-IP+Cfx_R2\\"
#Folder with WIG files contains normalized, smoothed +IP+Cfx data.
WIG_exp_NS_input="C:\Sutor\science\DNA-gyrase\Results\E_coli_synch_time-course_Topo-Seq\Seq_results\Coverage_analysis\WIG_files\\NS\+IP+Cfx_R2\\"

#Folder with WIG files contains normalized, smoothed, sequentially divided and filtered +IP+Cfx data.
WIG_exp_NSSDF_input="C:\Sutor\science\DNA-gyrase\Results\E_coli_synch_time-course_Topo-Seq\Seq_results\Coverage_analysis\WIG_files\Seq_ded\Seq_ded_filt\R2\\"

#Output folder for Div data.
Output_directory_div="C:\Sutor\science\DNA-gyrase\Results\E_coli_synch_time-course_Topo-Seq\Seq_results\Coverage_analysis\\"

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
    Total_C=0
    for line in wigin:
        line=line.rstrip().split(' ')
        if line[0] not in ['track', 'fixedStep']:
            C_values.append(float(line[0]))
            Total_C+=float(line[0])
    print('Total number of ends: ' + str(Total_C))
    wigin.close()
    return C_values, Total_C


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

def write_wig(input_array, fileout_path, flag, strain_id):
    fileout=open(fileout_path, 'w+')
    fileout.write('track type=wiggle_0 name="'+str(flag)+'" autoScale=off viewLimits=0.0:25.0'+'\n'+'fixedStep chrom='+str(strain_id)+' start=1 step=1\n')   
    for i in range(len(input_array)):
        fileout.write(str(input_array[i])+'\n')
    fileout.close()    
    return

#######
#Wrapp functions for normalization, smoothing, filtering and plotting together.
#######

def NSF_wrapper(wig_dir, window, diap, del_path, strain_id, path_out):
    #Read deletions info.
    deletions=deletions_info(del_path)
    #Read data.
    wig_data=read_wig_files(wig_dir)
    #Identify minimal coverage.
    total_coverage_ar=[]
    for sample in wig_data:
        total_coverage_ar.append(sample[1])
    Min_total_coverage=min(total_coverage_ar)
    print('Min_total_coverage: ' + str(Min_total_coverage))    
    #Normalize samples coverage (add pseudocount +1).
    print("Samples are normalizing...")
    wig_data_norm=[]
    i=0
    for sample in wig_data:
        print("Now sample " + str(i+1) + " out of " + str(len(wig_data)) + " is normalizing...")
        sample_norm=[1.0 * (x + 1) * Min_total_coverage/sample[1] for x in sample[0]]
        wig_data_norm.append(sample_norm)
        i+=1
        
    #Smooth normalized coverage.
    print("Samples are smoothing...")
    wig_data_norm_sm=[]
    i=0
    for sample_norm in wig_data_norm:
        print("Now sample " + str(i+1) + " out of " + str(len(wig_data)) + " is smoothing...")
        sample_norm_sm=Smoothing(sample_norm, deletions, window)
        wig_data_norm_sm.append(sample_norm_sm)
        i+=1
    
    #Data Fourier-filtration.
    print("Samples are filtering...")
    ff=Fourier()
    X=np.arange(0,len(wig_data_norm_sm[0]))
    
    retained_d={}
    for hnum in range(diap):
        retained_d[int(hnum)]=None 
    retained=list(retained_d)
    
    i=0
    for sample_norm_sm in wig_data_norm_sm:
        print("Now sample " + str(i+1) + " out of " + str(len(wig_data)) + " is filtering...")
        ft=ff.fourier(np.array(sample_norm_sm))
        ft_reduced=ff.reduce(ft, retained)
        print(ft_reduced)
        sample_norm_filtered=ff.recover(X, ft_reduced)
        print(X[:10])
        print(sample_norm_sm[:10])
        print(sample_norm_filtered[:10])
        #Plot the data.
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print("Sample is plotting...")
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        plt.suptitle("E. coli synchronization data filtration 2000bp sm 20h filtering", fontsize=20)
        plot1=plt.subplot() 
        plot1.plot(X, sample_norm_sm, '--', label='Original', color='#550000', linewidth=1)
        plot1.plot(X, sample_norm_filtered, '--', label='Filtered 20 fh', color='#D17E7E', linewidth=1)  
        plot1.set_xlabel('Genome position, nt', size=17)
        plot1.set_ylabel('Coverage depth', size=17)
        plot1.legend(loc='upper left')
        plt.show()
        plt.savefig(path_out + "Figures\+IP+Cfx_R2_" + str(i) + "_1000bp_sm_20fh_filt.png", dpi=300, figsize=(16, 8))
        plt.close()
        #Write the WIG files.
        write_wig(sample_norm_sm, path_out + 'WIG_files\\' + str(i) + "_+IP+Cfx_R2_normalized_smoothed_1000bp.wig", i, strain_id)
        write_wig(sample_norm_filtered, path_out + 'WIG_files\\' + str(i) + "_+IP+Cfx_R2_normalized_smoothed_1000bp_filtered_20fh.wig", i, strain_id)
        i+=1
        
        
    
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
    return

#NSF_wrapper(WIG_input, 1000, 20, Deletions, "NC_007779.1_w3110_Mu", Output_directory)

#######
#Wrapp functions for normalization of +IP+Cfx samples by corresponding -IP+Cfx and plotting.
#######

def NP_wrapper(wig_dir_control, wig_dir_exp, del_path, strain_id, path_out):
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
    #Data plotting.
    X=np.arange(0,len(wig_data_exp_div[0]))
    for i in range(len(wig_data_exp_div)):
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print("Sample is plotting...")
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        plt.suptitle("E. coli synchronization data divided +IP+Cfx/-IP+Cfx", fontsize=20)
        plot1=plt.subplot() 
        plot1.plot(X, wig_data_exp_div[i], '-', label='+IP+Cfx_NS/-IP+Cfx_NSF', color='#550000', linewidth=0.2)
        plot1.plot(X, np.array(wig_data_controls[i][0])/550, '-', label='-IP+Cfx NSF', color='#D17E7E', linewidth=3)
        plot1.set_xlabel('Genome position, nt', size=17)
        plot1.set_ylabel('Coverage depth normalized', size=17)
        plot1.legend(loc='upper left')
        plt.show()
        plt.savefig(path_out  + str(i) + "_+IP+Cfx_R2_NS_div_-IP+Cfx_R2_NSF_-IP+Cfx_R2_NSF.png", dpi=300, figsize=(16, 8))
        plt.close()
    return

#NP_wrapper(WIG_control_NSF_input, WIG_exp_NS_input, Deletions, "NC_007779.1_w3110_Mu", Output_directory_div)

#######
#Wrapp functions for sequential division of time-points t(i+1)/t(i).
#######

def Seq_div_wrapper(wig_dir_control, wig_dir_exp, diap, del_path, strain_id, path_out):
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
    #Sequental division t(i+1)/t(i)
    wig_data_seq_div=[]
    #for i in range(len(wig_data_exp_div)-1): #For real sequential division t(i+1)/t(i).
    for i in range(len(wig_data_exp_div)):  #For division on some fixed time-point (e.g. 10 min).
        seq_div=[]
        for j in range(len(wig_data_exp_div[i])):
            #j_seq_div=wig_data_exp_div[i+1][j]/wig_data_exp_div[i][j] #For real sequential division t(i+1)/t(i).
            j_seq_div=wig_data_exp_div[i][j]/wig_data_exp_div[4][j] #For division on some fixed time-point (e.g. 15 min).
            seq_div.append(j_seq_div)
        wig_data_seq_div.append(seq_div)
        
    #Data Fourier-filtration.
    print("Samples are filtering...")
    ff=Fourier()
    X=np.arange(0,len(wig_data_seq_div[0]))
    
    retained_d={}
    for hnum in range(diap):
        retained_d[int(hnum)]=None 
    retained=list(retained_d)
    
    #Data plotting.
    for i in range(len(wig_data_seq_div)):
        #Fourier transformations.
        print("Now sample " + str(i+1) + " out of " + str(len(wig_data_seq_div)) + " is filtering...")
        ft=ff.fourier(np.array(wig_data_seq_div[i]))
        ft_reduced=ff.reduce(ft, retained)
        print(ft_reduced)
        wig_data_seq_div_filtered=ff.recover(X, ft_reduced)
        print(X[:10])
        print(wig_data_seq_div[i][:10])
        print(wig_data_seq_div_filtered[:10])        
        #Some hack to avoid some bug in matplotlib (OverflowError: In draw_path: Exceeded cell block limit)
        #See: https://stackoverflow.com/questions/37470734/matplotlib-giving-error-overflowerror-in-draw-path-exceeded-cell-block-limit
        print("Sample is plotting...")
        mpl.rcParams['agg.path.chunksize']=10000
        plt.figure(figsize=(16, 8), dpi=100)
        plt.suptitle("E. coli synchronization data sequential division i/15 min", fontsize=20)
        plot1=plt.subplot() 
        plot1.plot(X, wig_data_seq_div[i], '-', label='(+IP+Cfx_NS/-IP+Cfx_NSF)i+1/(+IP+Cfx_NS/-IP+Cfx_NSF)15 min Original', color='#550000', linewidth=0.2)
        plot1.plot(X, wig_data_seq_div_filtered, '--', label='(+IP+Cfx_NS/-IP+Cfx_NSF)i+1/(+IP+Cfx_NS/-IP+Cfx_NSF)15 min Filtered', color='#D17E7E', linewidth=3)
        plot1.set_xlabel('Genome position, nt', size=17)
        plot1.set_ylabel('Coverage depth normalized', size=17)
        plot1.legend(loc='upper left')
        plt.show()
        plt.savefig(path_out  + "Figures\Div_15_min\\" + str(i) + "_div_" +str('15 min')+ "_plus_20fh_R2.png", dpi=300, figsize=(16, 8))
        plt.close()
        #Write the WIG files.
        write_wig(wig_data_seq_div[i], path_out + "WIG_files\Div_15_min\\" + str(i) + "_div_" + str('15 min') + "_R2.wig", str(i)+"/"+str('15 min'), strain_id)
        write_wig(wig_data_seq_div_filtered, path_out + "WIG_files\Div_15_min\\" + str(i) + "_div_" + str('15 min') + "_20fh_R2.wig", str(i)+"/"+str('15 min'), strain_id)        
    return

Seq_div_wrapper(WIG_control_NSF_input, WIG_exp_NS_input, 20, Deletions, "NC_007779.1_w3110_Mu", Output_directory_div)


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
#Wrapp functions for replisome and gyrase enrichment positions identification.
#######

def Repl_Gyr_ident(wig_dir_control_NFS, wig_dir_exp_seq_div, diap, del_path, strain_id, path_out):
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

#Repl_Gyr_ident(WIG_control_NSF_input, WIG_exp_NSSDF_input, 20, Deletions, "NC_007779.1_w3110_Mu", Output_directory_div)