#!bin/bash
echo '
###########
#Sutormin Dmitry, 2019
#Equalization of coverage depth for Time-course gyrase Topo-Seq data.
###########
'

head -28710948 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_9_S89_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_9_S89_eq_R1_001.fastq
head -28710948 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_9_S89_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_9_S89_eq_R2_001.fastq
head -28452108 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_10_S90_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_10_S90_eq_R1_001.fastq
head -28452108 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_10_S90_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_10_S90_eq_R2_001.fastq
head -30147488 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_11_S91_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_11_S91_eq_R1_001.fastq
head -30147488 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_11_S91_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_11_S91_eq_R2_001.fastq
head -28527256 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_12_S92_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_12_S92_eq_R1_001.fastq
head -28527256 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_12_S92_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_12_S92_eq_R2_001.fastq
head -31287640 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_13_S93_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_13_S93_eq_R1_001.fastq
head -31287640 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_13_S93_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_13_S93_eq_R2_001.fastq
head -29973224 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_14_S94_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_14_S94_eq_R1_001.fastq
head -29973224 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_14_S94_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_14_S94_eq_R2_001.fastq
head -29841900 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_15_S95_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_15_S95_eq_R1_001.fastq
head -29841900 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_15_S95_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_15_S95_eq_R2_001.fastq
head -32661124 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_16_S96_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_16_S96_eq_R1_001.fastq
head -32661124 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_16_S96_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_16_S96_eq_R2_001.fastq
head -28347504 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_17_S97_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_17_S97_eq_R1_001.fastq
head -28347504 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_17_S97_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_17_S97_eq_R2_001.fastq
head -28326876 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_18_S98_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_18_S98_eq_R1_001.fastq
head -28326876 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_18_S98_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_18_S98_eq_R2_001.fastq
echo '10 samples passed!'
head -28777500 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_19_S99_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_19_S99_eq_R1_001.fastq
head -28777500 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_19_S99_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_19_S99_eq_R2_001.fastq
head -28267752 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_20_S100_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_20_S100_eq_R1_001.fastq
head -28267752 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_20_S100_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_20_S100_eq_R2_001.fastq
head -28651252 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_21_S101_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_21_S101_eq_R1_001.fastq
head -28651252 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_21_S101_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_21_S101_eq_R2_001.fastq
head -29780980 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_22_S102_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_22_S102_eq_R1_001.fastq
head -29780980 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_22_S102_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_22_S102_eq_R2_001.fastq
head -31533248 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_23_S103_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_23_S103_eq_R1_001.fastq
head -31533248 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_23_S103_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_23_S103_eq_R2_001.fastq
head -32452096 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_24_S104_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_24_S104_eq_R1_001.fastq
head -32452096 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_24_S104_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_24_S104_eq_R2_001.fastq
head -28738032 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_33_S105_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_33_S105_eq_R1_001.fastq
head -28738032 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_33_S105_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_33_S105_eq_R2_001.fastq
head -28725728 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_34_S106_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_34_S106_eq_R1_001.fastq
head -28725728 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_34_S106_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_34_S106_eq_R2_001.fastq
head -29572188 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_35_S107_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_35_S107_eq_R1_001.fastq
head -29572188 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_35_S107_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_35_S107_eq_R2_001.fastq
head -28505820 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_36_S108_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_36_S108_eq_R1_001.fastq
head -28505820 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_36_S108_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_36_S108_eq_R2_001.fastq
echo '20 samples passed!'
head -28890516 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_37_S109_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_37_S109_eq_R1_001.fastq
head -28890516 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_37_S109_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_37_S109_eq_R2_001.fastq
head -29956780 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_38_S110_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_38_S110_eq_R1_001.fastq
head -29956780 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_38_S110_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_38_S110_eq_R2_001.fastq
head -31995044 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_39_S111_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_39_S111_eq_R1_001.fastq
head -31995044 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_39_S111_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_39_S111_eq_R2_001.fastq
head -32738160 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_40_S112_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_40_S112_eq_R1_001.fastq
head -32738160 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_40_S112_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_40_S112_eq_R2_001.fastq
head -28367156 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_41_S113_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_41_S113_eq_R1_001.fastq
head -28367156 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_41_S113_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_41_S113_eq_R2_001.fastq
head -28299856 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_42_S114_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_42_S114_eq_R1_001.fastq
head -28299856 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_42_S114_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_42_S114_eq_R2_001.fastq
head -28790024 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_43_S115_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_43_S115_eq_R1_001.fastq
head -28790024 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_43_S115_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_43_S115_eq_R2_001.fastq
head -28291336 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_44_S116_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_44_S116_eq_R1_001.fastq
head -28291336 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_44_S116_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_44_S116_eq_R2_001.fastq
head -28710692 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_45_S117_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_45_S117_eq_R1_001.fastq
head -28710692 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_45_S117_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_45_S117_eq_R2_001.fastq
head -29828612 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_46_S118_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_46_S118_eq_R1_001.fastq
head -29828612 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_46_S118_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_46_S118_eq_R2_001.fastq
echo '30 samples passed!'
head -31394476 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_47_S119_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_47_S119_eq_R1_001.fastq
head -31394476 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_47_S119_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_47_S119_eq_R2_001.fastq
head -32486204 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_48_S120_R1_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_48_S120_eq_R1_001.fastq
head -32486204 /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_unzipped/DSu_48_S120_R2_001.fastq > /home/cls01/Data_Hi-C/Topo-Seq_data/Time-course_gyrase/Raw_data_eq/DSu_48_S120_eq_R2_001.fastq
echo 'Done!'
