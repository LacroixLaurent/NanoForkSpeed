# NanoForkSpeed
### Laurent Lacroix (laurent.lacroix@inserm.fr)
### BrdU Basecalling

*basecalling_sample.sh* contains a example of the base calling procedure going from RawData to the megalodon mod_mappings.bam file.  
*basecalling_sample.r* proceeds in R from the mod_mappings.bam file to smoothed signal in a .rds format (to open with R) used for forks detection. In this script the bam file is first split into 50 000 reads bam files to allow better memory handling. This process requires samtools (we used v1.10).  
*Parsing_function4megalodon.r* contains the function used in the basecalling process going from the modified bam file created by megalodon to the parsed data used to detect th forks. The parsed data contain both the raw data and the smoothed data (see publication). This script also annotated gaps that can be introduced during the alignment step by megalodon. These gaps will be used during the forks detection process to exclude wrongly detected forks due to a signal discontinuity introduced by those gaps.      

*session_info()* : see *base_calling_session_info.txt*
