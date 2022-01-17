# Sample of the basecalling procedure used with our BrdU trained model
# This requires the ONT megalodon and guppy software
# We used guppy v4.4.1 GPU and megalogon v2.2.9
# either the ont_guppy folder is present in the working directory
# or a symbolic link for guppy should be present in the working directory
# ont-guppy -> ~/src/ont-guppy_4.4.1
ln -s ~/src/ont-guppy_4.4.1 ./ont-guppy
# the configuration file should point to the modified base model
# see BrdU_Configutation4megalodon.cfg
# (the model file is provided  as a compressed archive (model_adversial_bug.json.zip))
# the cfg file should be placed in the src/ont-guppy_4.4.1/data folder (not sure it is mandatory)
# a conda env was create using the yml file provided within miniconda3
conda env create -f megalodon.yml
# then within this env, we used megalodon with the following parameters
conda activate mega
REF=./data/S288CwExtrarDNA_ROMAN.fa
# we used a modified version of the S288C reference genome (R64.2.1)
# to which we add an extra chromosome (rDNA-10R) containing 10 rDNA tandem repats
CFG=BrdU_Configuration4megalodon.cfg
INPUT=./data/EXP_fast5
OUTPUT=./data/EXP_mega

nohup megalodon $INPUT --outputs mod_mappings --reference $REF  --output-dir $OUTPUT  --guppy-config $CFG --disable-mod-calibration --overwrite --device cuda:0 --processes 20 --guppy-timeout 90 --mod-min-prob 0  &
# please notice than to experiment with many long reads, we had to increase the guppy-timeout to 600
# also notice that this computer was equiped with a NVidia RTX2080Ti GPU

# this generate a megalon's style mod-mappings.bam file in the EXP_mega folder
# this bam file can be big and we had to split it to be able to process it with our parsing script
# the following step are done within R (v4.0.5) on the same computer used for the basecalling
# with the script basecalling_sample.r

