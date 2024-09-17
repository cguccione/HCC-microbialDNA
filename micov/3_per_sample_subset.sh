#!/bin/bash

#conda activate micov

micov per-sample-group \
	--qiita-coverages /panfs/cguccion/22_06_22_HCC_CRC_Amir/micov/rs210_clean/consolidated.tgz \
	--lengths /projects/wol/qiyun/rs210/length.map \
	--sample-metadata micov_meta_pangenome_rs210clean.tsv \
	--sample-metadata-column sample_type \
	--output /panfs/cguccion/22_06_22_HCC_CRC_Amir/micov/rs210_clean/all_microbes_/all_microbes

#	--features-to-keep of_interest_birdman.tsv 

echo 'done'

