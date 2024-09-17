#!/bin/bash -l

#conda activate micov

genome_len=/projects/wol/qiyun/rs210/length.map
path_list=list_compress_ouput.txt

micov consolidate \
	--lengths "$genome_len" \
	--paths "$path_list" \
	--output consolidated.tgz

echo 'done'
