#!/bin/bash -l

#SBATCH -J compress
#SBATCH --mail-type=ALL
#SBATCH --mail-user=cguccion@ucsd.edu
#SBATCH -t 1:00:00

conda activate micov

in_folder=/panfs/cguccion/22_06_22_HCC_CRC_Amir/woltka_rs210_clean/pese_pangenome_align-RS210_masked/sam_only
out_folder=/panfs/cguccion/22_06_22_HCC_CRC_Amir/micov/rs210_clean/compress_output

for file in "$in_folder"/*.sam; do
	if [ -e "$file" ]; then
		temp_fn=$(basename "$file" .sam)
		cat "$file" | micov compress > "$out_folder/${temp_fn}.cov"
	fi
done

echo 'done'

