#!/bin/bash -l
#SBATCH --partition=short
#SBATCH --mail-user="cguccion@ucsd.edu"
#SBATCH --mail-type=ALL
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --partition=short
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --array=1-20

#Change variables per run
group_fn=blood_HCC_v_CRC
#Everything else is set on this

main_path=/panfs/cguccion/22_06_22_HCC_CRC_Amir/HCC-microbialDNA/birdman/${group_fn}
cd ${main_path}

pwd; hostname; date

set -e

conda activate birdman

echo Chunk $SLURM_ARRAY_TASK_ID / $SLURM_ARRAY_TASK_MAX

TABLEID=${group_fn}_qiita15336_prep16181_pangenome_wol2_scrubbed_zebraFilter0.1
TABLE=/panfs/cguccion/22_06_22_HCC_CRC_Amir/HCC-microbialDNA/processed_data/biom/${TABLEID}/feature-table.biom
OUTDIR=${main_path}/inferences/${TABLEID}
LOGDIR=${main_path}/logs/${TABLEID}

mkdir -p $OUTDIR
mkdir -p $LOGDIR
mkdir -p ${main_path}/inferences-results

echo Starting Python script...
time python /panfs/cguccion/22_06_22_HCC_CRC_Amir/HCC-microbialDNA/birdman/src/${group_fn}_birdman_chunked.py \
    --table-path $TABLE \
    --inference-dir $OUTDIR \
    --num-chunks $SLURM_ARRAY_TASK_MAX \
    --chunk-num $SLURM_ARRAY_TASK_ID \
    --logfile "${LOGDIR}/chunk_${SLURM_ARRAY_TASK_ID}.log" && echo Finished Python script!

