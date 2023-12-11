qadabra add-dataset \
    --workflow-dest /home/cguccion/packages/qadabra/my_qadabra \
    --table adj_HCC_v_CRC_filtered.biom --metadata /panfs/cguccion/22_06_22_HCC_CRC_Amir/HCC-microbialDNA/processed_data/metadata/metadata_adj_HCC_v_CRC.tsv \
    --name adj_HCC_v_CRC_qiita15336_prep16181_pangenome_wol2_scrubbed_zebraFilter0.1 \
    --factor-name tumor_type \
    --target-level HCC \
    --reference-level CRC \
    --verbose
