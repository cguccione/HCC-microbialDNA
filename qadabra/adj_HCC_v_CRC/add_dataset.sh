qadabra add-dataset \
    --workflow-dest ~/Dropbox/current_Fall23/HCC_Amir/HCC-microbialDNA/qadabra/git_install/qadabra/qadabra \
    --table adj_HCC_v_CRC_filtered.biom \
    --metadata ~/Dropbox/current_Fall23/HCC_Amir/HCC-microbialDNA/processed_data/metadata/metadata_adj_HCC_v_CRC.tsv \
    --name adj_HCC_v_CRC_qiita15336_prep16181_pangenome_wol2_scrubbed_zebraFilter0.1 \
    --factor-name tumor_type \
    --target-level HCC \
    --reference-level CRC \
    --verbose
