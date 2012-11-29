parallel -u <<EOF
./computeKeysPdfSyst -m 750 -i ../Extractor2Dataset/MTT_Data_merged_2012.root --signal ../Extractor2Dataset/MTT_Zprime_750_Narrow_2012_08Nov_merged.root
./computeKeysPdfSyst -m 1000 -i ../Extractor2Dataset/MTT_Data_merged_2012.root --signal ../Extractor2Dataset/MTT_Zprime_1000_Narrow_2012_08Nov_merged.root
./computeKeysPdfSyst -m 1250 -i ../Extractor2Dataset/MTT_Data_merged_2012.root --signal ../Extractor2Dataset/MTT_Zprime_1250_Narrow_2012_08Nov_merged.root
./computeKeysPdfSyst -m 1500 -i ../Extractor2Dataset/MTT_Data_merged_2012.root --signal ../Extractor2Dataset/MTT_Zprime_1500_Narrow_2012_08Nov_merged.root
EOF
