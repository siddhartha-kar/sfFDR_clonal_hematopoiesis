#CREDITS: https://github.com/Y-Wentao/mdS2G?tab=readme-ov-file

library("mds2g")

cs2g_link(
    pheno_name = c("CH"),
    SNP_data = "./sffdr_ch_gwas_causal.csv",                
    s2g_folder = "./S2G_bed/",    
    save_folder = "./output/",                      
    keep_intermediate_file = TRUE,
    intermediate_file_folder = "./output/"
)