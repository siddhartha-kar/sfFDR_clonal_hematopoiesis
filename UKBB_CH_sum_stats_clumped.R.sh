### Set working directory and load required packages:
setwd("/Users/cs21140/Downloads/ch_cross_trait")
library(vroom)
library(tidyverse)

### Load data:
overall <- vroom("overall_GCST90165261_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
overall_inc <- vroom("overall_GCST90165267_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
dnmt3a <- vroom("dnmt3a_GCST90165271_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
tet2 <- vroom("tet2_GCST90165281_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")
asxl1 <- vroom("asxl1_GCST90165259_buildGRCh37_selectedcols_0.5pcMAF_rsid.tsv.gz")

### Format data for clumping:
overall <- overall %>%
  rename(SNPID = SNP,
         SNP = rsid,
         P = pval_overall)
vroom_write(overall, file = "overall_for_clumping.tsv.gz")

overall_inc <- overall_inc %>%
  rename(SNPID = SNP,
         SNP = rsid,
         P = pval_overall)
vroom_write(overall_inc, file = "overall_inc_for_clumping.tsv.gz")

dnmt3a <- dnmt3a %>%
  rename(SNPID = SNP,
         SNP = rsid,
         P = pval_dnmt3a)
vroom_write(dnmt3a, file = "dnmt3a_for_clumping.tsv.gz")

tet2 <- tet2 %>%
  rename(SNPID = SNP,
         SNP = rsid,
         P = pval_tet2)
vroom_write(tet2, file = "tet2_for_clumping.tsv.gz")

asxl1 <- asxl1 %>%
  rename(SNPID = SNP,
         SNP = rsid,
         P = pval_asxl1)
vroom_write(asxl1, file = "asxl1_for_clumping.tsv.gz")

### Go to terminal and clump:
cd Downloads/ch_cross_trait
chmod +x plink2
./plink2 \
--bfile g1000_eur \
--clump-p1 1 \
--clump-r2 0.02 \
--clump-kb 1000 \
--clump overall_for_clumping.tsv.gz  \
--clump-snp-field SNP \
--clump-field P \
--out overall_clumped
awk 'NR!=1{print $3}' overall_clumped.clumps > overall_clumped.valid.snp 

./plink2 \
--bfile g1000_eur \
--clump-p1 1 \
--clump-r2 0.02 \
--clump-kb 1000 \
--clump overall_inc_for_clumping.tsv.gz  \
--clump-snp-field SNP \
--clump-field P \
--out overall_inc_clumped
awk 'NR!=1{print $3}' overall_inc_clumped.clumps > overall_inc_clumped.valid.snp 

./plink2 \
--bfile g1000_eur \
--clump-p1 1 \
--clump-r2 0.02 \
--clump-kb 1000 \
--clump dnmt3a_for_clumping.tsv.gz \
--clump-snp-field SNP \
--clump-field P \
--out dnmt3a_clumped
awk 'NR!=1{print $3}' dnmt3a_clumped.clumps > dnmt3a_clumped.valid.snp 

./plink2 \
--bfile g1000_eur \
--clump-p1 1 \
--clump-r2 0.02 \
--clump-kb 1000 \
--clump tet2_for_clumping.tsv.gz  \
--clump-snp-field SNP \
--clump-field P \
--out tet2_clumped
awk 'NR!=1{print $3}' tet2_clumped.clumps > tet2_clumped.valid.snp 

./plink2 \
--bfile g1000_eur \
--clump-p1 1 \
--clump-r2 0.02 \
--clump-kb 1000 \
--clump asxl1_for_clumping.tsv.gz  \
--clump-snp-field SNP \
--clump-field P \
--out asxl1_clumped
awk 'NR!=1{print $3}' asxl1_clumped.clumps > asxl1_clumped.valid.snp 