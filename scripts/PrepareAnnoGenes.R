#!/usr/bin/env Rscript
#
# Build pre-joined dataframes for easy analytics access. Produces the "merged_vcf_data.csv" file.
#
# Usage:
#
#   ./PrepareAnnoGenes.R <vcf_dir> <merged_df_filename.csv>
#

args = commandArgs(trailingOnly=TRUE)

vcf_dir <- args[1];
merged_df_filename <- args[2];

################################################################
## EXTRACT ALL VCF DATA
################################################################

file_names <- list.files(vcf_dir);

vcf_df <- NULL;
for (file_name in file_names) {
  vcf <- read.vcfR(paste(vcf_dir, file_name, sep='/'), verbose = FALSE);
  
  tidyVcf <- vcfR2tidy(vcf);
  
  d1 <- data.frame("CellName" = rep(file_name, nrow(vcf)));
  d2 <- tidyVcf$fix;
  d3 <- tidyVcf$gt;
  df <- cbind(d1, d2, d3);
  rm(d1, d2, d3)
  
  if( is.null(vcf_df)) {
    vcf_df <- df;
  } else {
    vcf_df <- rbind(vcf_df, df);
  }
}

filename_pattern <- "^(\\d+)([BT])Treg.*";
proto <- data.frame(DonorNum=integer(), DonorType=character());
captured_df <- strcapture(filename_pattern, vcf_df$CellName, proto);

vcf_df$DonorNum <- captured_df$DonorNum;
vcf_df$DonorType <- captured_df$DonorType;

vcf_df_dedup <- subset(vcf_df, select=which(!duplicated(names(vcf_df)))) 

write.csv(vcf_df_dedup, file = merged_df_filename, row.names = F);

