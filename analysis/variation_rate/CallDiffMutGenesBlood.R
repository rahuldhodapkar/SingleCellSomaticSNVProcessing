#!/usr/bin/env Rscript
#
# Call differentially mutated genes from high quality variant calls of coding variants
# between MS patients and healthy controls in blood.
#

library(dplyr)
library(RPostgreSQL)
library(ggplot2)

# From `by_hand_queries.sql`
TOTAL_NUM_CELLS_MS = 758;
TOTAL_NUM_CELLS_HC = 763;

TOTAL_NUM_DONORS_MS = 8;
TOTAL_NUM_DONORS_HC = 8;

# loads the PostgreSQL driver
drv <- dbDriver("PostgreSQL")
# creates a connection to the postgres database
# note that "con" will be used later in each connection to the database
con <- dbConnect(drv, dbname = "ms",
                 host = "localhost", port = 5432)

vdf_ms <- dbGetQuery(con, '
  SELECT * from variants_by_gene
  WHERE 
      "Group" = \'MS\'
      AND "TissueType"=\'B\';
');

vdf_hc <- dbGetQuery(con, '
  SELECT * from variants_by_gene
  WHERE 
      "Group" = \'HC\'
      AND "TissueType"=\'B\';
');

# Full Join dataframes
vdf_merge <- vdf_ms %>%
  full_join(vdf_hc, by = 'Gene', suffix=c('.ms', '.hc'))

################################################################
## Differences By Number of Cells
################################################################

df_cell <- data.frame(
  Gene=         vdf_merge$Gene,
  NumCells.ms=  vdf_merge$NumCells.ms,
  NumCells.hc=  vdf_merge$NumCells.hc,
  NumDonors.ms= vdf_merge$NumDonors.ms,
  NumDonors.hc= vdf_merge$NumDonors.hc
);

df_cell[is.na(df_cell)] <- 0

# Generate p-value for each gene
p_vals <- c()
estimates <- c()

for (i in 1:nrow(df_cell)) {
  x <- fisher.test(matrix(c(
    df_cell$NumCells.ms[i], TOTAL_NUM_CELLS_MS - df_cell$NumCells.ms[i],
    df_cell$NumCells.hc[i], TOTAL_NUM_CELLS_HC - df_cell$NumCells.hc[i]
  ), nrow = 2, ncol = 2))
  
  p_vals <- c(p_vals, x$p.value)
  estimates <- c(estimates, x$estimate)
}

df_cell$p <- p_vals;
df_cell$or_est <- estimates;

df_cell$log2or <- log2(df_cell$or_est);
df_cell$negLog10p <- -log10(df_cell$p)

write.csv(df_cell, './ad_hoc/diff_cell_blood.csv', row.names = F);

################################################################
## Differences By Number of Donors
################################################################

df_donor <- data.frame(
  Gene=         vdf_merge$Gene,
  NumCells.ms=  vdf_merge$NumCells.ms,
  NumCells.hc=  vdf_merge$NumCells.hc,
  NumDonors.ms= vdf_merge$NumDonors.ms,
  NumDonors.hc= vdf_merge$NumDonors.hc
);

df_donor[is.na(df_donor)] <- 0

# Generate p-value for each gene
p_vals <- c()
estimates <- c()

for (i in 1:nrow(df_donor)) {
  x <- fisher.test(matrix(c(
    df_donor$NumDonors.ms[i], TOTAL_NUM_DONORS_MS - df_donor$NumDonors.ms[i],
    df_donor$NumDonors.hc[i], TOTAL_NUM_DONORS_HC - df_donor$NumDonors.hc[i]
  ), nrow = 2, ncol = 2))
  
  p_vals <- c(p_vals, x$p.value)
  estimates <- c(estimates, x$estimate)
}

df_donor$p <- p_vals;
df_donor$or_est <- estimates;

df_donor$log2or <- log2(df_donor$or_est);
df_donor$negLog10p <- -log10(df_donor$p)

write.csv(df_donor, './ad_hoc/diff_donor_blood.csv', row.names = F);

