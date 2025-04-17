# Quantile Normalization on my W0 TPM values
# 4/16/25

## **** NOT WORKING ******

# Following Coppola et al. 2021
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.763364/full#h13

source("Import_data.R")

# Average_tpm_W0Sputum

################################################
###### FUNCTION FROM COPPOLA ET AL. 2021 #######

quantile_normalisation <- function(df)
  
{df_rank <- apply (df,2, rank, ties.method="min")

df_sorted <- data.frame (apply (df, 2, sort))

df_mean <- apply (df_sorted, 1, mean)

index_to_mean <- function (my_index, my_mean) {return (my_mean [my_index])}

df_final <- apply (df_rank, 2, index_to_mean, my_mean = df_mean)

rownames (df_final) <- rownames(df)

return(df_final)}

test <- quantile_normalisation(Average_tpm_W0Sputum)
# This looks the same.....

