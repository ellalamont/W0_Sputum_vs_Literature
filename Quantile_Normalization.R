# Quantile Normalization on my W0 TPM values
# 4/16/25

# Following Coppola et al. 2021
# https://www.frontiersin.org/journals/immunology/articles/10.3389/fimmu.2021.763364/full#h13

source("Import_data.R")

# my_tpm_W0Sputum
# Won't work if starting with averages! Needs multiple columns

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

test <- as.data.frame(quantile_normalisation(my_tpm_W0Sputum))
class(test)

#########################################################
############# AVERAGE THE NEW RANK SCORES ###############

# Get the average tpm for my 3 week 0s
RANK_Average_tpm_W0Sputum <- test %>% 
  mutate(RANK_Average = rowMeans(across(where(is.numeric)))) %>%
  select(RANK_Average)



#########################################################
################### MAKE IT 0 TO 100 ####################
# Coppola2021 def did this but didn't seem to say how.......

# ChatGPT made this
rank_cols_to_percentile <- function(df) {
  apply(df, 2, function(col) {
    ranks <- rank(col, ties.method = "min")
    percentiles <- (ranks - 1) / (length(ranks) - 1) * 100
    return(percentiles)
  }) |> as.data.frame()
}

test_2 <- rank_cols_to_percentile(test)

# Get the average for my 3 week 0s
RANK_Average_tpm_W0Sputum <- test_2 %>% 
  mutate(RANK_Average = rowMeans(across(where(is.numeric)))) %>%
  select(RANK_Average)





# TESTING BELOW.....
#########################################################
###### QUANTILE NORMALISATION OUTSIDE OF FUNCTION #######

# Give all the genes a rank
df_rank <- apply(my_tpm_W0Sputum, 2, rank, ties.method="min")

# Sort all the genes based on original data (ascending is default)
df_sorted <- data.frame(apply(my_tpm_W0Sputum, 2, sort))

# Calculate the average TPM for each row
df_mean <- apply(df_sorted, 1, mean)

# New Function: takes the original ranks and replaces them with the mean values we just calculated
index_to_mean <- function(my_index, my_mean) {return (my_mean [my_index])}

# Replaces the values in each column using their original rank and substitute in the mean value for that rank.
# So if a gene was ranked 10th in a sample, it now gets the average of all the 10th-ranked values across samples.
df_final <- apply(df_rank, 2, index_to_mean, my_mean = df_mean)

rownames(df_final) <- rownames(df)
