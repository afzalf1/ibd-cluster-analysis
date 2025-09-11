#install.packages("readr")
library(readr)
library(MASS)


all_data <- read_tsv("Insert_file_name_here")

#Subset of data that has no unknown complications/phenotypes (loses 10%)
data <- subset(all_data, `Crohn s disease phenotype` != "Unknown")

#Renaming columns with spaces
names(data)[names(data) == "Crohn s disease phenotype"] <- "Phenotype"
names(data)[names(data) == "IBD surgery final"] <- "Surgery"

#Collapse clusters with less than 5 patients into a Miscellaneous cluster
cluster_counts <- table(data$Cluster)

data$Cluster <- as.character(data$Cluster)
data$Cluster[data$Cluster %in% names(cluster_counts[cluster_counts < 5])] <- "Miscellaneous"
data$Cluster <- factor(data$Cluster)

#Creating the 4-way contingency table
contingency_table <- table(
  Cluster = data$Cluster,
  Diagnosis = data$Diagnosis,
  Complication = data$Phenotype,
  Surgery = data$Surgery
)

#Find whether there is any dependence/independence in the whole table

#Finds the expected values
independence_model <- loglin(contingency_table, list(1, 2, 3, 4), fit = TRUE)
expected_values <- independence_model$fit

#Gets the chi-square statistic, the degrees of freedom and the p-value between observed and expected values
chi_square <- sum((contingency_table - expected_values)^2 / expected_values)
print(paste("Pearson Chi-square statistic:", chi_square))

df = independence_model$df
print(paste("Degrees of freedom:", df))

p_value <- 1 - pchisq(chi_square, df)
print(paste("P-value:", p_value))

#Computes the ratios observed/expected and ln(ratios)
ratios <- contingency_table / expected_values
log_ratios <- log(ratios)

#Find deviations
miss_cluster_diagnosis <- loglin(contingency_table, list(c(1,3), c(1,4), c(2,3), c(2,4), c(3, 4)), fit=TRUE)
miss_cluster_complication <- loglin(contingency_table, list(c(1,2), c(1,4), c(2,3), c(2,4), c(3,4)), fit=TRUE)
miss_cluster_surgery <- loglin(contingency_table, list(c(1,2), c(1,3), c(2,3), c(2,4), c(3,4)), fit=TRUE)



#FOR KD PAIRS ONLY



#Collapse over complication and surgery to get just observed K x D table
obs_KD_table <- table(data$Cluster, data$Diagnosis)

#Get expected KD_table
chi_result_KD <- chisq.test(obs_KD_table)
exp_KD_table <- chi_result_KD$expected

#Gets the chi-square statistic, the degrees of freedom and the p-value between observed and expected values
chi_square_KD <- sum((obs_KD_table - exp_KD_table)^2 / exp_KD_table)
print(paste("Pearson Chi-square statistic for cluster and diagnosis:", chi_square_KD))

df_KD = prod(dim(obs_KD_table) - 1)
print(paste("Degrees of freedom for cluster and diagnosis:", df_KD))

p_value_KD <- 1 - pchisq(chi_square_KD, df_KD)
print(paste("P-value for whether cluster and diagnosis are independent overall:", p_value_KD))

#Computes the ratios observed/expected and ln(ratios) to find how much more the cluster-diagnosis pair occurs
ratios_KD <- obs_KD_table / exp_KD_table
log_ratios_KD <- log(ratios_KD)

#Now get the table of the standard deviations of each cell
#Finds the difference in observed and expected values
difference_KD <- obs_KD_table - exp_KD_table

#Computes row and column proportions of the table
row_proportions_KD <- rowSums(obs_KD_table) / sum(obs_KD_table)
column_proportions_KD <- colSums(obs_KD_table) / sum(obs_KD_table)

#Gets variance of each cell
var_table_KD <- exp_KD_table * (1 - row_proportions_KD) %o% (1 - column_proportions_KD)

#Gets standard deviation of each cell
sd_table_KD <- sqrt(var_table_KD)

#Gets adjusted residuals of each cell
adj_resid_KD <- difference_KD / sd_table_KD

#Calculate p-value of each cell (before multiple test correction)
#First calculate the probability of each cell being any value up to the positive of the actual value in the cell
probability_up_to_value_KD <- pnorm(abs(adj_resid_KD))

#Calculate the probability of the cell being at least that extreme on BOTH sides (un adjusted p-values)
unadj_pvals_KD <- 2 * (1 - probability_up_to_value_KD)

#Now perform multiple test correction (FDR)
#FDR ranks the p-values from smallest to largest and changes the threshold for significance based upon how large
#the p-value is compared to the others (larger p-value means larger threshold)
#Turn unadjusted p-values into vector form
unadj_pvals_vec_KD <- as.vector(unadj_pvals_KD)

#Sort from smallest to largest
order_index_KD <- order(unadj_pvals_vec_KD)
sorted_unadj_p_KD <- unadj_pvals_vec_KD[order_index_KD]

#Total number of tests
total_KD <- length(unadj_pvals_vec_KD)

#Compute for each sorted p value the BH critical value
q_vals_KD <- sorted_unadj_p_KD * total_KD / seq_along(sorted_unadj_p_KD)

#Calculate cumulative min going backwards
sorted_adj_KD <- rev(cummin(rev(q_vals_KD)))

#Put adjusted values back in order
adj_pvals_KD <- numeric(total_KD)
adj_pvals_KD[order_index_KD] <- sorted_adj_KD

#Convert back into table
adj_pvals_table_KD <- matrix(adj_pvals_KD, nrow = nrow(obs_KD_table), ncol = ncol(obs_KD_table), dimnames = dimnames(obs_KD_table))
print("P-values For Cluster-Diagnosis Cells:")
print(adj_pvals_table_KD)

#Print the significant cells
sig_cells_KD <- which(adj_pvals_table_KD < 0.05, arr.ind = TRUE)
if (length(sig_cells_KD) == 0) {
  print("No significant p-values (< 0.05)")
} else {
  sig_results_KD <- data.frame(
    Cluster = rownames(adj_pvals_table_KD)[sig_cells_KD[, 1]],
    Diagnosis = colnames(adj_pvals_table_KD)[sig_cells_KD[, 2]],
    Adj_Pval = adj_pvals_table_KD[sig_cells_KD],
    Ratio_OE = ratios_KD[sig_cells_KD]
  )
  print("Significant Cluster-Diagnosis Associations (p < 0.05):")
  print(sig_results_KD)
}

#Print O/E ratios
print("Observed/Expected (O/E) Ratios:")
print(ratios_KD)



#FOR KC PAIRS ONLY



#Collapse over complication and complication/phenotype to get just observed K x C table
obs_KC_table <- table(data$Cluster, data$Phenotype)

#Get expected KC_table
chi_result_KC <- chisq.test(obs_KC_table)
exp_KC_table <- chi_result_KC$expected

#Gets the chi-square statistic, the degrees of freedom and the p-value between observed and expected values
chi_square_KC <- sum((obs_KC_table - exp_KC_table)^2 / exp_KC_table)
print(paste("Pearson Chi-square statistic for cluster and diagnosis:", chi_square_KC))

df_KC = prod(dim(obs_KC_table) - 1)
print(paste("Degrees of freedom for cluster and complication:", df_KC))

p_value_KC <- 1 - pchisq(chi_square_KC, df_KC)
print(paste("P-value for whether cluster and complication are independent overall:", p_value_KC))

#Computes the ratios observed/expected and ln(ratios) to find how much more the cluster-diagnosis pair occurs
ratios_KC <- obs_KC_table / exp_KC_table
log_ratios_KC <- log(ratios_KC)

#Now get the table of the standard deviations of each cell
#Finds the difference in observed and expected values
difference_KC <- obs_KC_table - exp_KC_table

#Computes row and column proportions of the table
row_proportions_KC <- rowSums(obs_KC_table) / sum(obs_KC_table)
column_proportions_KC <- colSums(obs_KC_table) / sum(obs_KC_table)

#Gets variance of each cell
var_table_KC <- exp_KC_table * (1 - row_proportions_KC) %o% (1 - column_proportions_KC)

#Gets standard deviation of each cell
sd_table_KC <- sqrt(var_table_KC)

#Gets adjusted residuals of each cell
adj_resid_KC <- difference_KC / sd_table_KC

#Calculate p-value of each cell (before multiple test correction)
#First calculate the probability of each cell being any value up to the positive of the actual value in the cell
probability_up_to_value_KC <- pnorm(abs(adj_resid_KC))

#Calculate the probability of the cell being at least that extreme on BOTH sides (un adjusted p-values)
unadj_pvals_KC <- 2 * (1 - probability_up_to_value_KC)

#Now perform multiple test correction (FDR)
#FDR ranks the p-values from smallest to largest and changes the threshold for significance based upon how large
#the p-value is compared to the others (larger p-value means larger threshold)
#Turn unadjusted p-values into vector form
unadj_pvals_vec_KC <- as.vector(unadj_pvals_KC)

#Sort from smallest to largest
order_index_KC <- order(unadj_pvals_vec_KC)
sorted_unadj_p_KC <- unadj_pvals_vec_KC[order_index_KC]

#Total number of tests
total_KC <- length(unadj_pvals_vec_KC)

#Compute for each sorted p value the BH critical value
q_vals_KC <- sorted_unadj_p_KC * total_KC / seq_along(sorted_unadj_p_KC)

#Calculate cumulative min going backwards
sorted_adj_KC <- rev(cummin(rev(q_vals_KC)))

#Put adjusted values back in order
adj_pvals_KC <- numeric(total_KC)
adj_pvals_KC[order_index_KC] <- sorted_adj_KC

#Convert back into table
adj_pvals_table_KC <- matrix(adj_pvals_KC, nrow = nrow(obs_KC_table), ncol = ncol(obs_KC_table), dimnames = dimnames(obs_KC_table))
print("P-values For Cluster-Complications Cells:")
print(adj_pvals_table_KC)

#Print the significant cells
sig_cells_KC <- which(adj_pvals_table_KC < 0.05, arr.ind = TRUE)
if (length(sig_cells_KC) == 0) {
  print("No significant p-values (< 0.05)")
} else {
  sig_results_KC <- data.frame(
    Cluster = rownames(adj_pvals_table_KC)[sig_cells_KC[, 1]],
    Complication = colnames(adj_pvals_table_KC)[sig_cells_KC[, 2]],
    Adj_Pval = adj_pvals_table_KC[sig_cells_KC],
    Ratio_OE = ratios_KC[sig_cells_KC]
  )
  print("Significant Cluster-Complication Associations (p < 0.05):")
  print(sig_results_KC)
}

#Print O/E ratios
print("Observed/Expected (O/E) Ratios:")
print(ratios_KC)



#FOR KS PAIRS ONLY



#Collapse over complication and surgery to get just observed K x S table
obs_KS_table <- table(data$Cluster, data$Surgery)

#Get expected KS_table
chi_result_KS <- chisq.test(obs_KS_table)
exp_KS_table <- chi_result_KS$expected

#Gets the chi-square statistic, the degrees of freedom and the p-value between observed and expected values
chi_square_KS <- sum((obs_KS_table - exp_KS_table)^2 / exp_KS_table)
print(paste("Pearson Chi-square statistic for cluster and diagnosis:", chi_square_KS))

df_KS = prod(dim(obs_KS_table) - 1)
print(paste("Degrees of freedom for cluster and surgery:", df_KS))

p_value_KS <- 1 - pchisq(chi_square_KS, df_KS)
print(paste("P-value for whether cluster and surgery are independent overall:", p_value_KS))

#Computes the ratios observed/expected and ln(ratios) to find how much more the cluster-surgery pair occurs
ratios_KS <- obs_KS_table / exp_KS_table
log_ratios_KS <- log(ratios_KS)

#Now get the table of the standard deviations of each cell
#Finds the difference in observed and expected values
difference_KS <- obs_KS_table - exp_KS_table

#Computes row and column proportions of the table
row_proportions_KS <- rowSums(obs_KS_table) / sum(obs_KS_table)
column_proportions_KS <- colSums(obs_KS_table) / sum(obs_KS_table)

#Gets variance of each cell
var_table_KS <- exp_KS_table * (1 - row_proportions_KS) %o% (1 - column_proportions_KS)

#Gets standard deviation of each cell
sd_table_KS <- sqrt(var_table_KS)

#Gets adjusted residuals of each cell
adj_resid_KS <- difference_KS / sd_table_KS

#Calculate p-value of each cell (before multiple test correction)
#First calculate the probability of each cell being any value up to the positive of the actual value in the cell
probability_up_to_value_KS <- pnorm(abs(adj_resid_KS))

#Calculate the probability of the cell being at least that extreme on BOTH sides (un adjusted p-values)
unadj_pvals_KS <- 2 * (1 - probability_up_to_value_KS)

#Now perform multiple test correction (FDR)
#FDR ranks the p-values from smallest to largest and changes the threshold for significance based upon how large
#the p-value is compared to the others (larger p-value means larger threshold)
#Turn unadjusted p-values into vector form
unadj_pvals_vec_KS <- as.vector(unadj_pvals_KS)

#Sort from smallest to largest
order_index_KS <- order(unadj_pvals_vec_KS)
sorted_unadj_p_KS <- unadj_pvals_vec_KC[order_index_KS]

#Total number of tests
total_KS <- length(unadj_pvals_vec_KS)

#Compute for each sorted p value the BH critical value
q_vals_KS <- sorted_unadj_p_KS * total_KS / seq_along(sorted_unadj_p_KS)

#Calculate cumulative min going backwards
sorted_adj_KS <- rev(cummin(rev(q_vals_KS)))

#Put adjusted values back in order
adj_pvals_KS <- numeric(total_KS)
adj_pvals_KS[order_index_KS] <- sorted_adj_KS

#Convert back into table
adj_pvals_table_KS <- matrix(adj_pvals_KS, nrow = nrow(obs_KS_table), ncol = ncol(obs_KS_table), dimnames = dimnames(obs_KS_table))
print("P-values For Cluster-Surgery Cells:")
print(adj_pvals_table_KS)

#Print the significant cells
sig_cells_KS <- which(adj_pvals_table_KS < 0.05, arr.ind = TRUE)
if (length(sig_cells_KS) == 0) {
  print("No significant p-values (< 0.05)")
} else {
  sig_results_KS <- data.frame(
    Cluster = rownames(adj_pvals_table_KS)[sig_cells_KS[, 1]],
    Surgery = colnames(adj_pvals_table_KS)[sig_cells_KS[, 2]],
    Adj_Pval = adj_pvals_table_KS[sig_cells_KS],
    Ratio_OE = ratios_KS[sig_cells_KS]
  )
  print("Significant Cluster-Surgery Associations (p < 0.05):")
  print(sig_results_KS)
}

#Print O/E ratios
print("Observed/Expected (O/E) Ratios:")
print(ratios_KS)