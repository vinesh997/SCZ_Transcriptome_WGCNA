library(ppcor)
library(ggplot2)
library(reshape2)
library(tibble)
library(dplyr)
library(tidyr)

setwd("C:/Vinesh/PhD_Analysis/RNA-seq/WGCNA_23_10_2024/Correlation")


Data<- read.csv("PANSS_Profile_TScore.csv", check.names = F) #dim 19;25
Data<- column_to_rownames(Data, var = "Scale")
Data<- as.data.frame(t(Data))

# Define variable pairs (all combinations of Pair 1 with Pair 2)
PANSS <- c("Positive", "Negative", "General", "Anergia", 
            "Thought Disturbance", "Activation", "Paranoid/Belligerence", 
            "Depression", "Composite", "Total")

Hemogram <- c("Lymphocytes", "Neutrophils", 
                "Monocytes", "Eosinophils", "Basophils" , "Platelet", "Haemoglobin", 
                "RBC count", "Tobacco_consumption", 
              "Alcohol_consumption", "Cannabis_consumption")

var_pairs <- expand.grid(PANSS, Hemogram)
colnames(var_pairs) <- c("PANSS", "Hemogram")

# Predefine the cofactors
cofactors <- c("Sex", "Age", "Duration")

# Ensure variable names are character strings
var_pairs$PANSS <- as.character(var_pairs$PANSS)
var_pairs$Hemogram <- as.character(var_pairs$Hemogram)


# Initialize an empty dataframe to store results
partial_corr_results <- data.frame(PANSS = character(), 
                                   Hemogram = character(), 
                                   PartialCorrelation = numeric(), 
                                   PValue = numeric(), 
                                   Statistic = numeric(), 
                                   N = numeric(), 
                                   GP = numeric(), 
                                   Method = character(), 
                                   stringsAsFactors = FALSE)

# Loop through the variable pairs
for (i in 1:nrow(var_pairs)) {
  x <- var_pairs[i, "PANSS"]
  y <- var_pairs[i, "Hemogram"]
  
  cat("\nIteration:", i, 
      "| PANSS:", x, 
      "| Hemogram:", y, "\n")
  
  # Perform partial correlation inside a tryCatch block to handle potential errors
  result <- tryCatch({
    ppcor::pcor.test(
      x = Data[[x]],
      y = Data[[y]],
      z = Data[, cofactors],
      method = "spearman"
    )
  }, error = function(e) {
    cat("Error in iteration", i, ":", e$message, "\n")
    return(NULL)  # Return NULL in case of an error
  })
  
  # Skip this iteration if result is NULL (i.e., an error occurred)
  if (is.null(result)) next
  
  # Ensure all necessary components exist before storing results
  method_value <- if (!is.null(result$Method)) result$Method else "Unknown"
  
  # Print the results for this iteration
  cat("Estimate:", result$estimate, 
      "| P-Value:", result$p.value, 
      "| Statistic:", result$statistic, 
      "| N:", result$n, 
      "| GP:", result$gp, 
      "| Method:", method_value, "\n")
  
  # Append the results
  partial_corr_results <- rbind(partial_corr_results, data.frame(
    PANSS = x,
    Hemogram = y,
    PartialCorrelation = result$estimate,
    PValue = result$p.value,
    Statistic = result$statistic,
    N = result$n,
    GP = result$gp,
    Method = method_value,  # Use default value if missing
    stringsAsFactors = FALSE
  ))
}



# Reshape data for plotting (correlation matrix style)
correlation_matrix <- dcast(partial_corr_results, PANSS ~ Hemogram, value.var = "PartialCorrelation")

# Replace NAs with zeros
correlation_matrix[is.na(correlation_matrix)] <- 0

# Convert to a long format for ggplot2
melted_corr <- melt(correlation_matrix, id.vars = "PANSS", variable.name = "Hemogram", value.name = "PartialCorrelation")

# Add labels for partial correlation and p-value
melted_corr <- merge(
  melted_corr, 
  partial_corr_results[, c("PANSS", "Hemogram", "PValue")],
  by = c("PANSS", "Hemogram")
)

# Format labels to include both correlation and p-value
melted_corr$Label <- paste0(
  "r = ", round(melted_corr$PartialCorrelation, 2), 
  "\n(p = ", formatC(melted_corr$PValue), ")"
)

# Create the correlation heatmap with values and p-values
ggplot(melted_corr, aes(x = PANSS, y = Hemogram, fill = PartialCorrelation)) +
  geom_tile(color = "white") +
  geom_text(aes(label = Label), color = "black", size = 3) +  # Add correlation and p-value
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Partial Corr") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "PANSS syndrome and clusters", y = "Hemogram count")


#verification of individual correlation
ppcor::pcor.test(x=Data$Negative,y=Data$Platelet, z=Data[,c( "Sex", "Age", "Duration")], 
                                        method = "spearman")

#Filtering the associations based on SCZ and HC mean count

#reading the hemogram count of SCZ and HC participants
hem_scz_hc<- read.csv("Hemogram_SCZ_HC.csv", check.names = F)

#calculating mean based on groups 
group_means <- hem_scz_hc %>%
  group_by(Groups) %>%
  summarise(across(where(is.numeric), mean, na.rm = TRUE)) %>%
  pivot_longer(-Groups, names_to = "Hemogram", values_to = "Mean") %>%
  pivot_wider(names_from = Groups, values_from = Mean, names_prefix = "Mean_")


group_means$Hemogram <- recode(group_means$Hemogram,
                               "RBC" = "RBC count", 
                               "Platelets" = "Platelet", 
                               "Haemoglobin" = "Haemoglobin")


melted_corr <- merge(melted_corr, group_means, by = "Hemogram")


melted_corr$DirectionalMatch <- with(melted_corr,
                                     ifelse(
                                       (PartialCorrelation > 0 & Mean_SCZ > Mean_CNT) |
                                         (PartialCorrelation < 0 & Mean_SCZ < Mean_CNT),
                                       TRUE, FALSE
                                     )
)

ggplot(melted_corr, aes(x = PANSS, y = Hemogram, fill = PartialCorrelation)) +
  geom_tile(color = "white") +
  geom_text(
    aes(label = ifelse(PValue < 0.05 & DirectionalMatch, Label, "")),
    color = "black", size = 3
  ) +
  scale_fill_gradient2(
    low = "blue", high = "red", mid = "white", midpoint = 0,
    limit = c(-1, 1), name = "Partial Corr"
  ) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "PANSS syndrome and clusters", y = "Hemogram count")


