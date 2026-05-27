#Loading required libraries
library(ppcor)
library(ggplot2)
library(reshape2)
library(tibble)
library(dplyr)
library(tidyr)

#Setting working directory
setwd(" ")

#Reading immune cell proportion data
Data<- read.csv(" ", check.names = F)

#Filtering only SCZ samples
Data<- Data[Data$Diagnosis == "SCZ", ]

#Setting sample names as rownames
Data <- Data %>% remove_rownames() %>% column_to_rownames("Samples")

#Removing diagnosis column
Data<- Data[,-c(1)]

# Define variable pairs (all combinations of Pair 1 with Pair 2)

#PANSS symptom domains
PANSS <- c("P","N","G","T")

#Immune cell subtype proportions
Hemogram <- c("B cells naive", "B cells memory", "Plasma cells",	"T cells CD8",	
              "T cells CD4 naive", "T cells CD4 memory resting",	
              "T cells CD4 memory activated",	"T cells follicular helper",	"T cells regulatory (Tregs)",
              "T cells gamma delta", "NK cells resting",	"NK cells activated",	
              "Monocytes",	"Macrophages M0",	"Macrophages M1",	"Macrophages M2",
              "Dendritic cells resting",	"Dendritic cells activated")

#Creating all possible PANSS–immune cell combinations
var_pairs <- expand.grid(PANSS, Hemogram)
colnames(var_pairs) <- c("PANSS", "Hemogram")

# Predefine the cofactors
#Covariates used for partial correlation analysis
cofactors <- c("Sex", "Age")

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
  
  #Extract current PANSS and immune cell variables
  x <- var_pairs[i, "PANSS"]
  y <- var_pairs[i, "Hemogram"]
  
  #Print current iteration details
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
    
    #Print error message if analysis fails
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

#partial_corr_results<- mutate_all(partial_corr_results$PValue, as.numeric)

# Reshape data for plotting (correlation matrix style)
correlation_matrix <- dcast(partial_corr_results, PANSS ~ Hemogram, value.var = "PartialCorrelation")

# Replace NAs with zeros (if needed)
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
  geom_text(aes(label = Label), color = "black", size = 2) +  # Add correlation and p-value
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Partial Corr") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)) +
  labs(x = "PANSS syndrome", y = "Immune cell subtype")
