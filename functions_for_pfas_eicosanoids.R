
#PFAS summary

# Function to calculate summary stats for a single column
pfas_col_summary <- function(data, col, cohort_col = "cohort") {
  # Create corresponding BLOD column name
  blod_col <- gsub("^conc_", "blod_", col)
  
  # Calculate summary statistics
  summary <- data %>%
    group_by(!!sym(cohort_col)) %>%
    summarise(
      N = sum(!is.na(.data[[col]])),
      Detection_Rate = 100 - mean(.data[[blod_col]], na.rm = TRUE)*100,
      Median = median(.data[[col]], na.rm = TRUE),
      Q1 = quantile(.data[[col]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(PFAS = sub("^conc_", "", col)) %>%
    mutate_if(is.numeric, round, 2)
  
  return(summary)
}

#Eicosanoids summary

# Function to calculate summary stats for a single column
eic_col_summary <- function(data, col, cohort_col = "cohort") {
 
  # Calculate summary statistics
  summary <- data %>%
    group_by(!!sym(cohort_col)) %>%
    summarise(
      Median = median(.data[[col]], na.rm = TRUE),
      Q1 = quantile(.data[[col]], 0.25, na.rm = TRUE),
      Q3 = quantile(.data[[col]], 0.75, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    mutate(Eicosanoids = col) %>%
    mutate_if(is.numeric, round, 2)
  
  return(summary)
}









#Creating correlation plots

create_heatplot <- function(data, variables, cohort_variable) {
  # Ensure the cohort variable is in the data
  if (!cohort_variable %in% colnames(data)) {
    
    # Compute correlation matrix
    cor_matrix <- cor(data[, variables], method = "spearman", use = "pairwise.complete.obs")
    
    # Convert correlation matrix to long format for ggplot
    cor_long <- as.data.frame(as.table(cor_matrix))
    colnames(cor_long) <- c("Var1", "Var2", "Correlation")
    # cor_long$Cohort <- "All"

    # Append to heatplot data
    heatplot_data_combined <- cor_long
    
    
    # Create heatplot
    ggplot(heatplot_data_combined, aes(x = Var1, y = Var2, fill = Correlation)) +
      geom_tile(color = "white") +
      geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
      scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
      labs(title = "Correlation Heatplot for all", fill = "Correlation") +
      xlab("") + ylab("")
  }
  
  else{
  
  # Unique cohorts
  unique_cohorts <- unique(data[[cohort_variable]])
  
  # List to store heatplot data
  heatplot_data <- list()
  
  # Loop through each cohort and compute correlation
  for (c in unique_cohorts) {
    # Subset data for the cohort
    subset_data <- data %>% filter(.data[[cohort_variable]] == c)
    
    # Compute correlation matrix
    cor_matrix <- cor(subset_data[, variables], method = "spearman", use = "pairwise.complete.obs")
    
    # Convert correlation matrix to long format for ggplot
    cor_long <- as.data.frame(as.table(cor_matrix))
    colnames(cor_long) <- c("Var1", "Var2", "Correlation")
    cor_long$Cohort <- c
    
    # Append to heatplot data
    heatplot_data[[c]] <- cor_long
  }
  
  # Combine all cohort data
  heatplot_data_combined <- do.call(rbind, heatplot_data)


  # Create heatplot
  ggplot(heatplot_data_combined, aes(x = Var1, y = Var2, fill = Correlation)) +
    geom_tile(color = "white") +
    geom_text(aes(label = round(Correlation, 2)), color = "black", size = 3) +
    facet_wrap(~ Cohort, scales = "free") +
    scale_fill_gradient2(low = "blue", high = "red", mid = "white", midpoint = 0, limit = c(-1, 1)) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Correlation Heatplot by Cohort", fill = "Correlation") +
    xlab("") + ylab("")
}  }

# Example usage
# Assuming pfas_eic_combined is your dataset
# variables_short is the vector of variables you want to correlate
# "cohort" is the column indicating cohorts

# create_heatplot(pfas_eic_combined, exposures_short, "cohort")
# create_heatplot(pfas_eic_combined, exposures_short, "")
# 
# 
# #Now create correlation plot between eicosanoids 
# 
# create_heatplot(pfas_eic_combined, eic_all, "")

#For eicosanoids with group

cor_eic_func <- function(data, cohort_variable){
  
  if (cohort_variable == "All"){
    subset_data <- data}
  else {
    subset_data <- data %>% filter(cohort == cohort_variable)
  }
  
  cor_eic <- cor(subset_data[, eic_all], method = "spearman", use = "pairwise.complete")
  
  # Reorder correlation matrix to match sorted annotation
  cor_eic_ordered <- cor_eic[eic_class$vars, eic_class$vars]
  
  p <- pheatmap(cor_eic_ordered,
                annotation_row = annotation_row_label,
                annotation_col = annotation_row_label,
                display_numbers = TRUE,
                number_format = "%.2f",
                show_rownames = TRUE, 
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                annotation_legend = TRUE,
                annotation_names_row = TRUE,
                fontsize = 6,
                main = ""
  ) 
  out <- list(p, cor_eic_ordered)
  return(out)
}

#filtering based on threshold
cor_eic_func2 <- function(data, cohort_variable, threshold = 0.3) {
  
  # Subset data based on cohort_variable
  if (cohort_variable == "All") {
    subset_data <- data
  } else {
    subset_data <- data %>% filter(cohort == cohort_variable)
  }
  
  # Calculate correlation matrix
  cor_eic <- cor(subset_data[, eic_all], method = "spearman", use = "pairwise.complete")
  
  # Reorder correlation matrix to match sorted annotation
  cor_eic_ordered <- cor_eic[eic_class$vars, eic_class$vars]
  
  # Create a matrix of numbers to display, applying the threshold
  display_numbers <- ifelse(abs(cor_eic_ordered) >= threshold, 
                            sprintf("%.2f", cor_eic_ordered), 
                            "")
  
  # Create heatmap using filtered correlations
  p <- pheatmap(cor_eic_ordered,
                annotation_row = annotation_row_label,
                annotation_col = annotation_row_label,
                display_numbers = display_numbers,  # Use the custom display matrix
                number_format = "%.2f",
                show_rownames = TRUE, 
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                annotation_legend = TRUE,
                annotation_names_row = TRUE,
                fontsize = 8,
                main = paste(cohort_variable, "(Filtered, >", threshold, ")")
  ) 
  p 
  
  return(cor_eic_ordered)
}




#filtering based on threshold
cor_eic_pfas_func <- function(data, cohort_variable, pfas_vec) {
  
  # Subset data based on cohort_variable
  if (cohort_variable == "All") {
    subset_data <- data
  } else {
    subset_data <- data %>% filter(cohort == cohort_variable)
  }
  
  # Calculate correlation matrix
  cor_eic <- cor(log(subset_data[, eic_all]), log(subset_data[, pfas_vec]), 
                 method = "spearman", use = "pairwise.complete")
  print(cor_eic)
  
  grouping_info <- eic_class %>% arrange(group, vars)  # Assuming `Group` is the grouping column in eic_class
  reordered_vars <- grouping_info$vars
  # 
  # Subset and reorder correlation matrix rows based on grouping
  cor_eic_ordered <- cor_eic[reordered_vars, ]
  

  
  # Create heatmap using filtered correlations
  p <- pheatmap(cor_eic_ordered,
                annotation_row = annotation_row_label,
                # display_numbers = display_numbers,  # Use the custom display matrix
                display_numbers = TRUE,
                number_format = "%.2f",
                show_rownames = TRUE, 
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                annotation_legend = TRUE,
                annotation_names_row = FALSE,
                fontsize = 8,
                main = ""
  ) 
  output <- list(p, cor_eic_ordered) 
  
  return(output)
}

#Correlation between eicosanoids and OS


#filtering based on threshold
cor_eic_osinf_func <- function(data, cohort_variable, osinf) {
  
  # Subset data based on cohort_variable
  if (cohort_variable == "All") {
    subset_data <- data
  } else {
    subset_data <- data %>% filter(cohort == cohort_variable)
  }
  
  # Calculate correlation matrix
  cor_eic <- cor(log(subset_data[, eic_all]), log(subset_data[, osinf]), 
                 method = "spearman", use = "pairwise.complete")
  print(cor_eic)
  
  grouping_info <- eic_class %>% arrange(group, vars)  # Assuming `Group` is the grouping column in eic_class
  reordered_vars <- grouping_info$vars
  # 
  # Subset and reorder correlation matrix rows based on grouping
  cor_eic_ordered <- cor_eic[reordered_vars, ]
  
  
  
  # Create heatmap using filtered correlations
  p <- pheatmap(cor_eic_ordered,
                annotation_row = annotation_row_label,
                # display_numbers = display_numbers,  # Use the custom display matrix
                display_numbers = TRUE,
                number_format = "%.2f",
                show_rownames = TRUE, 
                cluster_rows = FALSE,
                cluster_cols = FALSE,
                annotation_legend = TRUE,
                annotation_names_row = FALSE,
                fontsize = 8,
                main = paste(cohort_variable)
  ) 
  p 
  
  return(cor_eic_ordered)
}


#Adjusted models for PFAS and eicosanoids
pfas_eic_adjmod <- function(data, exposures, outcomes, covariates, weg = 1) {
  library(broom)
  library(dplyr)
  
  results <- list()
  
  for (exposure in exposures) {
    for (outcome in outcomes) {
      # Define the formula
      
      if (substr(exposure,1,4) == "alod"){
        
      formula <- as.formula(
      paste("log(", outcome, ") ~",  exposure, "+ ", paste(covariates, collapse = " + ")))  
        
      }
      else{
      
      formula <- as.formula(
        paste("log(", outcome, ") ~ log(", exposure, ") + ", paste(covariates, collapse = " + ")))
      }
      
      # Fit the model

      if (weg == 1){
        model <- glm(formula, data = data, family = gaussian())
      }
      else{
        model <- glm(formula, data = data, family = gaussian(), weights = data[,weg])
      }
      
      # Use broom::tidy() to extract the summary
      model_tidy <- tidy(model, conf.int = TRUE)
      
      # Filter for the coefficient of the exposure
      exposure_row <- model_tidy %>%
        filter(term == paste0("log(", exposure, ")")|term == exposure)
      
      # Store results
      results[[paste(exposure, outcome, sep = "_")]] <- exposure_row %>%
        mutate(Exposure = exposure, Outcome = outcome) %>%
        select(Exposure, Outcome, Estimate = estimate, `95% CI Lower` = conf.low, `95% CI Upper` = conf.high, `P-value` = p.value)
    }
  }
  
  # Combine results into a single data frame
  final_results <- bind_rows(results)
  return(final_results)
}

#Function to chekc the effect of year

#want to check whether the effects of year are significant
pfas_eic_adjmod_year <- function(data, exposures, outcomes, covariates, weg = 1) {
  
  results <- list()
  
  for (exposure in exposures) {
    for (outcome in outcomes) {
      # Define the formula
      formula <- as.formula(
        paste("log(", outcome, ") ~ log(", exposure, ") + ", paste(covariates, collapse = " + "))
      )
      
      # Fit the model
      
      if (weg == 1){
        model <- glm(formula, data = data, family = gaussian())
      }
      else{
        model <- glm(formula, data = data, family = gaussian(), weights = data[,weg])
      }
      
      # Use broom::tidy() to extract the summary
      model_tidy <- tidy(model, conf.int = TRUE)
      
      # Filter for the coefficient of the exposure
      exposure_row <- model_tidy %>%
        filter(term == "year_t1")
      
      # Store results
      results[[paste(exposure, outcome, sep = "_")]] <- exposure_row %>%
        mutate(Exposure = exposure, Outcome = outcome) %>%
        select(Exposure, Outcome, Estimate = estimate, `95% CI Lower` = conf.low, `95% CI Upper` = conf.high, `P-value` = p.value)
    }
  }
  
  # Combine results into a single data frame
  final_results <- bind_rows(results)
  return(final_results)
}


#Plotting results from regression


plot_adj_result <- function(data, pfas){
  if (pfas == "conc_pfhxs"|pfas == "alod_pfhxs") {pfas2 = "PFHxS"}
  else if (pfas == "conc_pfoa"|pfas == "alod_pfoa") {pfas2 = "PFOA"}
  else if (pfas == "conc_pfua"|pfas == "alod_pfua") {pfas2 = "PFUA"}
  else {pfas2 = toupper(sub("conc_", "", pfas))}
  
  data %>%
    filter(exposure == pfas) %>%
    ggplot(aes(x = estimate, xmin = conf.low, xmax = conf.high, y = outcome_format, col = group)) +
    geom_pointrange() +
    theme(legend.position = "none") +
    geom_vline(xintercept = 0) +
    xlab(pfas2) +
    facet_grid(group ~ ., scales = "free_y", switch = "y") +  # Facet by group, switch facet labels to y-axis
    theme(strip.placement = "outside", strip.text.y = element_text(angle = 0),
          axis.title.y = element_blank()) }















