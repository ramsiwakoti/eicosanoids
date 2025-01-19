# Load necessary libraries
library(qgcomp)
library(dplyr)

# Define the function for quantile gcomp model
quantile_gcomp_model <- function(data, exposures, outcomes, covariates, weight = NULL, q = 4) {
  final_results <- list()
  pos_weights_list <- list()
  neg_weights_list <- list()
  
  # Define the formula for exposures and covariates
  exposure_formula <- paste(exposures, collapse = " + ")
  covariate_formula <- paste(covariates, collapse = " + ")
  
  # Loop over each outcome
  for (outcome in outcomes) {
    # Define the full formula
    formula <- as.formula(paste("log(", outcome, ") ~ ", exposure_formula, " + ", covariate_formula, sep = ""))
    
    # Fit the quantile gcomp model
    if (!is.null(weight)) {
      model <- qgcomp.glm.noboot(
        formula, 
        expnms = exposures,
        data = data, 
        weights = data[[weight]], 
        family = gaussian(), 
        q = q
      )
    } else {
      model <- qgcomp.glm.noboot(
        formula, 
        expnms = exposures,
        data = data, 
        family = gaussian(), 
        q = q
      )
    }
    
    # Extract results
    model_summary <- summary(model)
    
    # Store model summary
    model_summary_tab <- as.data.frame(model_summary$coefficients)
    model_summary_tab$term <- rownames(model_summary_tab)
    model_summary_tab$outcome <- outcome
    model_summary_tab2 <- model_summary_tab %>% filter(term == "psi1")
    
    # Extract and store positive weights
    if (length(model$pos.weights) > 0) {
      pos_weight_df <- data.frame(term = names(model$pos.weights), pos_weight = model$pos.weights, 
                                  outcome = outcome, weight = "pos", stringsAsFactors = FALSE)
      names(pos_weight_df)  <- c("term", "weight", "outcome", "weightdir")
      pos_weights_list[[outcome]] <- pos_weight_df
    }
    
    # Extract and store negative weights
    if (length(model$neg.weights) > 0) {
      neg_weight_df <- data.frame(term = names(model$neg.weights), neg_weight = model$neg.weights, 
                                  outcome = outcome, weight = "neg", stringsAsFactors = FALSE)
      names(neg_weight_df)  <- c("term", "weight", "outcome", "weightdir")
      
      neg_weights_list[[outcome]] <- neg_weight_df
    }
    
    # Append the results
    final_results[[outcome]] <- model_summary_tab2
  }
  
  # Combine results into single data frames
  final_results_df <- bind_rows(final_results)
  pos_weights_df <- bind_rows(pos_weights_list)
  neg_weights_df <- bind_rows(neg_weights_list)
  
  return(list(final_results = final_results_df, pos_weights = pos_weights_df, neg_weights = neg_weights_df))
}

# Apply to LC
results_qgcomp_lc <- quantile_gcomp_model(data = pfas_eic_lc, exposures = exposures_lc, 
                                          outcomes = eic_all2, covariates = covars_list_lc, q = 4)

# Print the results
qgcomp_lc <- results_qgcomp_lc$final_results
rownames(qgcomp_lc) <- NULL
poswt_lc <- results_qgcomp_lc$pos_weights
rownames(poswt_lc) <- NULL
negwt_lc <- results_qgcomp_lc$neg_weights
rownames(negwt_lc) <- NULL

posnegwt_lc <- rbind.data.frame(poswt_lc, negwt_lc) %>% mutate(vars = outcome) %>% inner_join(eic_class, by = "vars") %>%
  arrange(group)

#Only get weights for outcomes with significant associations
posnegwt_lc_sig <- qgcomp_lc %>% filter(`Pr(>|t|)` <= 0.1) %>% inner_join(posnegwt_lc, by = "outcome")

#now plot

results_qgcomp_lc2 <- qgcomp_lc %>% mutate(vars = outcome) %>% inner_join(eic_class, by = "vars") %>%
  arrange(group) 


ggplot(results_qgcomp_lc2, aes(x = Estimate, xmin = `Lower CI`, xmax = `Upper CI`, y = outcome, col = group)) +
  geom_pointrange(position = position_dodge(width = 0.5), fatten = 1.5) +
  geom_vline(xintercept = 0) +
  facet_wrap(~group, scales = "free_y", nrow = 4)

#Plot the weights

#Plot weights for significant outcomes
posnegwt_lc_sig %>% mutate(weight2 = ifelse(weightdir == "neg", -1*weight, weight)) %>%
  ggplot(aes(y = weight2, x= term.y, fill = weightdir)) +
  geom_bar(stat = "identity") +
  facet_wrap(~outcome, scales = "free")


#Apply to PR
results_qgcomp_pr <- quantile_gcomp_model(data = pfas_eic_pr, exposures = exposures_short, 
                                          outcomes = eic_all2, covariates = covars_list_pr, q = 4)

# Print the results
qgcomp_pr <- results_qgcomp_pr$final_results
rownames(qgcomp_pr) <- NULL
poswt_pr <- results_qgcomp_pr$pos_weights
rownames(poswt_pr) <- NULL
negwt_pr <- results_qgcomp_pr$neg_weights
rownames(negwt_pr) <- NULL


posnegwt_pr <- rbind.data.frame(poswt_pr, negwt_pr) %>%  mutate(vars = outcome) %>% inner_join(eic_class, by = "vars") %>%
  arrange(group)


#Only get weights for outcomes with significant associations
posnegwt_pr_sig <- qgcomp_pr %>% filter(`Pr(>|t|)` <= 0.1) %>% inner_join(posnegwt_pr, by = "outcome")


#now plot

results_qgcomp_pr2 <- qgcomp_pr %>% mutate(vars = outcome) %>% inner_join(eic_class, by = "vars") %>%
  arrange(group)


ggplot(results_qgcomp_pr2, aes(x = Estimate, xmin = `Lower CI`, xmax = `Upper CI`, y = outcome, col = group)) +
  geom_pointrange(position = position_dodge(width = 0.5), fatten = 1.5) +
  geom_vline(xintercept = 0) +
  facet_wrap(~group, scales = "free_y", nrow = 4)

#Plot weights for significant outcomes
posnegwt_pr_sig %>% mutate(weight2 = ifelse(weightdir == "neg", -1*weight, weight)) %>%
  ggplot(aes(y = weight2, x= term.y, fill = weightdir)) +
  geom_bar(stat = "identity") +
  facet_wrap(~outcome, scales = "free")





