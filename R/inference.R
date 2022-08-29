# load packages
library(tidyverse)
library(MatchIt)

# calculate bootstrap confidence interval
bootstrap_ci <- function(match, dat, treat_variable = "treat",
                         outcome_variable = "y_obs",
                         trajectory_id = "personID",
                         ci_level = 0.95,
                         reps = 1000, outcome_model = NA){
  # add matching weights to data
  dat <- matrix_to_weights(match$cells, dat)
  # Sample only matched data
  matched_dat <- filter(dat, match_wt != 0)
  # get number of matched treated units
  num_matched_treated <- nrow(filter(matched_dat, !!sym(treat_variable) == 1))
  # initialize boostrap estimates
  boots_bias_correct <- rep(0, reps)
  # unique trajectory IDs
  matched_trajectory_ids <- unique(matched_dat[[trajectory_id]])
  # iterate over bootstrap replicates
  for(i in 1:reps){
    sample_trajectory_ids <- data.frame(trajectory_id = sample(matched_trajectory_ids, 
                                                          length(matched_trajectory_ids), 
                                                          replace = T))
    names(sample_trajectory_ids)[1] <- trajectory_id
    sample_dat <- left_join(sample_trajectory_ids, matched_dat, by = trajectory_id)
    # get bias correction values
    if(is.na(outcome_model) == F){
      preds <- predict(outcome_model, sample_dat)
    } else {
      preds <- rep(0, nrow(sample_dat))
    }
    # get bootstrap value
    boots_bias_correct[i] <- (1/num_matched_treated) * sum((2*sample_dat[[treat_variable]] - rep(1, nrow(sample_dat)))*(sample_dat$match_wt*(sample_dat[[outcome_variable]] - preds)))
  }
  # get confidence interval
  alpha <- 1 - ci_level
  boot_bias_correct_ci <- quantile(boots_bias_correct, c(alpha / 2, 1 - alpha / 2))
  # return confidence interval
  return(boot_bias_correct_ci)
}

# run falsification test
control_test <- function(control_data, control_variables, outcome_model = NA, 
                         outcome_variable, time_variable, trajectory_id,
                         time_points, rand_reps = 100,
                         caliper = 0.2, reps = 1000){
  # filter down to just the time points
  control_data <- control_data %>%
    filter(!!sym(time_variable) %in% time_points)
  
  # apply bias correction, if applicable
  if(is.na(outcome_model) == F){
    control_data$bias_correct <- predict(outcome_model, control_data)
  } else {
    control_data$bias_correct <- 0
  }
  
  # create matching formula
  matching_formula <- paste("new_treat ~ ", paste(control_variables, collapse = ' + '))
  
  # randomly split the data many times and take split with best balance
  balance_value <- 100000000
  for(i in 1:rand_reps){
    # get control timepoints
    new_control_ids <- sample(unique(control_data[[trajectory_id]]), 
                              length(unique(control_data[[trajectory_id]])) / 2, replace = F)
    # get treated timepoints
    new_treat_ids <- setdiff(unique(control_data[[trajectory_id]]), new_control_ids)
    
    # second timepoint is new control and first timepoint is new treatment
    new_control_df <- control_data %>% 
      filter(!!sym(trajectory_id) %in% new_control_ids, 
             !!sym(time_variable) == time_points[2]) %>%
      mutate(new_treat = 0)
    
    new_treat_df <- control_data %>% 
      filter(!!sym(trajectory_id) %in% new_treat_ids, 
             !!sym(time_variable) == time_points[1]) %>%
      mutate(new_treat = 1)
    
    # Calculate mahalonbis distance
    X1 <- new_control_df[control_variables]
    X2 <- new_treat_df[control_variables]
    X <- rbind(X1, X2)
    n <- nrow(X)
    p <- nrow(X2) / nrow(X)
    X_mean_diff <- colMeans(X2) - colMeans(X1)
    mahalanobis_value <- n * p * (1 - p) * t(X_mean_diff) %*% solve(cov(X)) %*% (X_mean_diff)
    
    # check if balanced improved
    if(mahalanobis_value < balance_value){
      # run match
      control_data_temp <- bind_rows(new_control_df, new_treat_df)
      mout <- matchit(formula(matching_formula), 
                           control_data_temp, replace = F, caliper = caliper)
      new_control_data <- control_data_temp
      balance_value <- mahalanobis_value
    }
  }

  control_data <- new_control_data
  
  matched_diffs <- (control_data[[outcome_variable]][control_data["new_treat"] == 1] - 
                      control_data$bias_correct[control_data["new_treat"] == 1]) - 
    (control_data[[outcome_variable]][as.numeric(mout$match.matrix)] - 
       control_data$bias_correct[as.numeric(mout$match.matrix)])

  matched_diffs <- matched_diffs[is.na(matched_diffs) == F]
  ATT.match <- mean(matched_diffs)
  
  # permutation test
  perm_coefs <- rep(0, reps)
  for(j in 1:reps){
    perm_coefs[j] <- mean(matched_diffs * sample(c(-1, 1), length(matched_diffs), replace = T))
  }
  # calculate p-value
  p_value <- mean(abs(perm_coefs) >= abs(ATT.match))
  
  return(list(p_value, balance_value))
}

