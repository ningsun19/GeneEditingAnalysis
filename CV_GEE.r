# Example in R for aggregated data (y1, y0 format)
library(stats)

# data has columns: prior_base, after_base, position, y1, y0, cell_type

# Unique cell types
cell_types <- unique(m_wide$cell)

# Initialize lists to store performance metrics
deviance_results <- list()
brie_results <- list()

for (type in cell_types) {
  # Split the data
  training_set <- subset(m_wide, cell != type)
  test_set <- subset(m_wide, cell == type)
  
  # Fit model on training data
  # Using weights for the binomial counts
  model <- geepack::geeglm(cbind(y1, y0) ~ m1 + p1 + group, id = as.factor(cell) , corstr ="exchangeable" , family = "binomial", data = training_set)
  
  # Predict on test set
  test_set$predicted_prob <- predict(model, newdata = test_set, type = "response")
  
  # Calculate model fit (deviance) on test set
  # Clamp probabilities to avoid log(0)
  eps   <- 1e-15
  p_hat <- pmin(pmax(test_set$predicted_prob, eps), 1 - eps)
  
  # ---- Null model: estimate from TRAINING data ----
  train_prob <- with(training_set, sum(y1) / sum(y1 + y0))  # overall event rate
  p_null     <- pmin(pmax(train_prob, eps), 1 - eps)
  
  # ---- Log-likelihoods on TEST data ----
  model_ll <- with(test_set, sum(y1 * log(p_hat) + y0 * log(1 - p_hat)))
  null_ll  <- with(test_set, sum(y1 * log(p_null) + y0 * log(1 - p_null)))
  
  # "Deviance explained" style metric based on log-likelihood
  deviance_explained <- 1 - (model_ll / null_ll)
  deviance_results[[type]] <- deviance_explained
  
  # Weighted Brier score on aggregated data
  obs_prob <- with(test_set, y1 / (y1 + y0))
  w        <- with(test_set, y1 + y0)
  
  weighted_brier <- sum((p_hat - obs_prob)^2 * w) / sum(w)
  brie_results[[type]] <- weighted_brier
}

# Calculate average performance across folds
mean_deviance_explained <- mean(unlist(deviance_results))
sd_deviance_explained <- sd(unlist(deviance_results))
print(paste("Cross-validated deviance explained:", mean_deviance_explained))


model <- geepack::geeglm(cbind(y1, y0) ~ m1 + p1 + group, id = as.factor(cell) , corstr ="exchangeable" , family = "binomial", data = m_wide)
m_wide$phat_train <- predict(model, type = "response")

eps <- 1e-15
p_hat <- pmin(pmax(m_wide$phat_train, eps), 1 - eps)

# Null probability from the same training data
p_null <- with(m_wide, sum(y1) / sum(y1 + y0))
p_null <- pmin(pmax(p_null, eps), 1 - eps)

ll_model <- with(m_wide, sum(y1 * log(p_hat) + y0 * log(1 - p_hat)))
ll_null  <- with(m_wide, sum(y1 * log(p_null) + y0 * log(1 - p_null)))

dev_exp_train <- 1 - (ll_model / ll_null)
dev_exp_train

mean(unlist(brie_results))
