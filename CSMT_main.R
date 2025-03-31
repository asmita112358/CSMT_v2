##CSMT main function, Mar 29
library(parallel)
library(caret)

divide_into_subsamples <- function(data, K) {
  n <- nrow(data)  # Total number of samples
  cat_vars <- names(data)[sapply(data, is.factor)]
  
  if (length(cat_vars) == 0) {
    ref_var <- data$predictor
  } else if (length(cat_vars) == 1) {
    ref_var <- data[[cat_vars]]  # Use the single categorical variable
  } else {
    interaction_var <- interaction(data[, cat_vars], sep = "_")
    ref_var <- interaction_var  # Use interaction of categorical variables
  }
  
  folds <- caret::createFolds(ref_var, k = K, list = TRUE)
  
  subsamples <- lapply(folds, function(idx) data[idx, , drop = FALSE])
  return(subsamples)
}

CSMT_pval <- function(model1, model2, data,intercept = FALSE, 
                      family1 = gaussian(), family2 = gaussian(),
                      K = NULL, num.subsamples = 500, wts = "uniform") {
  
  data <- na.omit(data)
  n <- nrow(data)
  # Define the number of partitions
  if (is.null(K)) {
    K <- floor(0.5 * sqrt(n))
  }
  
  min_size <- floor(n / K)
  single_iteration <- function(iter) {
    # Divide data into K stratified subsamples
    subsamples <- divide_into_subsamples(data = data, K)
    
    S_K <- numeric(K)
    
    for (i in seq_len(K)) {
      data_s <- subsamples[[i]]
      
      # Fit model1 (mediator ~ predictor + other_covariates)
      fit1 <- glm(model1, data = data_s, family = family1)
      
      # Fit model2 (outcome ~ mediator + predictor + other_covariates)
      fit2 <- glm(model2, data = data_s, family = family2)
      
      # Extract t-statistics of predictor coefficients
      if(intercept == FALSE){
        t1 <- summary(fit1)$coefficients[1,3]
        t2 <- summary(fit2)$coefficients[1,3]
      }else{
        t1 <- summary(fit1)$coefficients[2,3]
        t2 <- summary(fit2)$coefficients[2,3]
      }
      
      
      # Compute S_K statistic
      S_K[i] <- (t1 * t2) / sqrt(t1^2 + t2^2)
    }
    
    # Compute t-statistic and p-value
    tstat_iter <- sqrt(K) * mean(S_K) / sqrt(var(S_K))
    pval_iter <- 2 * pt(-abs(tstat_iter), df = K - 1)
    
    return(list(tstat = tstat_iter, pval = pval_iter))
  }
}


