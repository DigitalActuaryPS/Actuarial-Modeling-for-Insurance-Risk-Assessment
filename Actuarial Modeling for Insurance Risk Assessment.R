#Reading entire dataset
Data <- read.csv("C:/Users/Margashi/OneDrive - City, University of London/Desktop/CS2/CS2_Coursework/CW2_Group_SMM048_23_24_Data.csv")

# Setting seed for reproducibility
set.seed(123)

# Extracting the primary data for Group 4 - dataset4
dataset_4 = Data[, c(1, 5)]

# Converting string data into numeric 
dataset_4$dataset4 <- as.numeric(gsub(",", "", dataset_4$dataset4))

#A1----------------------------------------------------------------------------->
# Non-parametric estimation of pure-premium based on the given estimator
# na.rm is used to remove any missing values or NaNs in the dataset

pure_prem_1 <- 0.65 * mean(dataset_4$dataset4, na.rm = TRUE) + 
  0.20 * sqrt(var(dataset_4$dataset4, na.rm = TRUE)) + 
  0.15 * quantile(dataset_4$dataset4, 0.9, na.rm = TRUE, type = 1)

# Removing the name associated with the quartile value
names(pure_prem_1) <- NULL

pure_prem_1

#A2----------------------------------------------------------------------------->
# Bootstrapping
# Assigning the primary claims data to claim_sizes
claim_sizes = c(dataset_4$dataset4) 

set.seed(123)
# Choosing number of bootstrap samples by trial and error, with the aim to minimize the MSE
num_bootstrap_samples = 9000

# Creating a numeric vector with 'num_bootstrap_samples' elements, all initialized to 0
pure_prem_results = numeric(num_bootstrap_samples)

# For loop to perform bootstrap sampling based on chosen number of bootstrap samples
for (i in 1:num_bootstrap_samples) {
  # Generating a single bootstrap sample by randomly sampling with replacement from the primary data
  bootstrap_sample = sample(claim_sizes, replace = TRUE)
  
  # Calculating the pure premium estimate for the bootstrap samples as per the given estimator
  pure_prem_results[i] = 0.65* mean(bootstrap_sample) + 0.20 *var(bootstrap_sample)**(0.5)+ 0.15*quantile(bootstrap_sample , 0.9, type=1) 
}


# Initalizing MSE to 0
MSE = 0

# Creating a numeric vector with 'num_bootstrap_samples' elements, all initialized to 0
a = numeric(num_bootstrap_samples)

# For loop to calculate MSE for the bootstrap samples
for (j in 1:num_bootstrap_samples) {
# The MSE of the pure premium given as MSE(ˆt) := 1/m * Σ_{k=1 to m} (tk - t)^2 formula is implemented below
a[j] = ((pure_prem_results[j] - pure_prem_1)^2)/num_bootstrap_samples
MSE = MSE + a[j]
}
# Using the MSE result for various number of bootstrap samples was narrowed down to 9000


#A3----------------------------------------------------------------------------->
# Installing and loading required packages
if (!require("MASS")) install.packages("MASS")
library(MASS)

# Assigning the primary claims data to claim_sizes
claim_sizes <- dataset_4$dataset4 

# Setting a list of distributions to be fitted from five parametric families via MLE and choosing the ‘best‘ possible fit
distribution_families <- c("normal", "log-normal", "weibull", "logistic", "exponential")

# Initializing variables for storing the best fit
best_fit_name <- NULL
best_fit_aic <- Inf
best_fit_model <- NULL

# Defining a function that fits a distribution to the data and calculates AIC
# Arguments:
#   - dist_name: Name of the distribution to fit (e.g., "normal", "exponential", etc.)
#   - data: The dataset for fitting the distribution

fit_distribution <- function(dist_name, data) {
  if (dist_name == "lognormal") {
    #If the distribution to be fitted is "lognormal" then log-transform data
    data <- log(data[data > 0])
  }
  
  # Attempting to fit the distribution
  fit <- try(fitdistr(data, dist_name), silent = TRUE)
  
  # Checking if the fit is a 'try-error'
  if(class(fit) == "try-error") {
    return(list(aic = Inf))
  }
  
  # Calculating AIC and return the fit and AIC
  aic <- AIC(fit)
  return(list(fit = fit, aic = aic))
}

# Initializing a data frame to store the fit results for each distribution
fit_details <- data.frame(Distribution = character(), AIC = numeric(), Parameters = character())

# For Loop to fit each distribution and calculate AIC over the five list of distributions
for (dist_name in distribution_families) {
  fit_result <- fit_distribution(dist_name, claim_sizes)
  
  # Check if we successfully got a fit object
  if(!is.na(fit_result$aic)) {
    # Extracting parameter estimates and their values
    param_estimates <- as.character(fit_result$fit$estimate)
    # Combining them into a single string
    params_str <- paste(names(fit_result$fit$estimate), param_estimates, sep = "=", collapse = "; ")
    # Save the fit details in the data frame
    fit_details <- rbind(fit_details, data.frame(Distribution = dist_name, AIC = fit_result$aic, Parameters = params_str))
  }
  
  # Checking if AIC is smaller than the current best AIC
  if (fit_result$aic < best_fit_aic) {
    # Updating the best AIC and the corresponding best-fit distribution
    best_fit_name <- dist_name
    best_fit_aic <- fit_result$aic
    best_fit_model <- fit_result$fit
  }
}

# Print the fit details for each distribution
print(fit_details)

# Output the best fit results
cat("Best-fit distribution:", best_fit_name, "\n")
cat("Best-fit AIC:", best_fit_aic, "\n")
print(best_fit_model)

# Initializing pure_premium variable
pure_premium <- NA

# Calculating the pure premium based on the best-fit distribution
if (!is.null(best_fit_name)) { # Check if best_fit_name is not NULL
  if (best_fit_name == "normal") {  # If the best-fit distribution is normal
    mean_normal <- best_fit_model$estimate["mean"] # Extracting mean parameter
    sd_normal <- best_fit_model$estimate["sd"]# Extracting standard deviation parameter
    # Calculating pure premium using weighted sum of mean, standard deviation, and 90th percentile
    pure_premium <- 0.65 * mean_normal + 0.20 * sd_normal + 0.15 * qnorm(0.9, mean = mean_normal, sd = sd_normal)
    
  } else if (!is.null(best_fit_model) && "meanlog" %in% names(best_fit_model$estimate) && "sdlog" %in% names(best_fit_model$estimate)) {
    # If the best-fit distribution is log-normal and the necessary parameters are present
    
    meanlog <- best_fit_model$estimate["meanlog"] # Extracting meanlog parameter
    sdlog <- best_fit_model$estimate["sdlog"]# Extracting sdlog parameter
    # Ensure the parameters are numeric
    if (is.numeric(meanlog) && is.numeric(sdlog)) {
  
# Computing mean, standard deviation, and 90th percentile for the log-normal distribution
      
    mean_lognormal <- exp(meanlog + (sdlog^2) / 2)
    sd_lognormal <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
    q90_lognormal <- qlnorm(0.9, meanlog = meanlog, sdlog = sdlog)
    
    pure_premium <- 0.65 * mean_lognormal + 0.20 * sd_lognormal + 0.15 * q90_lognormal
    }
    
  } else if (best_fit_name == "logistic") { # If the best-fit distribution is logistic
    location <- best_fit_model$estimate["location"] # Extracting location parameter
    scale <- best_fit_model$estimate["scale"] # Extracting scale parameter
    
    # Calculating mean, standard deviation, and 90th percentile for the logistic distribution
    mean_logistic <- location
    sd_logistic <- scale * pi / sqrt(3)
    q90_logistic <- qlogis(0.9, location = location, scale = scale)
    pure_premium <- 0.65 * mean_logistic + 0.20 * sd_logistic + 0.15 * q90_logistic
    
  } else if (!is.null(best_fit_model) && "shape" %in% names(best_fit_model$estimate) && "scale" %in% names(best_fit_model$estimate)) {
    # If the best-fit distribution is Weibull and the necessary parameters are present
    
    shape <- best_fit_model$estimate["shape"] # Extracting shape parameter
    scale <- best_fit_model$estimate["scale"] # Extracting scale parameter
    
    # Calculating mean, standard deviation, and 90th percentile for the Weibull distribution
    mean_weibull <- scale * gamma(1 + 1 / shape)
    sd_weibull <- sqrt(scale^2 * (gamma(1 + 2 / shape) - (gamma(1 + 1 / shape))^2))
    q90_weibull <- qweibull(0.9, shape = shape, scale = scale)
    pure_premium <- 0.65 * mean_weibull + 0.20 * sd_weibull + 0.15 * q90_weibull
    
  } else if (best_fit_name == "exponential") { # If the best-fit distribution is exponential
    rate <- best_fit_model$estimate["rate"] # Extracting rate parameter
    pure_premium <- 0.65 / rate + 0.20 * sqrt(1 / (rate^2)) + 0.15 * qexp(0.9, rate = rate)
  }
}

# Print the pure premium and verifying the if it is similar to the non-parametric estimate 
cat("Estimated Pure Premium based on the best-fit distribution:", pure_premium, "\n")


#B1----------------------------------------------------------------------------->
#Redoing above steps for the new dataset

# Extracting additional data for Group 4 - extradataset4
extra_dataset_4 <- Data[, c(1, 20)]

# Converting string data into numeric
extra_dataset_4$extra_dataset4 <- as.numeric(gsub(",", "", extra_dataset_4$extra_dataset4))


# Non-parametric estimation of pure-premium based on the given estimator
# na.rm is used to remove any missing values or NaNs from the dataset

pure_prem_2 <- 0.65 * mean(extra_dataset_4$extra_dataset4, na.rm = TRUE) + 
  0.20 * sqrt(var(extra_dataset_4$extra_dataset4, na.rm = TRUE)) + 
  0.15 * quantile(extra_dataset_4$extra_dataset4, 0.9, na.rm = TRUE, type = 1)

# Remove the name associated with the quartile value
names(pure_prem_2) <- NULL

pure_prem_2

# Bootstrapping
# Assigning the additional claims data to claim_sizes
claim_sizes = c(extra_dataset_4$extra_dataset4) 

set.seed(123)
# Choosing number of bootstrap samples by trial and error, with the aim to minimize the MSE
num_bootstrap_samples = 9000

# Creating a numeric vector with 'num_bootstrap_samples' elements, all initialized to 0
pure_prem_results = numeric(num_bootstrap_samples)

# For loop to perform bootstrap sampling based on chosen number of bootstrap samples
for (i in 1:num_bootstrap_samples) {
  # Generating a single bootstrap sample by randomly sampling with replacement from the primary data
  bootstrap_sample = sample(claim_sizes, replace = TRUE)
  
  # Calculating the pure premium estimate for the bootstrap samples as per the given estimator
  pure_prem_results[i] = 0.65* mean(bootstrap_sample ) + 0.20 *var(bootstrap_sample )**(0.5)+ 0.15*quantile(bootstrap_sample , 0.9, type=1) 
}


# Initalizing MSE to 0
MSE = 0
 
# Creating a numeric vector with 'num_bootstrap_samples' elements, all initialized to 0
a = numeric(num_bootstrap_samples)

# For loop to calculate MSE for the bootstrap samples
for (j in 1:num_bootstrap_samples) {
# The MSE of the pure premium given as MSE(ˆt) := 1/m * Σ_{k=1 to m} (tk - t)^2 formula is implemented below
    a[j] = ((pure_prem_results[j] - pure_prem_2)^2)/num_bootstrap_samples
    MSE = MSE + a[j]
}
# Using the MSE result for various number of bootstrap samples was narrowed down to 9000


# Installing and loading required packages
if (!require("MASS")) install.packages("MASS")
library(MASS)

# Assigning the primary claims data to claim_sizes
claim_sizes <- extra_dataset_4$extra_dataset4 


# Setting a list of distributions to be fitted from five parametric families via MLE and choosing the ‘best‘ possible fit
distribution_families <- c("normal", "log-normal", "weibull", "logistic", "exponential")

# Initializing variables for storing the best fit
best_fit_name <- NULL
best_fit_aic <- Inf
best_fit_model <- NULL

# Defining a function that fits a distribution to the data and calculates AIC
# Arguments:
#   - dist_name: Name of the distribution to fit (e.g., "normal", "exponential", etc.)
#   - data: The dataset for fitting the distribution
fit_distribution <- function(dist_name, data) {
  if (dist_name == "lognormal") {
    #If the distribution to be fitted is "lognormal" then log-transform data
    data <- log(data[data > 0])
  }
  
  # Attempt to fit the distribution
  fit <- try(fitdistr(data, dist_name), silent = TRUE)
  
  # Check if fit is a 'try-error'
  if(class(fit) == "try-error") {
    return(list(aic = Inf))
  }
  
  # Calculating AIC and return the fit and AIC
  aic <- AIC(fit)
  return(list(fit = fit, aic = aic))
}

# Initializing a data frame to store the fit results for each distribution
fit_details <- data.frame(Distribution = character(), AIC = numeric(), Parameters = character())


# For Loop to fit each distribution and calculate AIC over the five list of distributions
for (dist_name in distribution_families) {
  fit_result <- fit_distribution(dist_name, claim_sizes)
  
  # Check if we successfully got a fit object
  if(!is.na(fit_result$aic)) {
    # Extracting parameter estimates and their values
    param_estimates <- as.character(fit_result$fit$estimate)
    # Combining them into a single string
    params_str <- paste(names(fit_result$fit$estimate), param_estimates, sep = "=", collapse = "; ")
    # Save the fit details in the data frame
    fit_details <- rbind(fit_details, data.frame(Distribution = dist_name, AIC = fit_result$aic, Parameters = params_str))
  }
  
  # Checking if AIC is smaller than the current best AIC
  if (fit_result$aic < best_fit_aic) {
    # Updating the best AIC and the corresponding best-fit distribution
    best_fit_name <- dist_name
    best_fit_aic <- fit_result$aic
    best_fit_model <- fit_result$fit
  }
}

# Print the fit details for each distribution
print(fit_details)

# Output the best fit results
cat("Best-fit distribution:", best_fit_name, "\n")
cat("Best-fit AIC:", best_fit_aic, "\n")
print(best_fit_model)

# Initializing pure_premium variable
pure_premium <- NA

# Calculating the pure premium based on the best-fit distribution
if (!is.null(best_fit_name)) { # Check if best_fit_name is not NULL
  if (best_fit_name == "normal") {  # If the best-fit distribution is normal
    mean_normal <- best_fit_model$estimate["mean"] # Extracting mean parameter
    sd_normal <- best_fit_model$estimate["sd"]# Extracting standard deviation parameter
    # Calculating pure premium using weighted sum of mean, standard deviation, and 90th percentile
    pure_premium <- 0.65 * mean_normal + 0.20 * sd_normal + 0.15 * qnorm(0.9, mean = mean_normal, sd = sd_normal)
    
  } else if (!is.null(best_fit_model) && "meanlog" %in% names(best_fit_model$estimate) && "sdlog" %in% names(best_fit_model$estimate)) {
    # If the best-fit distribution is log-normal and the necessary parameters are present
    
    meanlog <- best_fit_model$estimate["meanlog"] # Extracting meanlog parameter
    sdlog <- best_fit_model$estimate["sdlog"]# Extracting sdlog parameter
    # Ensure the parameters are numeric
    if (is.numeric(meanlog) && is.numeric(sdlog)) {
      
      # Computing mean, standard deviation, and 90th percentile for the log-normal distribution
      
      mean_lognormal <- exp(meanlog + (sdlog^2) / 2)
      sd_lognormal <- sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
      q90_lognormal <- qlnorm(0.9, meanlog = meanlog, sdlog = sdlog)
      
      pure_premium <- 0.65 * mean_lognormal + 0.20 * sd_lognormal + 0.15 * q90_lognormal
    }
    
  } else if (best_fit_name == "logistic") { # If the best-fit distribution is logistic
    location <- best_fit_model$estimate["location"] # Extracting location parameter
    scale <- best_fit_model$estimate["scale"] # Extracting scale parameter
    
    # Calculating mean, standard deviation, and 90th percentile for the logistic distribution
    mean_logistic <- location
    sd_logistic <- scale * pi / sqrt(3)
    q90_logistic <- qlogis(0.9, location = location, scale = scale)
    pure_premium <- 0.65 * mean_logistic + 0.20 * sd_logistic + 0.15 * q90_logistic
    
  } else if (!is.null(best_fit_model) && "shape" %in% names(best_fit_model$estimate) && "scale" %in% names(best_fit_model$estimate)) {
    # If the best-fit distribution is Weibull and the necessary parameters are present
    
    shape <- best_fit_model$estimate["shape"] # Extracting shape parameter
    scale <- best_fit_model$estimate["scale"] # Extracting scale parameter
    
    # Calculating mean, standard deviation, and 90th percentile for the Weibull distribution
    mean_weibull <- scale * gamma(1 + 1 / shape)
    sd_weibull <- sqrt(scale^2 * (gamma(1 + 2 / shape) - (gamma(1 + 1 / shape))^2))
    q90_weibull <- qweibull(0.9, shape = shape, scale = scale)
    pure_premium <- 0.65 * mean_weibull + 0.20 * sd_weibull + 0.15 * q90_weibull
    
  } else if (best_fit_name == "exponential") { # If the best-fit distribution is exponential
    rate <- best_fit_model$estimate["rate"] # Extracting rate parameter
    pure_premium <- 0.65 / rate + 0.20 * sqrt(1 / (rate^2)) + 0.15 * qexp(0.9, rate = rate)
  }
}

# Print the pure premium
cat("Estimated Pure Premium based on the best-fit distribution:", pure_premium, "\n")


#B2----------------------------------------------------------------------------->
# Setting parameters for the insurer's policy underwriting in 2024
num_policies <- 5000
prob_no_claim <- 0.4  # Not used directly, but defined for clarity
prob_one_claim <- 0.6

#A step through of one simulation 
# Simulating claim occurrences for each policy
# rbinom(n, size, prob) where 'n' is the number of observations, 'size' is the number of trials (1 for Bernoulli), and 'prob' is the probability of success
claims_occurrences <- rbinom(num_policies, size = 1, prob = prob_one_claim)

# The 'claims_occurrences' vector now contains a 0 for no claim and 1 for a claim for each policy

# Assigning the additional claim amounts data to combined_claims
combined_claims <- c(extra_dataset_4$extra_dataset4)

# Initializing a vector to store the claim amounts for the policies
claim_amounts <- numeric(num_policies)  # Start with all zeros

# For each policy that has a claim, randomly assign a claim amount from the combined set
has_claim_indices <- which(claims_occurrences == 1)  # Indices of policies with a claim
claim_amounts[has_claim_indices] <- sample(combined_claims, length(has_claim_indices), replace = TRUE)
# 'claim_amounts' now contains the claim amounts for each policy from the claims amounts dataset, with 0 for policies without a claim

# Get statistics for the assigned claim amounts
claim_stats <- summary(claim_amounts[has_claim_indices])
print(claim_stats)
print(length(claim_amounts[has_claim_indices]))

# Aggregate the claim amounts to get the total claims
total_claims <- sum(claim_amounts)

# Print the total claims
print(total_claims)


# Define the number of simulations
num_simulations <- 100000

# Initializing a vector to store the total claims for each simulation
total_claims_simulations <- numeric(num_simulations)

set.seed(123)

# Run the simulations
for (i in 1:num_simulations) {
  # Simulating the claims occurrences
  claims_occurrences_sim <- rbinom(num_policies, size = 1, prob = 0.6)
  
  # Assign claim amounts to the policies with a claim
  claim_amounts_sim <- numeric(num_policies)
  has_claim_indices_sim <- which(claims_occurrences_sim == 1)
  claim_amounts_sim[has_claim_indices_sim] <- sample(combined_claims, length(has_claim_indices_sim), replace = TRUE)
  
  # Sum the claims to get the total claims for the simulation
  total_claims_simulations[i] <- sum(claim_amounts_sim)
}


# Create the histogram of the aggregate loss distribution
# Create the histogram, but suppress the axis labels for now
hist(total_claims_simulations, breaks = 50, main = "Histogram of Aggregate Loss Distribution",
     xlab = "Total Claims", ylab = "Frequency", xaxt = 'n', yaxt = 'n')

# Calculating the range for x and y axes
x_range <- pretty(range(total_claims_simulations), n = 10)
y_range <- pretty((hist(total_claims_simulations, plot = FALSE)$counts), n = 10)

# Format and add x-axis labels with commas
axis(1, at = x_range, labels = format(x_range, big.mark = ",", scientific = FALSE))

# Format and add y-axis labels with commas
axis(2, at = y_range, labels = format(y_range, big.mark = ",", scientific = FALSE))


#B3----------------------------------------------------------------------------->

# Extracting the parameters from the best-fit model calculated as in B1
shape <- best_fit_model$estimate["shape"]
scale <- best_fit_model$estimate["scale"]

# Number of policies with claims based on the binomial distribution
num_policies_with_claims <- rbinom(1, size = 5000, prob = 0.6)

# Simulating claim amounts using the weibull distribution parameters
simulated_claims_weibull <- rweibull(num_policies_with_claims, shape = shape, scale = scale)

# Get statistics for the assigned claim amounts
claim_stats <- summary(simulated_claims_weibull)
print(claim_stats)
print(length(simulated_claims_weibull))
print(sum(simulated_claims_weibull))

# Define the number of simulations for the aggregate loss distribution
num_simulations <- 100000
aggregate_losses <- numeric(num_simulations)
set.seed(123)
for (i in 1:num_simulations) {
  # Simulating the claims for 5,000 policies
  num_claims <- rbinom(5000, size = 1, prob = 0.6)
  # Simulating the claim amounts using the Weibull distribution
  claims <- rweibull(num_claims, shape = shape, scale = scale)
  # Aggregate losses for the current simulation
  aggregate_losses[i] <- sum(claims)
}

# Create a histogram to visualize the simulated aggregate loss distribution
hist(aggregate_losses, breaks = 50, main = "Aggregate Loss Distribution (Weibull)", 
     xlab = "Aggregate Loss", ylab = "Frequency", xaxt='n', yaxt='n')

# Adding formatted axis labels with commas for readability
# Define the breaks for x and y axis
x_breaks <- pretty(aggregate_losses, n = 10)
y_breaks <- pretty(hist(aggregate_losses, plot = FALSE)$counts, n = 10)

# Add formatted x and y axis labels
axis(1, at = x_breaks, labels = format(x_breaks, big.mark = ",", scientific = FALSE))
axis(2, at = y_breaks, labels = format(y_breaks, big.mark = ",", scientific = FALSE))



#C1----------------------------------------------------------------------------->

# Setting seed for reproducibility
set.seed(123)
extra_dataset_4 <- Data[, c(1, 20)]

# Converting string data into numeric
extra_dataset_4$extra_dataset4 <- as.numeric(gsub(",", "", extra_dataset_4$extra_dataset4))

#Method 1 - Using boxplot
# Identifying outliers using a boxplot
boxplot(extra_dataset_4$extra_dataset4, main = "Boxplot of Claim Amounts (Dataset 4)", ylab = "Claim Amount")

#Method 2 - Using Interquartile range
# Computing the first quartile (Q1) of the dataset
Q1 <- quantile(extra_dataset_4$extra_dataset4, 0.25, na.rm = TRUE, type = 1)

# Computing the third quartile (Q3) of the dataset
Q3 <- quantile(extra_dataset_4$extra_dataset4, 0.75, na.rm = TRUE, type = 1)

# Calculating the interquartile range (IQR)
IQR <- Q3 - Q1

# Define the upper bound for outliers
upper_bound <- Q3 + 1.5 * IQR

# Identifying the outliers
outliers <- extra_dataset_4$extra_dataset4[extra_dataset_4$extra_dataset4 > upper_bound]

# Removing the outliers identified above
extra_dataset_4_no_outliers <- extra_dataset_4[!(extra_dataset_4$extra_dataset4 > upper_bound), ]

# Assigning the additional claim amounts data (free of outliers) to combined_claims
combined_claims <- c(extra_dataset_4_no_outliers$extra_dataset4)

# Initializing a vector to store the claim amounts for the policies
claim_amounts <- numeric(num_policies)  # Start with all zeros

# For each policy that has a claim, randomly assign a claim amount from the combined set
has_claim_indices <- which(claims_occurrences == 1)  # Indices of policies with a claim
claim_amounts[has_claim_indices] <- sample(combined_claims, length(has_claim_indices), replace = TRUE)
# 'claim_amounts' now contains the claim amounts for each policy, with 0 for policies without a claim

# Get statistics for the assigned claim amounts
claim_stats <- summary(claim_amounts[has_claim_indices])
print(claim_stats)
print(length(claim_amounts[has_claim_indices]))

# Aggregate the claim amounts to get the total claims
total_claims <- sum(claim_amounts)

# Print the total claims
print(total_claims)

# Define the number of simulations
num_simulations <- 100000

# Initializing a vector to store the total claims for each simulation
total_claims_simulations <- numeric(num_simulations)

set.seed(123)
# Run the simulations
for (i in 1:num_simulations) {
  # Simulating the claims occurrences
  claims_occurrences_sim <- rbinom(num_policies, size = 1, prob = 0.6)

  # Assign claim amounts to the policies with a claim
  claim_amounts_sim <- numeric(num_policies)
  has_claim_indices_sim <- which(claims_occurrences_sim == 1)
  claim_amounts_sim[has_claim_indices_sim] <- sample(combined_claims, length(has_claim_indices_sim), replace = TRUE)
  
  # Add the claims to get the total claims for the simulation
  total_claims_simulations[i] <- sum(claim_amounts_sim)
}


# Create the histogram of the aggregate loss distribution
# Create the histogram, but suppress the axis labels for now
hist(total_claims_simulations, breaks = 50, main = "Histogram of Aggregate Loss Distribution",
     xlab = "Total Claims", ylab = "Frequency", xaxt = 'n', yaxt = 'n')

# Calculating the range for x and y axes
x_range <- pretty(range(total_claims_simulations), n = 10)
y_range <- pretty(0:max(hist(total_claims_simulations, plot = FALSE)$counts), n = 10)

# Format and add x-axis labels with commas
axis(1, at = x_range, labels = format(x_range, big.mark = ",", scientific = FALSE))

# Format and add y-axis labels with commas
axis(2, at = y_range, labels = format(y_range, big.mark = ",", scientific = FALSE))


#C2----------------------------------------------------------------------------->

#Method 1 - Using boxplot
# Identifying outliers using a boxplot
boxplot(extra_dataset_4$extra_dataset4, main = "Boxplot of additional claim amounts (Extradataset 4)", ylab = "Claim Amount")

#Method 2 - Using Interquartile range
# Computing the first quartile (Q1) of the dataset
Q1 <- quantile(extra_dataset_4$extra_dataset4, 0.25, na.rm = TRUE, type = 1)

# Computing the third quartile (Q3) of the dataset
Q3 <- quantile(extra_dataset_4$extra_dataset4, 0.75, na.rm = TRUE, type = 1)

# Calculating the interquartile range (IQR)
IQR <- Q3 - Q1

# Define the upper bound for outliers
upper_bound <- Q3 + 1.5 * IQR

# Identifying the outliers
outliers <- extra_dataset_4$extra_dataset4[extra_dataset_4$extra_dataset4 > upper_bound]

# Removing the outliers identified above
extra_dataset_4_no_outliers <- extra_dataset_4[!(extra_dataset_4$extra_dataset4 > upper_bound), ]

# Assigning the additional claim amounts data (free of outliers) to combined_claims
combined_claims <- c(extra_dataset_4_no_outliers$extra_dataset4)

# Initializing a vector to store the claim amounts for the policies
claim_amounts <- numeric(num_policies)  # Start with all zeros
# Extracting the parameters for the weibull distribution for claim amounts after removing outliers
model <- fit_distribution("weibull", combined_claims)
param_estimates <- as.character(model$fit$estimate)

#Storing the parameters
shape <- as.numeric(param_estimates[1]) 
scale <- as.numeric(param_estimates[2])

# Number of policies with claims based on the binomial distribution
num_policies_with_claims <- rbinom(1, size = 5000, prob = 0.6)

set.seed(123)
# Simulating claim amounts using the weibull distribution parameters
simulated_claims_weibull <- rweibull(num_policies_with_claims, shape = shape, scale = scale)

# Get statistics for the assigned claim amounts
claim_stats <- summary(simulated_claims_weibull)
print(claim_stats)
print(length(simulated_claims_weibull))
print(sum(simulated_claims_weibull))

# Define the number of simulations for the aggregate loss distribution
num_simulations <- 100000
aggregate_losses <- numeric(num_simulations)

set.seed(123)
for (i in 1:num_simulations) {
  # Simulating the claims for 5,000 policies
  num_claims <- rbinom(5000, size = 1, prob = 0.6)
  # Simulating the claim amounts using the Weibull distribution
  claims <- rweibull(num_claims, shape = shape, scale = scale)
  # Aggregate losses for the current simulation
  aggregate_losses[i] <- sum(claims)
}

# Create a histogram to visualize the simulated aggregate loss distribution
hist(aggregate_losses, breaks = 50, main = "Aggregate Loss Distribution (Weibull)", 
     xlab = "Aggregate Loss", ylab = "Frequency", xaxt='n', yaxt='n')

# Adding formatted axis labels with commas for readability
# Define the breaks for x and y axis
x_breaks <- pretty(aggregate_losses, n = 10)
y_breaks <- pretty(hist(aggregate_losses, plot = FALSE)$counts, n = 10)

# Add formatted x and y axis labels
axis(1, at = x_breaks, labels = format(x_breaks, big.mark = ",", scientific = FALSE))
axis(2, at = y_breaks, labels = format(y_breaks, big.mark = ",", scientific = FALSE))

