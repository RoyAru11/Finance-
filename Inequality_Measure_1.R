install.packages("RSD")
library(RSD)
library(ggplot2)
# Load necessary libraries (base R is sufficient)
set.seed(42)  # For reproducibility
num_env <- 27
num_crops <- 5
crops <- paste0('Crop', 1:num_crops)
environments <- paste0('Env', 1:num_env)
mean_yield <- 2500
sd_yield <- 600

# Generate random yields
yields <- matrix(rnorm(num_env * num_crops, mean = mean_yield, sd = sd_yield), 
                 nrow = num_env, ncol = num_crops)
yields <- pmax(yields, 100)  # Ensure positive yields
print(yields)
# Create data frame
df <- data.frame(yields, row.names = environments)
colnames(df) <- crops

# Print rounded data frame
print(round(df, 2))

# Calculate mean and standard deviation for each crop
crop_means <- colMeans(df)
crop_sds <- apply(df, 2, sd)
crop_cv <- crop_sds / crop_means

# Print CV
print(round(crop_cv, 4)) # Crop4 least variability and Crop1 most

# Recalculate CV for Crop1
crop1_sd <- sd(c(df[, "Crop1"]))
crop1_mean <- mean(c(df[, "Crop1"]))
crop1_cv <- crop1_sd / crop1_mean
print(round(crop1_cv, 4))
print(df[,"Crop1"])
# Simulate progressive transfer for Crop1
df_new <- df
df_new["Env1", "Crop1"] <- df["Env1", "Crop1"] - 100
df_new["Env2", "Crop1"] <- df["Env2", "Crop1"] + 100
print(df_new[,"Crop1"])

crop1_sd_new <- sd(c(df_new[, "Crop1"]))
crop1_mean_new <- mean(c(df_new[,"Crop1"]))
crop1_cv_new <- crop1_sd_new / crop1_mean_new
# New CV
print(round(crop1_cv_new,4))











