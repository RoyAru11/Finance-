install.packages("RSD")
library(RSD)
library(ggplot2)
# Parameters
set.seed(42)  # For reproducibility; change for new data
num_env <- 27
num_crops <- 5  # Or 4 for fewer crops
crops <- paste0("Crop", 1:num_crops)
environments <- paste0("Env", 1:num_env)
mean_yield <- 2500  # Average yield
sd_yield <- 600     # Variability

# Generate random yields
yields <- matrix(rnorm(num_env * num_crops, mean = mean_yield, sd = sd_yield),
                 nrow = num_env, ncol = num_crops)

# Ensure positive yields (>= 100)
yields[yields < 100] <- 100

# Create DataFrame
df <- as.data.frame(yields)
colnames(df) <- crops
rownames(df) <- environments

# Print rounded values for readability
print(round(df, 2))

# Optional: Save to CSV for use in RSD package
# write.csv(df, "crop_yields.csv", row.names = TRUE)
pr <- rep(1/num_env,num_env)
outcome1 <- df$Crop1
outcome2 <- df$Crop4
pair <- createStochasticDominance(outcome1, outcome2, pr, pr)

fsd.test(pair)
ssd.test(pair)

fsd.plot(pair, names=c("Crop1","Crop2"))+
  labs(
    x= "Yield",
    y="CDF",
    color="Crops")

ssd.plot(pair, names=c("Crop1","Crop2"))+
  labs(
    x= "Yield",
    y="Integral of CDF",
    color="Crops")
assd.test(pair, type = 'll')
assd.test(sd.obj, type = 'ths')
