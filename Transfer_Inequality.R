# Transfer Sensitivity, Lorenz curves and TSD check
set.seed(42)
install.packages("ggplot2")
install.packages("RSD")
library(RSD)
library(gridExtra)
library(ggplot2)
library(ineq)

num_people <- 20
num_groups <- 5
people <- paste0("Person", 1:num_people)
groups <- paste0("Group", 1:num_groups)

mean_income <- 5000
sd_income <- 600

incomes <- matrix(rnorm(num_people*num_groups, mean = mean_income,
                        sd = sd_income), nrow= num_people, 
                        ncol = num_groups)
incomes <- pmax(incomes, 1000)
df_income <- data.frame(incomes)
rownames(df_income) <- people
colnames(df_income) <- groups
print(df_income)

# Run the experiment on Group 1
x_orig <- df_income[,"Group1"]

# Coefficient of variation (CV)
cv <- function(x) sd(x)/mean(x)

# Gini 
gini <- function(x) ineq::Gini(x)

# Atkinson with epsilon parameter
atkinson <- function(x, epsilon = 0.5) ineq::Atkinson(x, parameter = epsilon)

# Lorenz curve data for plotting  
lorenz_df <- function(x){
  x_sorted <- sort(x)
  n <- length(x_sorted)
  cum_income <- cumsum(x_sorted)
  total <- sum(x_sorted)
  p <- (1:n)/n
  L <- cum_income / total #cumulative share
  # include the origin (0,0)
  df <- data.frame(p=c(0,p), L=c(0,L))
  return(df)
}

#Progressive transfer
progressive_transfer <- function(vec, from_idx, to_idx, amount){
  vec2 <- vec
  vec2[from_idx] <- vec2[from_idx] - amount
  vec2[to_idx] <- vec2[to_idx] + amount
  return(vec2)
}
richest_idx <- function(vec) which.max(vec)
poorest_idx <- function(vec) which.min(vec)

# S1(z) = sum_i (z - x_i)
S1_of_z <- function(x,z) sum(pmax(z-x, 0))
# S2(z) = sum_i (z - x_i)_+^2
S2_of_z <- function(x,z) sum((pmax(z-x, 0))^2)

compute_S1_S2_grid <- function(x, zgrid){
  s1 <- sapply(zgrid, function(z) S1_of_z(x,z))
  s2 <- sapply(zgrid, function(z) S2_of_z(x,z))
  data.frame(z=zgrid, S1=s1, S2=s2)
}
# Checking SSD between x and y
# Here we adopt the discrete S1 check: 
#x SSD y iff S1_y(z) >= S1_x(z) for all z and > for some z.

check_SSD <- function(before, after, zgrid){
  s_before <- compute_S1_S2_grid(before, zgrid)$S1
  s_after  <- compute_S1_S2_grid(after, zgrid)$S1
  le_all <- all((s_after - s_before) <= 1e-8)   # after ≤ before
  strict_some <- any((s_after - s_before) < -1e-8)
  list(le_all = le_all, strict_some = strict_some, diff = s_after - s_before)
}

check_TSD <- function(before, after, zgrid){
  s_before <- compute_S1_S2_grid(before, zgrid)$S2
  s_after  <- compute_S1_S2_grid(after, zgrid)$S2
  le_all <- all((s_after - s_before) <= 1e-8)   # after ≤ before
  strict_some <- any((s_after - s_before) < -1e-8)
  list(le_all = le_all, strict_some = strict_some, diff = s_after - s_before)
}

# --- Baseline inequality metrics --------------------------------------------
cat("Baseline inequality metrics for Group1:\n")
cat("  Mean:", round(mean(x_orig),2), " SD:", round(sd(x_orig),2), " CV:", round(cv(x_orig),4), "\n")
cat("  Gini:", round(gini(x_orig),4), "  Atkinson(eps=0.5):", round(atkinson(x_orig, 0.5),4), "\n\n")

# --- Lorenz curves before transfer ------------------------------------------
lor_x <- lorenz_df(x_orig)
# Apply a single progressive transfer --------------------------
donor <- richest_idx(x_orig)
recipient <- poorest_idx(x_orig)
transfer_amt <- 500

x_after <- progressive_transfer(x_orig, donor, recipient, transfer_amt)
cat("Single transfer: from", people[donor], "to", people[recipient], "amount =", transfer_amt, "\n")
cat("  New mean:", round(mean(x_after),2), " CV:", round(cv(x_after),4), " Gini:", round(gini(x_after),4), "\n\n")

# Sequence of multiple progressive transfers 
x_multi <- x_orig
steps <- 10
amt_each <- 100

cat("Simulating", steps, "iterative progressive transfers of",
    amt_each,"each:\n")
for(i in 1:steps){
  r <- richest_idx(x_multi)
  p <- poorest_idx(x_multi)
  amt <- min(amt_each, x_multi[r] - 1e-6)
  x_multi <- progressive_transfer(x_multi, r, p, amt)
  cat(sprintf("  Step %2d: donor=%s recipient=%s CV=%.4f Gini=%.4f\n",
              i, people[r], people[p], cv(x_multi), gini(x_multi)))
}
cat("\n")

# Lorenz curves after transfers ---------------------------------
lor_after1  <- lorenz_df(x_after)
lor_aftermul <- lorenz_df(x_multi)

p_lor <- ggplot() +
  geom_abline(slope=1, intercept=0, linetype="dotted", color="black") +
  geom_line(data=lor_x, aes(x=p, y=L, color="Original"), linewidth=1) +
  geom_line(data=lor_after1, aes(x=p, y=L, color="After 1 Transfer"), linewidth=1) +
  geom_line(data=lor_aftermul, aes(x=p, y=L, color="After Multiple Transfers"), linewidth=1) +
  labs(title="Lorenz Curves (closer to 45° line = more equality)",
       x="Cumulative population share", y="Cumulative income share", color="Legend") +
  theme_minimal()

# Compute S1 and S2 grids and check SSD/TSD
zmin <- min(c(x_orig, x_after, x_multi)) - 50
zmax <- max(c(x_orig, x_after, x_multi)) + 50
zgrid <- seq(zmin, zmax, length.out = 400)

ssd_check_1 <- check_SSD(x_orig, x_after, zgrid)  # is after1 SSD-dominant relative to orig?
tsd_check_1 <- check_TSD(x_orig, x_after, zgrid)  # is after1 TSD-dominant relative to orig?

cat("SSD check (after multiple vs original): S1_after >= S1_orig for all z? ", ssd_check_1$ge_all,
    " and strict for some z? ", ssd_check_1$strict_some, "\n")
cat("TSD check (after multiple vs original): S2_after >= S2_orig for all z? ", tsd_check_1$ge_all,
    " and strict for some z? ", tsd_check_1$strict_some, "\n\n")

# --- Plot S1 and S2 functions for visual inspection -------------------------
S_orig <- compute_S1_S2_grid(x_orig, zgrid)
S_after1 <- compute_S1_S2_grid(x_after, zgrid)
S_multi <- compute_S1_S2_grid(x_multi, zgrid)

df_S1 <- data.frame(z = rep(zgrid, 3),
                    S1 = c(S_orig$S1, S_after1$S1, S_multi$S1),
                    series = rep(c("orig", "after", "after_multi"),
                               each = length(zgrid)))
df_S2 <- data.frame(z = rep(zgrid, 3),
                   S2 = c(S_orig$S2, S_after1$S2, S_multi$S2),
                   series = rep(c("orig", "after", "after_multi"),
                                each = length(zgrid)))
p_S1 <- ggplot(df_S1, aes(x=z, y=S1, color=series)) + 
  geom_line(linewidth = 1) + theme_minimal() +
  labs(title = "S1(z) = Σ (z - x_i)_+  (used for SSD checks)",
       x = "z (income threshold)", y = "S1(z)")

p_S2 <- ggplot(df_S2, aes(x = z, y = S2, color = series)) +
  geom_line(linewidth = 1) + theme_minimal() +
  labs(title = "S2(z) = Σ (z - x_i)_+^2  (used for TSD checks)",
       x = "z (income threshold)", y = "S2(z)")

# Combine plots
gridExtra::grid.arrange(p_lor, p_S1, p_S2, ncol = 2)

# --- Print a short summary --------------------------------------------------
cat("Summary of results (Group1):\n")
cat(" - Original CV:", round(cv(x_orig),4), " Gini:", round(gini(x_orig),4), "\n")
cat(" - After 1 transfer CV:", round(cv(x_after),4), " Gini:", round(gini(x_after),4), "\n")
cat(" - After multiple transfers CV:", round(cv(x_multi),4), " Gini:", round(gini(x_multi),4), "\n\n")

cat("Interpretation notes:\n")
cat("  * A decrease in CV/Gini/Atkinson after progressive transfers indicates inequality falling.\n")
cat("  * SSD check uses S1: if S1_after(z) >= S1_orig(z) for all z (and > for some z), then 'after' SSD-dominates 'orig'.\n")
cat("  * TSD check uses S2: if S2_after(z) >= S2_orig(z) for all z (and > for some z), then 'after' TSD-dominates 'orig'.\n")


#Extra
pr <- rep(1/20,20)
outcome1 <- x_orig
outcome2 <- x_multi
pair <- createStochasticDominance(outcome1, outcome2, pr, pr)
ssd.test(pair)
ssd.plot(pair, names=c("Original","Multiple"))+
  labs(
    x= "Incomes",
    y="Integral of CDF",
    color="Group1")

