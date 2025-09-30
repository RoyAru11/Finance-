install.packages("ggplot2")
install.packages("RSD")
library(RSD)
library(gridExtra)
library(ggplot2)
library(ineq)
set.seed(123)

# --- Functions ---------------------------------------------------------

# phi(v; z) = sum (z - v_i)_+ (first-order shortfall)
phi <- function(v, z) {
  sum(pmax(z - v, 0))
}

# Phi(v; z) = sum (z - v_i)_+^2 (second-order shortfall)
Phi <- function(v, z) {
  sum(pmax(z - v, 0)^2)
}

# Difference functions for dominance
g <- function(x, y, z) { phi(y, z) - phi(x, z) }  # for SSD
G <- function(x, y, z) { Phi(y, z) - Phi(x, z) }  # for TSD

pr <- rep(1/20, 20)
x_orig <- round(runif(20, 1000, 10000), 0)  # 20 random incomes
cat("Original incomes:\n"); print(x_orig)

# Progressive transfer: richest -> poorest (with condition 4b: gap > amt)
progressive_transfer <- function(vec, amt) {
  v <- vec
  richest <- which.max(v)
  poorest <- which.min(v)
  gap <- v[richest] - v[poorest]
  amt <- min(amt, gap - 1e-8)  # Ensure gap > amt (condition 4b)
  v[richest] <- v[richest] - amt
  v[poorest] <- v[poorest] + amt
  return(v)
}

# Regressive transfer: poorest -> richest (with condition 4b)
regressive_transfer <- function(vec, amt) {
  v <- vec
  richest <- which.max(v)
  poorest <- which.min(v)
  gap <- v[richest] - v[poorest]
  amt <- min(amt, gap - 1e-8)  # Ensure gap > amt (condition 4b)
  v[richest] <- v[richest] + amt
  v[poorest] <- v[poorest] - amt
  return(v)
}

# Favorable Composite Transfer (FACT): Progressive + smaller regressive at higher level
composite_transfer <- function(vec, amt_prog, amt_reg) {
  v <- progressive_transfer(vec, amt_prog)
  v_sorted <- sort(v, index.return = TRUE)
  second_richest_idx <- v_sorted$ix[length(v_sorted$ix) - 1]
  second_poorest_idx <- v_sorted$ix[2]
  gap_reg <- v[second_richest_idx] - v[second_poorest_idx]
  amt_reg <- min(amt_reg, gap_reg - 1e-8)
  v[second_richest_idx] <- v[second_richest_idx] - amt_reg
  v[second_poorest_idx] <- v[second_poorest_idx] + amt_reg
  return(v)
}

# Check TSD dominance with improved tolerance
check_dom <- function(Gvals, tol = 1e-10) {
  all(Gvals >= -tol) & any(Gvals > tol)
}

# Sequential transfers
steps <- 10
x_seq_fact <- list(sort(x_orig))
print(x_seq_fact)

for (i in 1:steps) {
  x_seq_fact[[i + 1]] <- composite_transfer(x_seq_fact[[i]], amt_prog = 100, amt_reg = 100)
}
print(x_seq_fact)

# --- Diagnostics after multiple transfers ------------------------------
zgrid <- seq(min(unlist(x_seq_fact)) - 100, max(unlist(x_seq_fact)) + 100, length.out = 400)

# Compute G(z) differences after last step
Phi_x <- sapply(zgrid, function(z) Phi(x_orig, z))
Phi_fact <- sapply(zgrid, function(z) Phi(x_seq_fact[[steps + 1]], z))

G_fact <- Phi_fact - Phi_x

cat("\nTSD dominance after", steps, "FACT steps:\n")
cat("FACT transfers TSD? ", check_dom(G_fact), "\n")

df_G <- data.frame(z = zgrid, G = G_fact)
ggplot(df_G, aes(x = z, y = G)) +
  geom_line(color = "darkred", linewidth = 1) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  labs(title = "TSD Check: G(z) for Original vs FACT (Step 10)",
       x = "z (income threshold)", y = "G(z)") +
  theme_minimal()

# Checking SSD for composite favourable transfers (FACTs)
outcome1 <- x_orig
outcome2 <- x_seq_fact[[11]]
pair <- createStochasticDominance(outcome1, outcome2, pr, pr)
ssd.test(pair)
ssd.plot(pair, names = c("Original", "FACTs")) +
  labs(x = "Incomes", y = "Integral of CDF", color = "Distribution")

# Coefficient of variation (CV) = sd/mean
cv <- function(x) sd(x) / mean(x)

cv_orig <- cv(x_orig)
cv_fact <- cv(x_seq_fact[[11]])

cat("\nCoefficient of Variation:\n")
cat("Original CV =", round(cv_orig, 4), "\n")
cat("FACT CV     =", round(cv_fact, 4), "\n")
 
if (cv_fact < cv_orig) {
  cat("=> Inequality decreased after FACT transfers (CV fell)\n")
} else if (cv_fact > cv_orig) {
  cat("=> Inequality increased after FACT transfers (CV rose)\n")
} else {
  cat("=> No change in inequality (CV same)\n")
}


