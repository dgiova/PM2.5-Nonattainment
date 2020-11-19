compute_var <- function(Z, N, y) {
  y_bar <- sum(Z * y) / N
  var_Z <- 1/(N-1) * sum((Z * (y-y_bar)^2))
  return(var_Z)
}

dat <- subset(fulldat, outofrange==0)
df <- data.frame(na.omit(dat[,c("a","pmfu","pscat")]))

num_cat <- 4
strata_dim <- rep(0,num_cat)
strata_var <- rep(0,num_cat)
strata_weight <- rep(0,num_cat)

idx <- 1
for (cat in levels(df$pscat)) {
  z_strata <- df[df$pscat == cat, c("a", "pmfu")]
  N_1 <- sum(z_strata[,"a"])
  N_0 <- nrow(z_strata) - N_1
  Z <- z_strata[,"a"]
  y <- z_strata[,"pmfu"]
  
  strata_dim[idx] <- (1/N_1 * sum(Z * y) -
                        1/N_0 * sum((1 - Z) * y))
  
  V_1 <- compute_var(Z, N_1, y)
  V_0 <- compute_var(1-Z, N_0, y)
  strata_var[idx] <- V_1/N_1 + V_0/N_0
  
  strata_weight[idx] <- nrow(z_strata) / nrow(df)
  
  idx = idx + 1
}

tau_dim <- sum(strata_dim * strata_weight)
var_dim <- sum(strata_var * strata_weight^2)

cat("DIM per stratum:", strata_dim, "\n",
    "Variance per stratum:", strata_var, "\n",
    "Stratified DIM estimate of ATE =", tau_dim, "\n",
    "Estimate of the variance =", var_dim)
