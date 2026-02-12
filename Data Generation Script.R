# Data generation MC Poster
set.seed(123)

n  <- 200      # number of observations
p_true <- 2    # true predictors
p_noise <- 3   # irrelevant predictors
sigma <- 1     # noise SD

# geneating X Matrix
X_true <- matrix(rnorm(n * p_true), n, p_true)
colnames(X_true) <- paste0("X", 1:p_true)
# selecting betas
beta <- sample(as.numeric(-2:2), size = p_true)  # true coefficients
# generating error terms
epsilon <- rnorm(n, mean = 0, sd = sigma)

# calculating Y
y <- X_true %*% beta + epsilon

# generating noise variables
X_noise <- matrix(rnorm(n * p_noise), n, p_noise)
colnames(X_noise) <- paste0("XF", 1:p_noise)


X_full <- cbind(X_true, X_noise)

data_sim <- data.frame(y, X_full)

# true model fit
true_vars <- paste0("X", 1:p_true)
fit_true <- lm(
  as.formula(paste("y ~", paste(true_vars, collapse = " + "))),
  data = data_sim
)
summary(fit_true)

# full model fit
fit_full <- lm(y ~ ., data = data_sim)
summary(fit_full)



