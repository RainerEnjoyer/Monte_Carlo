# Datensatz generieren
set.seed(1)

# 1) n
n <- 100

# 2) beta (Parametervektor)
beta <- c(2, 5, 0, 0, 0)
p <- length(beta)  # = 5

# 3) X ~ N_p(0, I)
# Bei Sigma = I sind die Kovariaten unkorreliert (und bei Normalverteilung: unabhängig)
X <- matrix(rnorm(n * p), nrow = n, ncol = p)
colnames(X) <- paste0("X", 1:p)

# 4) Residuen epsilon ~ N(0,1)
eps <- rnorm(n, mean = 0, sd = 1)

# 5) y = X beta + eps
y <- as.vector(X %*% beta + eps)

# 6) y und X mergen
dat <- data.frame(y = y, X)

# Check
head(dat)

# Modellfit - Beispiel
fit_full <- lm(y ~ X1 + X2 + X3 + X4 + X5, data = dat)
summary(fit_full)
coef(fit_full)

# Prediciton Error
# Anzahl Testbeobachtungen
Ntest <- 100

# 1) X_new ~ N_p(0, I)
X_new <- matrix(rnorm(Ntest * p), nrow = Ntest, ncol = p)
colnames(X_new) <- paste0("X", 1:p)

# 2) y_new = X_new %*% beta + eps_new, eps_new ~ N(0,1)
eps_new <- rnorm(Ntest, mean = 0, sd = 1)
y_new  <- as.vector(X_new %*% beta + eps_new)

testdat <- data.frame(y = y_new, X_new)

# 3) Vorhersage mit einem trainierten Modell (z.B. fit_a oder fit_full)
y_hat <- predict(fit_full, newdata = testdat)

# 4) Prediction Error (MSE auf Testdaten)
PE_hat <- mean((testdat$y - y_hat)^2)
PE_hat

library(MASS)

set.seed(1)

# -------------------------
# DGP
# -------------------------
generate_data <- function(n, p, beta, sigma, Sigma = diag(p)) {
  X <- mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)
  eps <- rnorm(n, 0, sigma)
  y <- as.vector(X %*% beta + eps)
  data.frame(y = y, X)
}

# -------------------------
# Step 1: Screening (Top-K) via grobem Bootstrap-PE
# -------------------------
step1_screen <- function(dat, p, B1, K) {
  n <- nrow(dat)
  
  # alle nicht-leeren Modelle
  A <- unlist(lapply(1:p, function(k) combn(1:p, k, simplify = FALSE)), recursive = FALSE)
  
  boot_pe_one_model <- function(a) {
    rhs  <- paste0("X", a, collapse = " + ")
    form <- as.formula(paste("y ~", rhs))
    
    pe_b <- numeric(B1)
    for (b in 1:B1) {
      idx <- sample.int(n, size = n, replace = TRUE)
      dat_b <- dat[idx, , drop = FALSE]
      fit_b <- lm(form, data = dat_b)
      yhat  <- predict(fit_b, newdata = dat)
      pe_b[b] <- mean((dat$y - yhat)^2)
    }
    mean(pe_b)
  }
  
  PE1 <- sapply(A, boot_pe_one_model)
  A[order(PE1)[1:K]]
}

# -------------------------
# Step 2: m-out-of-n Bootstrap -> Selektionswahrscheinlichkeiten
# Gewinner pro Bootstrap-Run = Modell mit minimalem PE (MSE) auf Originaldaten
# -------------------------
selection_prob_m <- function(A_screen, dat, m, B2) {
  n <- nrow(dat)
  
  # model labels
  models_chr <- sapply(A_screen, function(a) paste0("{", paste(a, collapse = ","), "}"))
  wins <- setNames(integer(length(A_screen)), models_chr)
  
  # Formeln einmal bauen
  forms <- lapply(A_screen, function(a) {
    rhs <- paste0("X", a, collapse = " + ")
    as.formula(paste("y ~", rhs))
  })
  
  for (b in 1:B2) {
    idx <- sample.int(n, size = m, replace = TRUE)
    dat_b <- dat[idx, , drop = FALSE]
    
    pe_vec <- numeric(length(A_screen))
    for (j in seq_along(A_screen)) {
      fit_b <- lm(forms[[j]], data = dat_b)
      yhat  <- predict(fit_b, newdata = dat)
      pe_vec[j] <- mean((dat$y - yhat)^2)
    }
    
    j_star <- which.min(pe_vec)
    wins[j_star] <- wins[j_star] + 1L
  }
  
  prob <- wins / B2
  out <- data.frame(model = names(prob), sel_prob = as.numeric(prob), m = m, row.names = NULL)
  out[order(-out$sel_prob), ]
}

# -------------------------
# Monte Carlo Einbettung
# -------------------------
R <- 100
n <- 100
beta <- c(2, 5, 0, 0, 0)
p <- length(beta)
Sigma <- diag(p)
sigma <- 1

B1 <- 50
K  <- 6
B2 <- 300

m_grid <- c(50, 70, 85, 100)  # m4 = n

# Zählstruktur: pro m zählen wir, welches Modell in der Replikation "gewinnt"
win_counts <- setNames(vector("list", length(m_grid)), paste0("m=", m_grid))
for (nm in names(win_counts)) win_counts[[nm]] <- integer(0)

for (r in 1:R) {
  dat_r <- generate_data(n, p, beta, sigma, Sigma)
  
  # Step 1: Screening
  A_screen_r <- step1_screen(dat_r, p, B1, K)
  
  # Step 2: pro m Selektionswahrscheinlichkeiten; Replikations-Gewinner = max sel_prob
  for (m_val in m_grid) {
    sel_r <- selection_prob_m(A_screen_r, dat_r, m = m_val, B2 = B2)
    winner <- sel_r$model[1]  # weil sel_r bereits absteigend sortiert ist
    
    key <- paste0("m=", m_val)
    vec <- win_counts[[key]]
    vec[winner] <- if (winner %in% names(vec)) vec[winner] + 1L else 1L
    win_counts[[key]] <- vec
  }
}

# Ergebnis-Tabelle: Gewinnerwahrscheinlichkeit über Replikationen
mc_res <- do.call(rbind, lapply(m_grid, function(m_val) {
  key <- paste0("m=", m_val)
  counts <- win_counts[[key]]
  data.frame(
    m = m_val,
    model = names(counts),
    wins = as.integer(counts),
    win_prob_over_R = as.numeric(counts) / R,
    row.names = NULL
  )
}))

mc_res <- mc_res[order(mc_res$m, -mc_res$win_prob_over_R), ]
mc_res

mc_res$is_true <- mc_res$model == "{1,2}"
mc_res$model_lab <- ifelse(mc_res$is_true, "a0 = {1,2}", mc_res$model)

# m als Faktor (für saubere Achsen-Reihenfolge)
mc_res$m_f <- factor(mc_res$m, levels = sort(unique(mc_res$m)))
library(ggplot2)

ggplot(mc_res,
       aes(x = as.numeric(as.character(m_f)),
           y = win_prob_over_R,
           group = model,
           color = is_true,
           size = is_true,
           alpha = is_true)) +
  geom_line() +
  geom_point() +
  scale_color_manual(values = c("grey60", "red")) +
  scale_size_manual(values = c(0.7, 1.4)) +
  scale_alpha_manual(values = c(0.6, 1)) +
  labs(
    x = "Bootstrap sample size m (n = 100)",
    y = "Selection probability",
    color = "True model",
    title = "Model selection probabilities over m",
    subtitle = "True model: a0 = {1,2}"
  ) +
  theme_minimal(base_size = 13) +
  theme(legend.position = "none")


ggplot(mc_res,
       aes(x = m_f,
           y = win_prob_over_R,
           fill = is_true)) +
  geom_col(width = 0.75,
           color = "white") +
  geom_text(
    aes(label = paste0(round(100 * win_prob_over_R), "%")),
    position = position_stack(vjust = 0.5),
    size = 3,
    color = "black"
  ) +
  scale_fill_manual(
    values = c("grey75", "#C0392B"),
    guide = "none"
  ) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(
    x = "Bootstrap sample size m (n = 100)",
    y = "Selection probability",
    title = "Model selection probabilities",
    subtitle = "True model: a0 = {1,2}"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    panel.grid.major.x = element_blank(),
    panel.grid.minor = element_blank()
  )

ggplot(mc_res, aes(x = m_f, y = reorder(model, win_prob_over_R), fill = win_prob_over_R)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_tile(
    data = subset(mc_res, is_true),
    fill = NA, color = "black", linewidth = 1
  ) +
  labs(x = "m (n = 100)", y = "Model", fill = "Sel. prob.") +
  theme_minimal(base_size = 12)


