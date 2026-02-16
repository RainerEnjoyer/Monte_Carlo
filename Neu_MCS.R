library(MASS)
library(ggplot2)
library(dplyr)
library(forcats)

set.seed(1)

# -------------------------
# DGP
# -------------------------
n <- 100
beta <- c(2, 5, 0, 0, 0)
p <- length(beta)  
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
      idx_oob <- setdiff(1:n, unique(idx))
      if (length(idx_oob) < 5) { pe_b[b] <- NA_real_; next }
      
      dat_b   <- dat[idx, , drop = FALSE]
      dat_oob <- dat[idx_oob, , drop = FALSE]
      
      fit_b <- lm(form, data = dat_b)
      yhat  <- predict(fit_b, newdata = dat_oob)
      pe_b[b] <- mean((dat_oob$y - yhat)^2)
    }
    mean(pe_b, na.rm = TRUE)
  }
  
  PE1 <- sapply(A, boot_pe_one_model)
  A[order(PE1)[1:K]]
}

# -------------------------
# Step 2: m-out-of-n Bootstrap -> Selektionswahrscheinlichkeiten
# Gewinner pro Bootstrap-Run = Modell mit minimalem PE (MSE) auf Originaldaten
# -------------------------

selection_prob_m_PE <- function(A_screen, dat, m, B2) {
  n <- nrow(dat)
  
  models_chr <- sapply(A_screen, function(a) paste0("{", paste(a, collapse = ","), "}"))
  wins_PE <- setNames(integer(length(A_screen)), models_chr)
  
  forms <- lapply(A_screen, function(a) {
    rhs <- paste0("X", a, collapse = " + ")
    as.formula(paste("y ~", rhs))
  })
  
  B2_eff <- 0L
  
  for (b in 1:B2) {
    idx <- sample.int(n, size = m, replace = TRUE)
    idx_oob <- setdiff(1:n, unique(idx))
    if (length(idx_oob) < 5) next
    
    dat_b   <- dat[idx, , drop = FALSE]
    dat_oob <- dat[idx_oob, , drop = FALSE]
    
    pe_vec <- numeric(length(A_screen))
    
    for (j in seq_along(A_screen)) {
      fit_b <- lm(forms[[j]], data = dat_b)
      
      # OOB-MSE (kein Overlap)
      yhat_oob <- predict(fit_b, newdata = dat_oob)
      mse_oob  <- mean((dat_oob$y - yhat_oob)^2)
      
      # Wenn du "optimism-style" behalten willst:
      # yhat_boot <- predict(fit_b, newdata = dat_b)
      # mse_boot  <- mean((dat_b$y - yhat_boot)^2)
      # pe_vec[j] <- 2*mse_oob - mse_boot
      
      # Empfehlung: direkt OOB-MSE als PE
      pe_vec[j] <- mse_oob
    }
    
    j_star <- which.min(pe_vec)
    wins_PE[j_star] <- wins_PE[j_star] + 1L
    B2_eff <- B2_eff + 1L
  }
  
  if (B2_eff == 0L) stop("Kein gültiger OOB-Rep (idx_oob immer zu klein).")
  
  prob <- wins_PE / B2_eff
  data.frame(model = names(prob), sel_prob = as.numeric(prob), m = m, row.names = NULL) |>
    dplyr::arrange(dplyr::desc(sel_prob))
}

# -------------------------
# Monte Carlo Einbettung
# -------------------------
R <- 10
n <- 100
beta <- c(2, 5, 0, 0, 0)
p <- length(beta)
Sigma <- diag(p)
sigma <- 1

B1 <- 50
K  <- 6
B2 <- 50

# Anteil nahe 1, stark verdichtet im Bereich (0.7..1) und besonders nahe 1
m_frac <- 1 - 10^-seq(0.125, 3, by = 0.10)     # 1 - 10^{-0.5} ... 1 - 10^{-3}
m_grid <- sort(unique(pmax(5, pmin(n, ceiling(n * m_frac)))))
if (!(n %in% m_grid)) m_grid <- sort(unique(c(m_grid, n)))


# Zählstruktur: pro m zählen wir, welches Modell in der Replikation "gewinnt"


# Zählstruktur:
# PE: pro m
win_counts_PE <- setNames(vector("list", length(m_grid)), paste0("m=", m_grid))
for (key in names(win_counts_PE)) win_counts_PE[[key]] <- integer(0)

# BIC: 
win_counts_BIC <- integer(0)

for (r in 1:R) {
  dat_r <- generate_data(n, p, beta, sigma, Sigma)
  
  # Step 1: Screening (Top-K)
  A_screen_r <- step1_screen(dat_r, p, B1, K)
  
  # BIC-Referenz (einmal pro Replikation, nur über Top-K)
  bic_win <- bic_winner_on_screened(A_screen_r, dat_r)
  win_counts_BIC[bic_win] <- if (bic_win %in% names(win_counts_BIC)) win_counts_BIC[bic_win] + 1L else 1L
  
  # Step 2: PE(m) bootstrappen
  for (m_val in m_grid) {
    sel_r <- selection_prob_m_PE(A_screen_r, dat_r, m = m_val, B2 = B2)
    pe_win <- sel_r$model[1]  # max sel_prob
    
    key <- paste0("m=", m_val)
    vec <- win_counts_PE[[key]]
    vec[pe_win] <- if (pe_win %in% names(vec)) vec[pe_win] + 1L else 1L
    win_counts_PE[[key]] <- vec
  }
}

# Ergebnis-Tabelle: Gewinnerwahrscheinlichkeit über Replikationen (je m; PE)
mc_res <- do.call(rbind, lapply(m_grid, function(m_val) {
  key <- paste0("m=", m_val)
  counts <- win_counts_PE[[key]]
  
  data.frame(
    m = m_val,
    model = names(counts),
    wins = as.integer(counts),
    win_prob_over_R = as.numeric(counts) / R,
    row.names = NULL
  )
}))

mc_res <- mc_res %>%
  mutate(
    is_true = (model == "{1,2}"),
    model_lab = ifelse(is_true, "a0 = {1,2}", model),
    m_f = factor(m, levels = sort(unique(m)))
  ) %>%
  arrange(m, desc(win_prob_over_R))


mc_res

mc_bic <- data.frame(
  model = names(win_counts_BIC),
  wins  = as.integer(win_counts_BIC),
  win_prob_over_R = as.numeric(win_counts_BIC) / R,
  row.names = NULL
) %>%
  mutate(
    is_true = (model == "{1,2}"),
    model_lab = ifelse(is_true, "a0 = {1,2}", model)
  ) %>%
  arrange(desc(win_prob_over_R))

mc_bic


mc_bic_full <- mc_res %>%
  distinct(model) %>%
  left_join(mc_bic %>% select(model, win_prob_over_R), by = "model") %>%
  mutate(win_prob_over_R = ifelse(is.na(win_prob_over_R), 0, win_prob_over_R))

# PE-Teil
mc_pe <- mc_res %>%
  transmute(
    m,
    model,
    win_prob_over_R,
    method = "PE"
  )

bic_best_model <- mc_bic$model[1]

mc_bic_rep <- tibble(
  m = sort(unique(mc_res$m)),
  model = bic_best_model
) %>%
  left_join(
    mc_bic %>% select(model, win_prob_over_R),
    by = "model"
  ) %>%
  mutate(method = "BIC")

mc_plot <- bind_rows(mc_pe, mc_bic_rep) %>%
  mutate(
    is_true = (model == "{1,2}"),
    model_lab = ifelse(is_true, "a0 = {1,2}", model),
    m_f = factor(m, levels = sort(unique(m)))
  )


# Linienplot
ggplot(mc_plot,
       aes(x = m,
           y = win_prob_over_R,
           group = interaction(model, method),
           linetype = method,
           color = is_true,
           alpha = is_true)) +
  geom_line() +
  geom_point(size = 1) +
  scale_alpha_manual(values = c("TRUE" = 1, "FALSE" = 0.4)) +
  labs(
    x = "Bootstrap sample size m (n = 100)",
    y = "Selection probability over MC replications",
    title = "Selection probabilities across m",
    subtitle = "PE (m-out-of-n) vs BIC reference; True model: a0 = {1,2}"
  ) +
  theme_minimal(base_size = 13)

# Heatmap prep
mc_pe_hm <- mc_plot %>%
  filter(method == "PE") %>%
  group_by(model) %>%
  mutate(total_prob = sum(win_prob_over_R)) %>%
  ungroup() %>%
  mutate(model_f = fct_reorder(model, total_prob))

ggplot(mc_pe_hm,
       aes(x = m_f, y = model_f, fill = win_prob_over_R)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_tile(data = subset(mc_pe_hm, is_true),
            fill = NA, color = "black", linewidth = 1.0) +
  labs(x = "m", y = "Model", fill = "Selection prob.",
       title = "PE selection probabilities across m",
       subtitle = "True model outlined") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank())

# Gap-Plot
gap_df <- mc_plot %>%
  filter(method=="PE") %>%
  group_by(m) %>%
  summarise(
    p1 = sort(win_prob_over_R, decreasing=TRUE)[1],
    p2 = sort(win_prob_over_R, decreasing=TRUE)[2],
    gap = p1 - p2,
    .groups="drop"
  )

ggplot(gap_df, aes(x=m, y=gap)) +
  geom_line() + geom_point(size=1) +
  labs(x="m", y="p1 - p2",
       title="Separation between best and second-best model") +
  theme_minimal(base_size=13)

# Rank-Trajectory

rank_true <- mc_plot %>%
  filter(method=="PE") %>%
  group_by(m) %>%
  arrange(desc(win_prob_over_R)) %>%
  mutate(rank = row_number()) %>%
  ungroup() %>%
  filter(model=="{1,2}")

ggplot(rank_true, aes(x=m, y=rank)) +
  geom_line() + geom_point(size=1) +
  scale_y_reverse() +
  labs(x="m", y="Rank of true model (1=best)",
       title="Rank trajectory of the true model across m") +
  theme_minimal(base_size=13)
