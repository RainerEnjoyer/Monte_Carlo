library(MASS)
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)

library(future)
library(future.apply)
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)
# -------------------------
# DGP-Konfiguration
# -------------------------
get_scenario_params <- function(scenario, p = 5) {
  if (scenario == "baseline") {
    list(beta = c(2, 5, 0, 0, 0), sigma = 1, Sigma = diag(p))
  } else if (scenario == "high_noise") {
    list(beta = c(2, 5, 0, 0, 0), sigma = 3, Sigma = diag(p))
  } else if (scenario == "weak_signal") {   # = low signal
    list(beta = c(0.6, 0.8, 0, 0, 0), sigma = 1, Sigma = diag(p))
  } else if (scenario == "correlated") {
    rho <- 0.7
    Sigma <- outer(1:p, 1:p, function(i, j) rho^abs(i - j))
    list(beta = c(2, 5, 0, 0, 0), sigma = 1, Sigma = Sigma)
  } else {
    stop("Unknown scenario: ", scenario)
  }
}

generate_data <- function(n, p, beta, sigma, Sigma, scenario = "baseline") {
  X <- MASS::mvrnorm(n, mu = rep(0, p), Sigma = Sigma)
  colnames(X) <- paste0("X", 1:p)
  eps <- rnorm(n, 0, sigma)
  y <- as.vector(X %*% beta + eps)
  data.frame(y = y, X)
}

# -------------------------
# Hilfsfunktion: alle nicht-leeren Modelle (Full search space)
# -------------------------
all_models <- function(p) {
  unlist(lapply(1:p, function(k) combn(1:p, k, simplify = FALSE)), recursive = FALSE)
}

# -------------------------
# Step 1: Screening (Top-K) via grobem Bootstrap-PE
# -------------------------
step1_screen <- function(dat, p, B1, K) {
  n <- nrow(dat)
  
  # alle nicht-leeren Modelle
  A <- all_models(p)
  
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
# Gewinner pro Bootstrap-Run = Modell mit minimalem PE (MSE) auf OOB
# -------------------------
selection_prob_m_PE <- function(A_candidates, dat, m, B2) {
  n <- nrow(dat)
  
  models_chr <- sapply(A_candidates, function(a) paste0("{", paste(a, collapse = ","), "}"))
  wins_PE <- setNames(integer(length(A_candidates)), models_chr)
  
  forms <- lapply(A_candidates, function(a) {
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
    
    pe_vec <- numeric(length(A_candidates))
    
    for (j in seq_along(A_candidates)) {
      fit_b <- lm(forms[[j]], data = dat_b)
      
      yhat_oob <- predict(fit_b, newdata = dat_oob)
      mse_oob  <- mean((dat_oob$y - yhat_oob)^2)
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

bic_winner_on_candidates <- function(A_candidates, dat) {
  forms <- lapply(A_candidates, function(a) {
    rhs <- paste0("X", a, collapse = " + ")
    as.formula(paste("y ~", rhs))
  })
  bics <- sapply(forms, function(f) BIC(lm(f, data = dat)))
  j_star <- which.min(bics)
  paste0("{", paste(A_candidates[[j_star]], collapse = ","), "}")
}

# -------------------------
# Monte Carlo Einbettung
# -------------------------
run_scenario <- function(scenario,
                         n = 100, p = 5,
                         R = 100, B1 = 50, K = 6, B2 = 100,
                         m_grid = NULL,
                         seed = 1) {
  
  set.seed(seed)
  
  par <- get_scenario_params(scenario, p = p)
  beta  <- par$beta
  sigma <- par$sigma
  Sigma <- par$Sigma
  
  if (is.null(m_grid)) {
    m_frac <- 1 - 10^-seq(0.125, 3, by = 0.10)
    m_grid <- sort(unique(pmax(5, pmin(n, ceiling(n * m_frac)))))
    if (!(n %in% m_grid)) m_grid <- sort(unique(c(m_grid, n)))
  }
  
  # Full candidate set (ohne Screening)
  A_full <- all_models(p)
  
  # zwei Verfahren
  procedures_raw <- c("with_screen", "no_screen")
  procedure_factor <- factor(
    procedures_raw,
    levels = c("with_screen", "no_screen"),
    labels = c("with screening", "no screening")
  )
  
  # Zählstrukturen getrennt nach procedure
  win_counts_PE <- setNames(vector("list", length(procedures_raw)), procedures_raw)
  for (proc in procedures_raw) {
    win_counts_PE[[proc]] <- setNames(vector("list", length(m_grid)), paste0("m=", m_grid))
    for (key in names(win_counts_PE[[proc]])) win_counts_PE[[proc]][[key]] <- integer(0)
  }
  
  win_counts_BIC <- setNames(vector("list", length(procedures_raw)), procedures_raw)
  for (proc in procedures_raw) win_counts_BIC[[proc]] <- integer(0)
  
  for (r in 1:R) {
    dat_r <- generate_data(n, p, beta, sigma, Sigma, scenario = scenario)
    
    # Screening nur einmal berechnen
    A_screen_r <- step1_screen(dat_r, p, B1, K)
    
    for (proc in procedures_raw) {
      
      A_candidates <- if (proc == "with_screen") A_screen_r else A_full
      
      # BIC einmal pro Replikation & procedure
      bic_win <- bic_winner_on_candidates(A_candidates, dat_r)
      win_counts_BIC[[proc]][bic_win] <- if (bic_win %in% names(win_counts_BIC[[proc]]))
        win_counts_BIC[[proc]][bic_win] + 1L else 1L
      
      # PE über m
      for (m_val in m_grid) {
        sel_r <- selection_prob_m_PE(A_candidates, dat_r, m = m_val, B2 = B2)
        pe_win <- sel_r$model[1]  # max sel_prob
        
        key <- paste0("m=", m_val)
        vec <- win_counts_PE[[proc]][[key]]
        vec[pe_win] <- if (pe_win %in% names(vec)) vec[pe_win] + 1L else 1L
        win_counts_PE[[proc]][[key]] <- vec
      }
    }
  }
  
  # mc_res (PE) für beide procedures
  mc_res <- do.call(rbind, lapply(procedures_raw, function(proc) {
    do.call(rbind, lapply(m_grid, function(m_val) {
      key <- paste0("m=", m_val)
      counts <- win_counts_PE[[proc]][[key]]
      data.frame(
        scenario = scenario,
        procedure = proc,
        m = m_val,
        model = names(counts),
        wins = as.integer(counts),
        win_prob_over_R = as.numeric(counts) / R,
        row.names = NULL
      )
    }))
  })) %>%
    mutate(
      is_true = (model == "{1,2}"),
      model_lab = ifelse(is_true, "a0 = {1,2}", model),
      m_f = factor(m, levels = sort(unique(m))),
      procedure = factor(
        procedure,
        levels = c("with_screen", "no_screen"),
        labels = c("with screening", "no screening")
      )
    ) %>%
    arrange(procedure, m, desc(win_prob_over_R))
  
  # BIC table für beide procedures
  mc_bic <- do.call(rbind, lapply(procedures_raw, function(proc) {
    data.frame(
      scenario = scenario,
      procedure = proc,
      model = names(win_counts_BIC[[proc]]),
      wins  = as.integer(win_counts_BIC[[proc]]),
      win_prob_over_R = as.numeric(win_counts_BIC[[proc]]) / R,
      row.names = NULL
    )
  })) %>%
    mutate(
      is_true = (model == "{1,2}"),
      model_lab = ifelse(is_true, "a0 = {1,2}", model),
      procedure = factor(
        procedure,
        levels = c("with_screen", "no_screen"),
        labels = c("with screening", "no screening")
      )
    ) %>%
    group_by(scenario, procedure) %>%
    arrange(desc(win_prob_over_R), .by_group = TRUE) %>%
    ungroup()
  
  # Für Plot: PE-Linien + BIC-Referenz als horizontale Linie pro (Szenario, procedure)
  mc_bic_ref <- mc_bic %>%
    group_by(scenario, procedure) %>%
    slice_max(order_by = win_prob_over_R, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    select(scenario, procedure, model, win_prob_over_R) %>%
    tidyr::expand_grid(m = sort(unique(mc_res$m))) %>%
    mutate(
      method = "BIC",
      is_true = (model == "{1,2}"),
      model_lab = ifelse(is_true, "a0 = {1,2}", model),
      m_f = factor(m, levels = sort(unique(mc_res$m)))
    )
  
  mc_pe <- mc_res %>%
    transmute(
      scenario, procedure, m, model, win_prob_over_R,
      method = "PE", is_true, model_lab, m_f
    )
  
  mc_plot <- bind_rows(mc_pe, mc_bic_ref)
  
  # Heatmap-Prep (nur PE)
  mc_pe_hm <- mc_plot %>%
    filter(method == "PE") %>%
    group_by(scenario, procedure, model) %>%
    mutate(total_prob = sum(win_prob_over_R)) %>%
    ungroup() %>%
    group_by(scenario, procedure) %>%
    mutate(model_f = fct_reorder(model, total_prob)) %>%
    ungroup()
  
  # Gap-Plot (nur PE)
  gap_df <- mc_plot %>%
    filter(method == "PE") %>%
    group_by(scenario, procedure, m) %>%
    summarise(
      p1 = sort(win_prob_over_R, decreasing = TRUE)[1],
      p2 = sort(win_prob_over_R, decreasing = TRUE)[2],
      gap = p1 - p2,
      .groups = "drop"
    )
  
  # Rank trajectory (true model, PE)
  rank_true <- mc_plot %>%
    filter(method == "PE") %>%
    group_by(scenario, procedure, m) %>%
    arrange(desc(win_prob_over_R), .by_group = TRUE) %>%
    mutate(rank = row_number()) %>%
    ungroup() %>%
    filter(model == "{1,2}")
  
  list(
    mc_plot = mc_plot,
    mc_pe_hm = mc_pe_hm,
    gap_df = gap_df,
    rank_true = rank_true,
    mc_bic = mc_bic
  )
}

# -------------------------
# Run over scenarios
# -------------------------
scenarios <- c("baseline", "high_noise", "weak_signal", "correlated")

res_list <- future.apply::future_lapply(seq_along(scenarios), function(i) {
  run_scenario(scenarios[i], seed = 1 + i)
}, future.seed = TRUE)

future::plan(sequential)

mc_plot_all   <- bind_rows(lapply(res_list, [[, "mc_plot"))
mc_pe_hm_all  <- bind_rows(lapply(res_list, [[, "mc_pe_hm"))
gap_df_all    <- bind_rows(lapply(res_list, [[, "gap_df"))
rank_true_all <- bind_rows(lapply(res_list, [[, "rank_true"))
mc_bic_all    <- bind_rows(lapply(res_list, [[, "mc_bic"))

# -------------------------
# Plots: scenario × procedure (nebeneinander per facet_grid)
# -------------------------

ggplot(mc_plot_all,
       aes(x = m, y = win_prob_over_R,
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
  theme_minimal(base_size = 13) +
  theme(legend.position = "none") +
  facet_grid(procedure ~ scenario)

ggplot(mc_pe_hm_all,
       aes(x = m_f, y = model_f, fill = win_prob_over_R)) +
  geom_tile(color = "white", linewidth = 0.35) +
  geom_tile(data = subset(mc_pe_hm_all, is_true),
            fill = NA, color = "black", linewidth = 1.0) +
  labs(x = "m", y = "Model", fill = "Selection prob.",
       title = "PE selection probabilities across m",
       subtitle = "True model outlined") +
  theme_minimal(base_size = 13) +
  theme(panel.grid = element_blank()) +
  facet_grid(procedure ~ scenario)

ggplot(gap_df_all, aes(x = m, y = gap)) +
  geom_line() + geom_point(size = 1) +
  labs(x = "m", y = "p1 - p2",
       title = "Separation between best and second-best model (PE)") +
  theme_minimal(base_size = 13) +
  facet_grid(procedure ~ scenario)

ggplot(rank_true_all, aes(x = m, y = rank)) +
  geom_line() + geom_point(size = 1) +
  scale_y_reverse() +
  labs(x = "m", y = "Rank of true model (1=best)",
       title = "Rank trajectory of the true model (PE)") +
  theme_minimal(base_size = 13) +
  facet_grid(procedure ~ scenario)