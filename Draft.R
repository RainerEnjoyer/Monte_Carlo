library(MASS)
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)
library(stringr)

library(future)
library(future.apply)
workers <- max(1, parallel::detectCores() - 1)
plan(multisession, workers = workers)
# -------------------------
# DGP-Konfiguration
# -------------------------
get_scenario_params <- function(scenario, p = 9) {
  
  make_beta <- function(b1, b2) {
    beta <- numeric(p)
    beta[1] <- b1
    beta[2] <- b2
    beta
  }
  
  if (scenario == "baseline") {
    list(beta = make_beta(2, 5), sigma = 1, Sigma = diag(p))
  } else if (scenario == "weak_signal") {
    list(beta = make_beta(0.6, 0.8), sigma = 1, Sigma = diag(p))
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
all_models <- function(p, k_max = p) {
  k_max <- min(p, k_max)
  unlist(lapply(1:k_max, function(k) combn(1:p, k, simplify = FALSE)),
         recursive = FALSE)
}

# -------------------------
# Step 1: Screening (Top-K) via grobem Bootstrap-PE
# -------------------------
step1_screen <- function(dat, p, B1, K, k_max) {
  n <- nrow(dat)
  A <- all_models(p, k_max)
  
  boot_pe <- function(a) {
    form <- as.formula(paste("y ~", paste0("X", a, collapse = " + ")))
    pe <- numeric(B1)
    
    for (b in 1:B1) {
      idx <- sample.int(n, n, replace = TRUE)
      idx_oob <- setdiff(1:n, unique(idx))
      if (length(idx_oob) < 5) { pe[b] <- NA; next }
      fit <- lm(form, data = dat[idx, ])
      yhat <- predict(fit, newdata = dat[idx_oob, ])
      pe[b] <- mean((dat$y[idx_oob] - yhat)^2)
    }
    mean(pe, na.rm = TRUE)
  }
  
  PE <- sapply(A, boot_pe)
  A[order(PE)[1:K]]
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
                         n = 100, p = 9,
                         R = 100, B1 = 20, B2 = 100,
                         m_grid = NULL,
                         K_grid = NULL,
                         k_max = p,
                         seed = 1) {
  
  set.seed(seed)
  
  par <- get_scenario_params(scenario, p = p)
  beta  <- par$beta
  sigma <- par$sigma
  Sigma <- par$Sigma
  
  true_set <- which(beta != 0)
  true_lab <- paste0("{", paste(true_set, collapse = ","), "}")
  
  # -------------------------
  # (1) m-grid: 
  # -------------------------
  if (is.null(m_grid)) {
    m_frac <- c(0.25, 0.4, 0.65, 0.8, 0.9, 1.00)
    m_grid <- sort(unique(pmax(5L, pmin(n, as.integer(ceiling(n * m_frac))))))
    if (!(n %in% m_grid)) m_grid <- sort(unique(c(m_grid, n)))
  }
  
  # -------------------------
  # (2) K-grid fürs Screening
  # -------------------------
  if (is.null(K_grid)) {
    K_grid <- sort(unique(pmax(2L, pmin(p, c(2L, 3L, 5L, 7L, as.integer(ceiling(p/2)))))))
  }
  
  # Full candidate set (ohne Screening)
  A_full <- all_models(p, k_max = k_max)
  
  # Procedures: "no_screen" + mehrere "with_screen" (je K)
  procedures_raw <- c(paste0("with_screen_K", K_grid), "no_screen")
  
  # ---- Timing container (Sekunden) ----
  timing_rows <- vector("list", R * length(procedures_raw) * length(m_grid))
  row_i <- 0L
  
  # Zählstrukturen getrennt nach procedure
  win_counts_PE <- setNames(vector("list", length(procedures_raw)), procedures_raw)
  for (proc in procedures_raw) {
    win_counts_PE[[proc]] <- setNames(vector("list", length(m_grid)), paste0("m=", m_grid))
    for (key in names(win_counts_PE[[proc]])) win_counts_PE[[proc]][[key]] <- integer(0)
  }
  
  win_counts_BIC <- setNames(vector("list", length(procedures_raw)), procedures_raw)
  for (proc in procedures_raw) win_counts_BIC[[proc]] <- integer(0)
  
  # -------------------------
  # Monte Carlo
  # -------------------------
  for (r in 1:R) {
    dat_r <- generate_data(n, p, beta, sigma, Sigma, scenario = scenario)
    
    # Screening für alle K nur einmal pro Replikation berechnen
    A_screen_by_K <- setNames(vector("list", length(K_grid)), paste0("K", K_grid))
    t_screen_by_K <- setNames(numeric(length(K_grid)), paste0("K", K_grid))
    
    for (K_val in K_grid) {
      keyK <- paste0("K", K_val)
      t0 <- proc.time()[["elapsed"]]
      A_screen_by_K[[keyK]] <- step1_screen(dat_r, p, B1, K = K_val, k_max = k_max)
      t_screen_by_K[[keyK]] <- proc.time()[["elapsed"]] - t0
    }
    
    for (proc in procedures_raw) {
      
      is_screen <- proc != "no_screen"
      K_val <- if (is_screen) as.integer(sub("^with_screen_K", "", proc)) else NA_integer_
      keyK <- if (is_screen) paste0("K", K_val) else NA_character_
      
      A_candidates <- if (!is_screen) {
        A_full
      } else {
        A_screen_by_K[[keyK]]
      }
      
      t_screen <- if (is_screen) t_screen_by_K[[keyK]] else 0
      
      # BIC einmal pro Replikation & procedure
      t0_bic <- proc.time()[["elapsed"]]
      bic_win <- bic_winner_on_candidates(A_candidates, dat_r)
      t_bic  <- proc.time()[["elapsed"]] - t0_bic
      
      win_counts_BIC[[proc]][bic_win] <- if (bic_win %in% names(win_counts_BIC[[proc]]))
        win_counts_BIC[[proc]][bic_win] + 1L else 1L
      
      # PE über m
      for (m_val in m_grid) {
        t0_pe <- proc.time()[["elapsed"]]
        sel_r <- selection_prob_m_PE(A_candidates, dat_r, m = m_val, B2 = B2)
        t_pe  <- proc.time()[["elapsed"]] - t0_pe
        
        pe_win <- sel_r$model[1]
        
        key <- paste0("m=", m_val)
        vec <- win_counts_PE[[proc]][[key]]
        vec[pe_win] <- if (pe_win %in% names(vec)) vec[pe_win] + 1L else 1L
        win_counts_PE[[proc]][[key]] <- vec
        
        # ---- record timing row (pro rep x proc x m) ----
        row_i <- row_i + 1L
        timing_rows[[row_i]] <- data.frame(
          scenario = scenario,
          rep = r,
          procedure = proc,
          K = K_val,
          m = m_val,
          t_screen = t_screen,   # identisch über m innerhalb (r,K)
          t_bic = t_bic,         # identisch über m innerhalb (r,proc)
          t_pe = t_pe,           # variiert über m
          t_total = t_screen + t_bic + t_pe,
          row.names = NULL
        )
      }
    }
  }
  timing_df <- dplyr::bind_rows(timing_rows[seq_len(row_i)]) %>%
    dplyr::mutate(
      procedure_lab = dplyr::case_when(
        procedure == "no_screen" ~ "no screening",
        grepl("^with_screen_K", procedure) ~ paste0("with screening (K=", sub("^with_screen_K", "", procedure), ")"),
        TRUE ~ procedure
      )
    )
  
  # -------------------------
  # mc_res (PE) für alle procedures
  # -------------------------
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
    dplyr::mutate(
      is_true = (model == true_lab),
      model_lab = ifelse(is_true, paste0("a0 = ", true_lab), model),
      m_f = factor(m, levels = sort(unique(m))),
      procedure_lab = dplyr::case_when(
        procedure == "no_screen" ~ "no screening",
        grepl("^with_screen_K", procedure) ~ paste0("with screening (K=", sub("^with_screen_K", "", procedure), ")"),
        TRUE ~ procedure
      )
    ) %>%
    dplyr::arrange(procedure_lab, m, dplyr::desc(win_prob_over_R))
  
  # -------------------------
  # BIC table für alle procedures
  # -------------------------
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
    dplyr::mutate(
      is_true = (model == true_lab),
      model_lab = ifelse(is_true, paste0("a0 = ", true_lab), model),
      procedure_lab = dplyr::case_when(
        procedure == "no_screen" ~ "no screening",
        grepl("^with_screen_K", procedure) ~ paste0("with screening (K=", sub("^with_screen_K", "", procedure), ")"),
        TRUE ~ procedure
      )
    ) %>%
    dplyr::group_by(scenario, procedure, procedure_lab) %>%
    dplyr::arrange(dplyr::desc(win_prob_over_R), .by_group = TRUE) %>%
    dplyr::ungroup()
  
  # -------------------------
  # Plot-Prep: PE-Linien + BIC-Referenz als horizontale Linie pro (Szenario, procedure)
  # -------------------------
  mc_bic_ref <- mc_bic %>%
    dplyr::group_by(scenario, procedure, procedure_lab) %>%
    dplyr::slice_max(order_by = win_prob_over_R, n = 1, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::select(scenario, procedure, procedure_lab, model, win_prob_over_R) %>%
    tidyr::expand_grid(m = sort(unique(mc_res$m))) %>%
    dplyr::mutate(
      method = "BIC",
      is_true = (model == true_lab),
      model_lab = ifelse(is_true, paste0("a0 = ", true_lab), model),
      m_f = factor(m, levels = sort(unique(mc_res$m)))
    )
  
  mc_pe <- mc_res %>%
    dplyr::transmute(
      scenario, procedure, procedure_lab, m, model, win_prob_over_R,
      method = "PE", is_true, model_lab, m_f
    )
  
  mc_plot <- dplyr::bind_rows(mc_pe, mc_bic_ref)
  
  # Heatmap-Prep (nur PE)
  mc_pe_hm <- mc_plot %>%
    dplyr::filter(method == "PE") %>%
    dplyr::group_by(scenario, procedure_lab, model) %>%
    dplyr::mutate(total_prob = sum(win_prob_over_R)) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(scenario, procedure_lab) %>%
    dplyr::mutate(model_f = forcats::fct_reorder(model, total_prob)) %>%
    dplyr::ungroup()
  
  # Gap-Plot (nur PE)
  gap_df <- mc_plot %>%
    dplyr::filter(method == "PE") %>%
    dplyr::group_by(scenario, procedure_lab, m) %>%
    dplyr::summarise(
      p1 = sort(win_prob_over_R, decreasing = TRUE)[1],
      p2 = sort(win_prob_over_R, decreasing = TRUE)[2],
      gap = p1 - p2,
      .groups = "drop"
    )
  
  # Rank trajectory (true model, PE)
  rank_true <- mc_plot %>%
    dplyr::filter(method == "PE") %>%
    dplyr::group_by(scenario, procedure_lab, m) %>%
    dplyr::arrange(dplyr::desc(win_prob_over_R), .by_group = TRUE) %>%
    dplyr::mutate(rank = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::filter(model == true_lab)
  
  list(
    mc_plot = mc_plot,
    mc_pe_hm = mc_pe_hm,
    gap_df = gap_df,
    rank_true = rank_true,
    mc_bic = mc_bic,
    timing_df = timing_df
  )
}

# -------------------------
# Run over scenarios
# -------------------------
scenarios <- c("baseline", "weak_signal")

R0  <- 5; B10 <- 20; B20 <- 50

t0 <- system.time(
  invisible(run_scenario("baseline",
                         n = 100, p = 10,
                         R = R0, B1 = B10, B2 = B20,
                         k_max = 5, seed = 1))
)["elapsed"]

R  <- 50; B1 <- 20; B2 <- 100

t_pred <- t0 * (R/R0) * (B1/B10) * (B2/B20)
c(seconds = t_pred, minutes = t_pred/60, hours = t_pred/3600)

res_list <- future.apply::future_lapply(seq_along(scenarios), function(i) {
  run_scenario(scenarios[i], seed = 1 + i)
}, future.seed = TRUE)

future::plan(sequential)

mc_plot_all   <- bind_rows(lapply(res_list, "[[", "mc_plot"))
mc_pe_hm_all  <- bind_rows(lapply(res_list, "[[", "mc_pe_hm"))
gap_df_all    <- bind_rows(lapply(res_list, "[[", "gap_df"))
rank_true_all <- bind_rows(lapply(res_list, "[[", "rank_true"))
mc_bic_all    <- bind_rows(lapply(res_list, "[[", "mc_bic"))

timing_all <- bind_rows(lapply(res_list, "[[", "timing_df"))

timing_sum <- timing_all %>%
  group_by(scenario, procedure_lab, K, m) %>%
  summarise(
    mean_total = mean(t_total),
    mean_screen = mean(t_screen),
    mean_bic = mean(t_bic),
    mean_pe = mean(t_pe),
    .groups = "drop"
  )


# nur true model, beide Methoden
true_df <- mc_plot_all %>%
  filter(is_true) %>%
  mutate(
    K_val = ifelse(grepl("^with_screen_K", procedure),
                   as.integer(sub("^with_screen_K", "", procedure)),
                   NA_integer_),
    facet_row = case_when(
      procedure == "no_screen" ~ "no screening",
      !is.na(K_val) ~ paste0("screening (K=", K_val, ")"),
      TRUE ~ procedure
    )
  )

# numerische Reihenfolge: no_screen oben, dann K aufsteigend
K_levels <- sort(unique(true_df$K_val[!is.na(true_df$K_val)]))
true_df$facet_row <- factor(true_df$facet_row,
                            levels = c("no screening", paste0("screening (K=", K_levels, ")")))

ggplot(true_df, aes(x = m, y = win_prob_over_R, linetype = method, group = method)) +
  geom_line(linewidth = 0.7) +
  geom_point(size = 1.2) +
  labs(
    x = "Resample size m",
    y = expression(hat(pi)(a[0])~~"(true-model selection rate)"),
    title = "How often is the true model selected?",
    subtitle = "PE (m-out-of-n) across m; BIC shown as reference (constant in m)"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "top") +
  facet_grid(facet_row ~ scenario)

timing_all <- bind_rows(lapply(res_list, "[[", "timing_df"))

timing_sum <- timing_all %>%
  group_by(scenario, procedure_lab, K, m) %>%
  summarise(
    mean_screen = mean(t_screen),
    mean_bic    = mean(t_bic),
    mean_pe     = mean(t_pe),
    .groups = "drop"
  ) %>%
  mutate(
    K_lab = case_when(
      procedure_lab == "no screening" ~ "no screening",
      !is.na(K) ~ paste0("screening (K=", K, ")"),
      TRUE ~ procedure_lab
    ),
    m_f = factor(m, levels = sort(unique(m)))
  )

timing_long <- timing_sum %>%
  tidyr::pivot_longer(
    cols = c(mean_screen, mean_bic, mean_pe),
    names_to = "component",
    values_to = "mean_time"
  ) %>%
  mutate(
    component = recode(component,
                       mean_screen = "screening",
                       mean_bic    = "BIC",
                       mean_pe     = "PE (m-out-of-n)")
  )

ggplot(timing_long, aes(x = m_f, y = mean_time, fill = component)) +
  geom_col() +
  labs(
    x = "m",
    y = "Mean time [s]",
    title = "Runtime decomposition",
    subtitle = "Stacked: screening + BIC + PE"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text.y = element_text(angle = 0),
        legend.position = "top") +
  facet_grid(K_lab ~ scenario, scales = "free_y")