library(MASS)
library(ggplot2)
library(dplyr)
library(forcats)
library(tidyr)
library(stringr)
library(grid)
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
  
  # Full candidate set (ohne Screening)
  A_full <- all_models(p, k_max = k_max)
  
  # (2) K-grid fürs Screening (Top-K models)
  if (is.null(K_grid)) {
    K_grid <- sort(unique(pmax(2L, pmin(length(A_full), c(5L, 10L, 20L, 50L)))))
    K_grid <- K_grid[K_grid <= length(A_full)]  # safety
  }
  
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

# ------------------------------------------------------------
# 0) Konventionen: Szenarien & Farben
# ------------------------------------------------------------
scenarios_keep <- c("baseline", "weak_signal")

# colorblind-safe palette for K (incl. no screening)
pal_K <- c(
  "no screening" = "#4C4C4C",
  "K=5"  = "#1B9E77",
  "K=10" = "#D95F02",
  "K=20" = "#7570B3",
  "K=50" = "#E7298A"
)

col_true   <- "#0072B2"  # blue
col_comp   <- "#D55E00"  # orange
col_grey55 <- "grey55"

col_area_screen <- "#2C7FB8"
col_area_pe     <- "#F28E2B"

# ------------------------------------------------------------
# 1) Timing: einmal sauber zusammenfassen (Basis für alle Timing-Plots)
# ------------------------------------------------------------
timing_all <- bind_rows(lapply(res_list, "[[", "timing_df")) %>%
  filter(scenario %in% scenarios_keep)

timing_sum <- timing_all %>%
  group_by(scenario, procedure_lab, K, m) %>%
  summarise(
    mean_screen = mean(t_screen),
    mean_bic    = mean(t_bic),
    mean_pe     = mean(t_pe),
    mean_total  = mean(t_total),
    .groups = "drop"
  ) %>%
  mutate(
    # einheitliches K-Label für Legenden und Facets
    K_lab = case_when(
      procedure_lab == "no screening" ~ "no screening",
      !is.na(K) ~ paste0("K=", K),
      TRUE ~ NA_character_
    ),
    m_num = as.numeric(m)
  ) %>%
  filter(!is.na(K_lab))

# Faktor-Reihenfolge für K fixieren (no screening oben, dann numerisch)
K_levels <- c("no screening", paste0("K=", sort(unique(timing_sum$K[!is.na(timing_sum$K)]))))
timing_sum <- timing_sum %>%
  mutate(K_lab = factor(K_lab, levels = K_levels))

# ------------------------------------------------------------
# 2) Plot A: Runtime decomposition (screening + PE; BIC raus)
# ------------------------------------------------------------
timing_comp_long <- timing_sum %>%
  transmute(
    scenario, K_lab, m_num,
    mean_screen, mean_pe,
    mean_total = mean_screen + mean_pe
  ) %>%
  pivot_longer(
    cols = c(mean_screen, mean_pe),
    names_to = "component",
    values_to = "mean_time"
  ) %>%
  mutate(
    component = recode(component,
                       mean_screen = "screening",
                       mean_pe     = "PE (m-out-of-n)")
  ) %>%
  arrange(scenario, K_lab, component, m_num)

timing_total_sp <- timing_comp_long %>%
  distinct(scenario, K_lab, m_num, mean_total)

p_runtime_decomp <- ggplot() +
  geom_hline(yintercept = 6, linetype = "dashed",
             linewidth = 0.6, color = "grey40") +
  geom_area(
    data = timing_comp_long,
    aes(x = m_num, y = mean_time, fill = component, group = component),
    position = "stack",
    alpha = 0.9
  ) +
  geom_line(
    data = timing_total_sp,
    aes(x = m_num, y = mean_total, group = 1),
    linewidth = 1.1,
    color = "black"
  ) +
  scale_fill_manual(values = c("screening" = col_area_screen,
                               "PE (m-out-of-n)" = col_area_pe)) +
  scale_x_continuous(breaks = sort(unique(timing_total_sp$m_num))) +
  labs(
    x = "Resample size m",
    y = "Mean elapsed time [seconds]",
    title = "Runtime decomposition",
    subtitle = paste(
      "Elapsed wall-clock time per Monte Carlo replication.",
      "Black line: total; dashed line: 6-second threshold.",
      sep = "\n"
    )
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text.y = element_text(angle = 0, face = "bold"),
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold")
  ) +
  facet_grid(K_lab ~ scenario, scales = "free_y")

# ------------------------------------------------------------
# 3) Model-selection Plot: true vs strongest competitor + others + BIC ref
# ------------------------------------------------------------
mc_plot_clean <- mc_plot_all %>%
  filter(scenario %in% scenarios_keep) %>%
  filter(procedure == "no_screen" | grepl("^with_screen_K\\d+$", procedure)) %>%
  mutate(
    K_val = ifelse(grepl("^with_screen_K\\d+$", procedure),
                   as.integer(sub("^with_screen_K", "", procedure)),
                   NA_integer_),
    facet_row = case_when(
      procedure == "no_screen" ~ "no screening",
      !is.na(K_val) ~ paste0("K=", K_val),
      TRUE ~ procedure
    )
  )

# Facet-Reihenfolge fixieren
K_levels_mc <- c("no screening", paste0("K=", sort(unique(mc_plot_clean$K_val[!is.na(mc_plot_clean$K_val)]))))
mc_plot_clean <- mc_plot_clean %>%
  mutate(facet_row = factor(facet_row, levels = K_levels_mc))

# strongest competitor: über m aggregiert (PE only, excluding true)
top_comp <- mc_plot_clean %>%
  filter(method == "PE", !is_true) %>%
  group_by(scenario, procedure, facet_row, model) %>%
  summarise(score = sum(win_prob_over_R), .groups = "drop") %>%
  group_by(scenario, procedure, facet_row) %>%
  slice_max(score, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(scenario, procedure, facet_row, comp_model = model)

plot_df <- mc_plot_clean %>%
  left_join(top_comp, by = c("scenario","procedure","facet_row")) %>%
  mutate(
    model_class = case_when(
      method == "PE"  & is_true              ~ "true (PE)",
      method == "BIC" & is_true              ~ "true (BIC ref)",
      method == "PE"  & model == comp_model  ~ "top competitor (PE)",
      method == "PE"                          ~ "other (PE)",
      TRUE                                    ~ "other"
    )
  )

p_true_vs_comp <- ggplot() +
  geom_line(
    data = plot_df %>% filter(model_class == "other (PE)"),
    aes(x = m, y = win_prob_over_R, group = model),
    linewidth = 0.35, alpha = 0.18, color = col_grey55
  ) +
  geom_line(
    data = plot_df %>% filter(model_class == "top competitor (PE)"),
    aes(x = m, y = win_prob_over_R, group = model),
    linewidth = 0.85, color = col_comp
  ) +
  geom_point(
    data = plot_df %>% filter(model_class == "top competitor (PE)"),
    aes(x = m, y = win_prob_over_R),
    size = 1.25, color = col_comp
  ) +
  geom_line(
    data = plot_df %>% filter(model_class == "true (PE)"),
    aes(x = m, y = win_prob_over_R),
    linewidth = 1.25, color = col_true
  ) +
  geom_point(
    data = plot_df %>% filter(model_class == "true (PE)"),
    aes(x = m, y = win_prob_over_R),
    size = 1.35, color = col_true
  ) +
  geom_line(
    data = plot_df %>% filter(model_class == "true (BIC ref)"),
    aes(x = m, y = win_prob_over_R),
    linewidth = 0.95, color = "black"
  ) +
  facet_grid(facet_row ~ scenario) +
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,0.2)) +
  labs(
    x = "Resample size m",
    y = expression(hat(pi)(a[0])~~"(selection rate)"),
    title = "True-model selection rate (with strongest competitor)",
    subtitle = "Blue: true model (PE); orange: strongest competitor; grey \n other models; black: BIC reference"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    strip.text.y = element_text(angle = 0),
    panel.spacing.y = unit(0.7, "lines")
  )

# ------------------------------------------------------------
# 4) Trade-off Plot: selection vs runtime (dual axis; no BIC)
# ------------------------------------------------------------
# sel_df muss existieren: (scenario, m, K_lab, sel) – falls nicht, so erzeugen:
# sel_df <- plot_df %>% filter(model_class == "true (PE)") %>% transmute(scenario, m, K_lab = facet_row, sel = win_prob_over_R)

rt_df <- timing_sum %>%
  transmute(
    scenario,
    m,
    K_lab,
    t_total = mean_screen + mean_pe
  )

t_max <- max(rt_df$t_total, na.rm = TRUE)

p_tradeoff <- ggplot() +
  geom_line(
    data = sel_df,
    aes(x = m, y = sel, color = K_lab),
    linewidth = 1.3
  ) +
  geom_point(
    data = sel_df,
    aes(x = m, y = sel, color = K_lab),
    size = 2.2
  ) +
  geom_line(
    data = rt_df,
    aes(x = m, y = t_total / t_max, color = K_lab),
    linewidth = 0.9,
    linetype = "dashed",
    alpha = 0.6
  ) +
  scale_color_manual(values = pal_K) +
  scale_y_continuous(
    name = expression(hat(pi)(a[0])~"(selection rate)"),
    limits = c(0, 1),
    sec.axis = sec_axis(~ . * t_max, name = "Mean elapsed time [s]")
  ) +
  facet_grid(. ~ scenario) +
  labs(
    x = "Resample size m",
    title = "Accuracy–Runtime Trade-off",
    subtitle = "Solid: selection rate (PE). Dashed: runtime."
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    legend.title = element_blank(),
    panel.grid.minor = element_blank(),
    panel.grid.major.x = element_blank(),
    plot.title = element_text(face = "bold", size = 16),
    strip.text = element_text(face = "bold")
  )

# ------------------------------------------------------------
# 5) Print plots (oder speichern)
# ------------------------------------------------------------
p_runtime_decomp
p_true_vs_comp
p_tradeoff
