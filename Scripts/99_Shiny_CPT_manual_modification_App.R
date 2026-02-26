#############################################################################################
#                                                                                           #
#                                                                                           #
#####             Simple Shiny interface to modify CPTs with mouse cursors              #####
#                                                                                           #
#                                                                                           #
#############################################################################################

library(shiny)
library(DT)
library(ggplot2)
library(tidyr)
library(dplyr)

# =========================
# Utility functions
# =========================

`%||%` <- function(a, b) if (!is.null(a)) a else b

detect_cpt_type <- function(df) {
  if (all(c("State", "Probability") %in% names(df))) {
    return("no_parents")
  } else {
    return("with_parents")
  }
}

# This gets the index of the column, not the column name itself...
get_state_columns <- function(df) {
  res <- grep("^P_", names(df))
  if (length(res) == 0 && "Probability" %in% names(df)) {
    res <- which(names(df) == "Probability")
  }
  res
}

PROB_STEP   <- 0.01
PROB_MIN    <- PROB_STEP                 # enforce min prob 0.01
TOTAL_UNITS <- round(1 / PROB_STEP)      # 100

# Snap any probability vector to the grid with min PROB_MIN and exact sum 1.
# Useful to "clean" imported CPTs that contain 1/3, etc.
snap_probs_to_grid <- function(probs, step = PROB_STEP, min_prob = PROB_MIN) {
  probs <- as.numeric(probs)
  K <- length(probs)
  if (K == 0) return(probs)
  
  # Replace NA/Inf with min_prob
  probs[!is.finite(probs)] <- min_prob
  
  # Clamp to [min_prob, 1]
  probs <- pmax(probs, min_prob)
  probs <- pmin(probs, 1)
  
  # Convert to units, enforce min 1 unit each, enforce sum=TOTAL_UNITS
  U <- round(probs / step)
  U[U < 1] <- 1
  
  # If sum too big, remove units from largest entries
  while (sum(U) > TOTAL_UNITS) {
    j <- which.max(U)
    if (U[j] <= 1) break
    U[j] <- U[j] - 1L
  }
  
  # If sum too small, add units to largest entry (or any entry)
  while (sum(U) < TOTAL_UNITS) {
    j <- which.max(U)
    U[j] <- U[j] + 1L
  }
  
  U * step
}

# The renormalise function for keeping 0-1 between states when modifying a cursor
# Enforces:
#   - grid 0.01
#   - min prob 0.01 for every unlocked state
#   - sum exactly 1.00 on grid
renormalise <- function(probs, locked, changed, new_value, step = PROB_STEP) {
  
  n <- length(probs)
  locked <- as.logical(locked)
  total_units <- round(1 / step)
  
  probs <- as.numeric(probs)
  probs[!is.finite(probs)] <- step
  
  # integer units (grid)
  units <- round(probs / step)
  
  idx_all <- seq_len(n)
  locked_others <- locked & idx_all != changed
  free_idx <- which(!locked & idx_all != changed)
  
  # Minimum units: 1 for any unlocked state (including changed if not locked)
  min_units <- rep(0L, n)
  min_units[which(!locked)] <- 1L
  
  used_units_locked <- sum(units[locked_others])
  
  # Minimum required for everyone except changed
  min_required_except_changed <- sum(min_units[idx_all != changed])
  
  # Max feasible for changed
  max_units_changed <- total_units - used_units_locked - min_required_except_changed
  max_units_changed <- max(max_units_changed, min_units[changed]) # ensure >= min
  
  # Clamp changed to [min, max]
  target_units_changed <- round(new_value / step)
  target_units_changed <- max(target_units_changed, min_units[changed])
  target_units_changed <- min(target_units_changed, max_units_changed)
  
  units[changed] <- target_units_changed
  
  # Start free states at min
  units[free_idx] <- min_units[free_idx]
  
  # Remaining units to distribute among free states
  remaining_units <- total_units - used_units_locked - units[changed] - sum(units[free_idx])
  
  if (length(free_idx) > 0 && remaining_units > 0) {
    # weights above minimum from old probs
    w <- pmax(round(probs[free_idx] / step) - min_units[free_idx], 0L)
    sw <- sum(w)
    
    if (sw == 0) {
      base <- remaining_units %/% length(free_idx)
      extra <- remaining_units %% length(free_idx)
      units[free_idx] <- units[free_idx] + base
      if (extra > 0) units[free_idx][seq_len(extra)] <- units[free_idx][seq_len(extra)] + 1L
    } else {
      alloc <- floor(w / sw * remaining_units)
      rem <- remaining_units - sum(alloc)
      units[free_idx] <- units[free_idx] + alloc
      if (rem > 0) {
        ord <- order(-w)
        units[free_idx][ord[seq_len(rem)]] <- units[free_idx][ord[seq_len(rem)]] + 1L
      }
    }
  }
  
  # Final exact sum correction (adjust a free state if possible, otherwise changed)
  diff <- total_units - sum(units)
  if (diff != 0) {
    if (length(free_idx) > 0) {
      j <- free_idx[1]
      units[j] <- units[j] + diff
      if (units[j] < min_units[j]) {
        # fallback: push back to min and adjust changed
        d2 <- min_units[j] - units[j]
        units[j] <- min_units[j]
        units[changed] <- units[changed] - d2
      }
    } else {
      units[changed] <- units[changed] + diff
    }
  }
  
  # Guarantee mins again
  units[which(!locked)] <- pmax(units[which(!locked)], 1L)
  
  # Guarantee sum again
  diff2 <- total_units - sum(units)
  if (diff2 != 0) {
    j <- which(!locked)[1]
    units[j] <- units[j] + diff2
  }
  
  units * step
}

# Monotone CPT generation
.clean_txt <- function(x) trimws(gsub("_+", " ", as.character(x)))

.default_scores <- function(levels_vec) {
  m <- length(levels_vec)
  seq(-(m - 1) / 2, (m - 1) / 2, length.out = m)
}

infer_cpt_structure <- function(df) {
  child_cols <- grep("^P_", names(df))
  if (length(child_cols) == 0) stop("Expected P_* columns for a with-parents CPT.")
  parent_cols <- setdiff(seq_along(df), child_cols)
  child_states <- sub("^P_", "", names(df)[child_cols])
  list(parent_cols = parent_cols, child_cols = child_cols, child_states = child_states)
}

.probs_ordinal_logit <- function(score, K, temperature = 1, cut_span = 2.0) {
  temperature <- max(1e-6, temperature)
  theta <- seq(-cut_span, cut_span, length.out = K - 1)
  mu <- score / temperature
  cum <- plogis(theta - mu)
  
  p <- numeric(K)
  p[1] <- cum[1]
  if (K > 2) for (k in 2:(K - 1)) p[k] <- cum[k] - cum[k - 1]
  p[K] <- 1 - cum[K - 1]
  pmax(p, 0)
}

# Gaussian-center with decoupled saturation (center_temp) and smoothing width (sigma)
.probs_gaussian_center_cs <- function(score, K, center_temp = 1, sigma = 1) {
  center_temp <- max(1e-6, center_temp)
  sigma <- max(1e-6, sigma)
  
  # center position (1..K), saturates for large |score|
  center <- 1 + (K - 1) * plogis(score / center_temp)
  
  idx <- seq_len(K)
  w <- exp(-0.5 * ((idx - center) / sigma)^2)
  w / sum(w)
}

build_monotone_cpt <- function(df,
                               parent_scores = NULL,
                               parent_weights = NULL,
                               method = c("ordinal_logit", "gaussian_center"),
                               temperature = 1,
                               center_temp = 1,
                               sigma = 1,
                               floor = 0.01) {
  
  method <- match.arg(method)
  st <- infer_cpt_structure(df)
  parent_cols <- st$parent_cols
  child_cols  <- st$child_cols
  K <- length(child_cols)
  
  parent_names <- names(df)[parent_cols]
  
  if (is.null(parent_weights)) {
    parent_weights <- setNames(rep(1, length(parent_names)), parent_names)
  } else {
    missing <- setdiff(parent_names, names(parent_weights))
    if (length(missing) > 0) parent_weights[missing] <- 1
    parent_weights <- parent_weights[parent_names]
  }
  
  if (is.null(parent_scores)) {
    parent_scores <- lapply(parent_names, function(pn) {
      v <- df[[pn]]
      levs <- if (is.factor(v)) levels(v) else unique(as.character(v))
      levs <- as.character(levs)
      setNames(.default_scores(levs), levs)
    })
    names(parent_scores) <- parent_names
  } else {
    missing <- setdiff(parent_names, names(parent_scores))
    if (length(missing) > 0) {
      for (pn in missing) {
        v <- df[[pn]]
        levs <- if (is.factor(v)) levels(v) else unique(as.character(v))
        parent_scores[[pn]] <- setNames(.default_scores(levs), as.character(levs))
      }
    }
    parent_scores <- parent_scores[parent_names]
  }
  
  total_score <- numeric(nrow(df))
  for (pn in parent_names) {
    w <- as.numeric(parent_weights[pn])
    mp <- parent_scores[[pn]]
    vals <- as.character(df[[pn]])
    sc <- mp[vals]
    sc[is.na(sc)] <- 0
    total_score <- total_score + w * as.numeric(sc)
  }
  
  P <- matrix(0, nrow(df), K)
  for (i in seq_len(nrow(df))) {
    s <- total_score[i]
    p <- if (method == "ordinal_logit") {
      .probs_ordinal_logit(s, K, temperature = temperature)
    } else {
      .probs_gaussian_center_cs(s, K, center_temp = center_temp, sigma = sigma)
    }
    p <- pmax(p, 0)
    p <- p + floor
    p <- p / sum(p)
    
    # IMPORTANT: snap monotone output to your grid + min=0.01
    p <- snap_probs_to_grid(p, step = PROB_STEP, min_prob = PROB_MIN)
    
    P[i, ] <- p
  }
  
  out <- df
  out[, child_cols] <- P
  out
}

# =========================
# UI
# =========================

run_cpt_app <- function(pathProCpt,
                        folder = "01_CPT",
                        cpt_patterns = NULL) {
  
  ui <- fluidPage(
    
    titlePanel("Bayesian Network CPT Editor"),
    
    sidebarLayout(
      sidebarPanel(
        selectInput("cpt_file", "Select CPT", choices = NULL),
        tags$h4("Monotone CPT generator"),
        checkboxInput("show_contours", "Show contours", value = FALSE),
        uiOutput("mono_weights_ui"),
        uiOutput("mono_scores_ui"),
        selectInput("mono_method", "Method",
                    choices = c("Ordinal logit" = "ordinal_logit",
                                "Gaussian center" = "gaussian_center"),
                    selected = "gaussian_center"),
        sliderInput("mono_center_temp", "Center saturation (lower = stronger extremes)",
                    min = 0.1, max = 5, value = 1, step = 0.1),
        sliderInput("mono_sigma", "Smoothing width (higher = smoother transitions)",
                    min = 0.2, max = 4, value = 1, step = 0.1),
        numericInput("mono_floor", "Smoothness floor", value = 0.01, min = 0, step = 0.01),
        actionButton("apply_monotone", "Apply monotone CPT", class = "btn-primary"),
        hr(),
        uiOutput("row_selector"),
        uiOutput("parent_question"),
        hr(),
        uiOutput("sliders_ui"),
        hr(),
        uiOutput("row_sum_indicator"),
        actionButton("save", "Save CPT"),
        actionButton("quit", "Quit"),
        tags$style(HTML("
          tr.bad-row td {
            background-color: #f8d7da !important;
          }
          tr.active-row td {
            background-color: #5485b0 !important;
            color: white !important;
          }
        ")),
        tags$script(HTML("
        (function() {

          function applyHighlight(activeId) {
            var $tbl = $('#cpt_table table');
            if ($tbl.length === 0) return;
            if (!$.fn.dataTable.isDataTable($tbl)) return;

            var table = $tbl.DataTable();

            table.rows().every(function() {
              var d = this.data();
              var rid = parseInt(d[0]); // hidden .row_id
              var tr = this.node();
              if (rid === activeId) $(tr).addClass('active-row');
              else $(tr).removeClass('active-row');
            });
          }

          window.__activeRowId = 1;

          $(document).on('draw.dt', '#cpt_table table', function() {
            applyHighlight(window.__activeRowId);
          });

          Shiny.addCustomMessageHandler('activeRow', function(message) {
            window.__activeRowId = parseInt(message.row_id);
            applyHighlight(window.__activeRowId);
          });

        })();
        "))
      ),
      
      mainPanel(
        fluidRow(
          column(
            width = 12,
            plotOutput("cpt_curve", height = "420px"),
            uiOutput("mono_params_text")
          )
        ),
        fluidRow(
          column(
            width = 12,
            br(),
            plotOutput("prob_bar", height = "170px")
          )
        ),
        hr(),
        DTOutput("cpt_table")
      )
    )
  )
  
  # =========================
  # Server
  # =========================
  
  server <- function(input, output, session) {
    
    cpt_dir <- file.path(pathProCpt, folder)
    if (!dir.exists(cpt_dir)) stop(paste("CPT directory does not exist:", cpt_dir))
    
    params_file <- file.path(cpt_dir, "00_monotone_params_applied.csv")
    
    updating <- reactiveVal(FALSE)
    
    observe({
      files <- list.files(cpt_dir, pattern = "\\.csv$", full.names = FALSE)
      
      if (!is.null(cpt_patterns) && length(cpt_patterns) > 0) {
        rx <- paste(cpt_patterns, collapse = "|")
        files <- files[grepl(rx, files)]
      }
      
      if (length(files) == 0) {
        updateSelectInput(session, "cpt_file", choices = character(0), selected = character(0))
        return()
      }
      
      sel <- isolate(input$cpt_file)
      if (is.null(sel) || !(sel %in% files)) sel <- files[1]
      updateSelectInput(session, "cpt_file", choices = files, selected = sel)
    })
    
    cpt <- reactiveVal(NULL)
    cpt_type <- reactiveVal("with_parents")
    state_cols <- reactiveVal(integer(0))
    
    lock_store <- reactiveVal(NULL)
    current_row <- reactiveVal(1)
    
    slider_obs <- reactiveVal(list())
    lock_obs   <- reactiveVal(NULL)
    
    cpt_proxy <- DT::dataTableProxy("cpt_table")
    
    # --- bad row flags: OFF-GRID IS BAD (your requirement)
    row_bad_flags <- function(df, cpt_type, state_cols, step = PROB_STEP) {
      if (cpt_type != "with_parents") return(rep(FALSE, nrow(df)))
      sc <- state_cols
      K <- length(sc)
      if (K == 0) return(rep(FALSE, nrow(df)))
      
      P <- as.matrix(df[, sc, drop = FALSE])
      storage.mode(P) <- "numeric"
      
      bad_num <- apply(P, 1, function(x) any(!is.finite(x)))
      bad_range <- apply(P, 1, function(x) any(x < step | x > 1))
      
      U <- round(P / step)
      bad_grid <- apply(U, 1, function(u) any(u < 1) || sum(u) != TOTAL_UNITS)
      
      bad_num | bad_range | bad_grid
    }
    
    observe({
      session$sendCustomMessage("activeRow", list(row_id = current_row()))
    })
    
    destroy_observers <- function() {
      old <- slider_obs()
      if (length(old) > 0) {
        lapply(old, function(h) { try(h$destroy(), silent = TRUE); NULL })
      }
      slider_obs(list())
      
      lo <- lock_obs()
      if (!is.null(lo)) {
        try(lo$destroy(), silent = TRUE)
        lock_obs(NULL)
      }
    }
    
    make_slider_observers <- function() {
      req(cpt(), cpt_type(), !is.null(lock_store()))
      
      destroy_observers()
      
      df <- cpt()
      sc <- state_cols()
      n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(sc)
      req(n_controls > 0)
      
      handles <- vector("list", n_controls)
      
      for (i in seq_len(n_controls)) {
        local({
          idx <- i
          
          handles[[idx]] <<- observeEvent(input[[paste0("p_", idx)]], {
            if (updating()) return()
            
            df <- cpt()
            row <- current_row()
            sc  <- state_cols()
            
            n_controls2 <- if (cpt_type() == "no_parents") nrow(df) else length(sc)
            
            lock_vals <- vapply(
              seq_len(n_controls2),
              function(k) isTRUE(input[[paste0("lock_", k)]]),
              logical(1)
            )
            
            ls <- lock_store()
            if (cpt_type() == "no_parents") {
              lock_store(lock_vals)
            } else {
              ls[row, ] <- lock_vals
              lock_store(ls)
            }
            
            old_probs <- if (cpt_type() == "no_parents") {
              as.numeric(df$Probability)
            } else {
              as.numeric(df[row, sc])
            }
            
            # Clamp max so that every OTHER unlocked state can still be >= 0.01
            sum_locked_other <- sum(old_probs[lock_vals & seq_along(old_probs) != idx])
            min_mass_free_others <- sum((!lock_vals) & seq_along(old_probs) != idx) * PROB_STEP
            max_allowed <- max(0, 1 - sum_locked_other - min_mass_free_others)
            max_allowed <- max(max_allowed, PROB_STEP) # idx itself must be >= 0.01
            
            new_value_raw <- as.numeric(input[[paste0("p_", idx)]])
            new_value <- min(max(new_value_raw, PROB_STEP), max_allowed)
            
            if (!isTRUE(all.equal(new_value, new_value_raw))) {
              updating(TRUE)
              on.exit(updating(FALSE), add = TRUE)
              updateSliderInput(session, paste0("p_", idx), value = new_value)
            }
            
            updating(TRUE)
            on.exit(updating(FALSE), add = TRUE)
            
            new_probs <- renormalise(
              probs     = old_probs,
              locked    = lock_vals,
              changed   = idx,
              new_value = new_value
            )
            
            if (cpt_type() == "no_parents") {
              df$Probability <- new_probs
            } else {
              df[row, sc] <- new_probs
            }
            cpt(df)
            
            for (j in seq_along(new_probs)) {
              if (j != idx) updateSliderInput(session, paste0("p_", j), value = new_probs[j])
            }
          }, ignoreInit = TRUE)
        })
      }
      
      slider_obs(handles)
      
      lo <- observeEvent({
        df <- cpt()
        sc <- state_cols()
        n_controls2 <- if (cpt_type() == "no_parents") nrow(df) else length(sc)
        lapply(seq_len(n_controls2), function(i) input[[paste0("lock_", i)]])
      }, {
        req(!is.null(lock_store()))
        df <- cpt()
        sc <- state_cols()
        n_controls2 <- if (cpt_type() == "no_parents") nrow(df) else length(sc)
        row <- current_row()
        
        lock_vals <- vapply(
          seq_len(n_controls2),
          function(k) isTRUE(input[[paste0("lock_", k)]]),
          logical(1)
        )
        
        if (cpt_type() == "no_parents") {
          lock_store(lock_vals)
        } else {
          ls <- lock_store()
          ls[row, ] <- lock_vals
          lock_store(ls)
        }
      }, ignoreInit = TRUE)
      
      lock_obs(lo)
    }
    
    sanitize_for_dt <- function(df) {
      df <- as.data.frame(df, stringsAsFactors = FALSE, check.names = FALSE)
      df[] <- lapply(df, function(col) {
        if (is.factor(col)) col <- as.character(col)
        if (inherits(col, "POSIXlt")) col <- as.character(col)
        if (is.list(col)) col <- vapply(col, function(x) paste(x, collapse = ","), character(1))
        if (is.numeric(col)) col[!is.finite(col)] <- NA_real_
        col
      })
      df
    }
    
    # --- Helpers to compute X axis (parent score) and Y (expected child state)
    get_parent_names <- function(df, sc) {
      st <- infer_cpt_structure(df)
      names(df)[st$parent_cols]
    }
    
    get_child_positions <- function(K) {
      # y-axis meaning: expected child state position (1..K).
      # If you prefer centered scores, replace with .default_scores(...) or similar.
      seq_len(K)
    }
    
    # Read current UI weights into a named numeric vector
    current_parent_weights <- reactive({
      req(cpt(), cpt_type())
      if (cpt_type() != "with_parents") return(NULL)
      
      df <- cpt()
      st <- infer_cpt_structure(df)
      parent_names <- names(df)[st$parent_cols]
      
      w <- setNames(rep(1, length(parent_names)), parent_names)
      for (pn in parent_names) {
        val <- input[[paste0("w_", pn)]]
        if (!is.null(val) && is.finite(val)) w[pn] <- as.numeric(val)
      }
      w
    })
    
    # Read current UI scores into a named list of named numeric vectors
    current_parent_scores <- reactive({
      req(cpt(), cpt_type())
      if (cpt_type() != "with_parents") return(NULL)
      
      df <- cpt()
      st <- infer_cpt_structure(df)
      parent_names <- names(df)[st$parent_cols]
      
      ps <- list()
      for (pn in parent_names) {
        v <- df[[pn]]
        levs <- if (is.factor(v)) levels(v) else unique(as.character(v))
        levs <- as.character(levs)
        
        default_sc <- .default_scores(levs)
        scv <- numeric(length(levs))
        for (i in seq_along(levs)) {
          val <- input[[paste0("score_", pn, "_", i)]]
          scv[i] <- if (is.null(val) || !is.finite(val)) default_sc[i] else as.numeric(val)
        }
        ps[[pn]] <- setNames(scv, levs)
      }
      ps
    })
    
    # Compute total parent score per CPT row using current UI weights/scores
    compute_total_parent_score <- function(df, parent_scores, parent_weights) {
      st <- infer_cpt_structure(df)
      parent_names <- names(df)[st$parent_cols]
      x <- numeric(nrow(df))
      
      for (pn in parent_names) {
        w <- as.numeric(parent_weights[pn])
        mp <- parent_scores[[pn]]
        vals <- as.character(df[[pn]])
        sc <- mp[vals]
        sc[is.na(sc)] <- 0
        x <- x + w * as.numeric(sc)
      }
      x
    }
    
    # Expected child position per row (y = sum p_k * k)
    expected_child_position <- function(P_row) {
      P_row <- as.numeric(P_row)
      P_row[!is.finite(P_row)] <- 0
      s <- sum(P_row)
      if (s <= 0) return(NA_real_)
      P_row <- P_row / s
      pos <- get_child_positions(length(P_row))
      sum(P_row * pos)
    }
    
    # Serialize params for saving
    serialize_weights <- function(w) {
      paste(sprintf("%s=%g", names(w), as.numeric(w)), collapse = ";")
    }
    serialize_scores <- function(ps) {
      # parent{level=score;level=score}|parent{...}
      paste(vapply(names(ps), function(pn) {
        lev <- names(ps[[pn]])
        sc  <- as.numeric(ps[[pn]])
        inner <- paste(sprintf("%s=%g", lev, sc), collapse = ";")
        paste0(pn, "{", inner, "}")
      }, character(1)), collapse = "|")
    }
    
    # ---- Heatmap helpers ----
    
    # Long format: (x, child_index, prob)
    cpt_long_from_matrix <- function(x, Pmat) {
      K <- ncol(Pmat)
      data.frame(
        x = rep(x, each = K),
        child_index = rep(seq_len(K), times = length(x)),
        prob = as.numeric(t(Pmat)),
        stringsAsFactors = FALSE
      )
    }
    
    # Bin x to reduce overplotting (and average probabilities within bins)
    bin_and_average <- function(df_long, n_bins = 60) {
      xr <- range(df_long$x, finite = TRUE)
      if (!all(is.finite(xr)) || diff(xr) <= 0) {
        df_long$x_bin <- df_long$x
      } else {
        bw <- diff(xr) / n_bins
        bw <- max(bw, 1e-6)
        df_long$x_bin <- round(df_long$x / bw) * bw
      }
      
      # mean prob per (x_bin, child_index)
      aggregate(prob ~ x_bin + child_index, data = df_long, FUN = mean)
    }
    
    # Compute monotone preview probabilities on a dense x-grid (long format)
    monotone_preview_long <- function(xg, K, method, temp, floor) {
      out <- vector("list", length(xg))
      
      for (i in seq_along(xg)) {
        s <- xg[i]
        p <- if (method == "ordinal_logit") {
          .probs_ordinal_logit(s, K, temperature = temp)
        } else {
          .probs_gaussian_center(s, K, temperature = temp)
        }
        p <- pmax(p, 0)
        p <- p + floor
        p <- p / sum(p)
        
        out[[i]] <- data.frame(
          x = rep(s, K),
          child_index = seq_len(K),
          prob = p,
          stringsAsFactors = FALSE
        )
      }
      do.call(rbind, out)
    }
    
    # Row score sum to display
    row_score_sum <- reactive({
      req(cpt())
      df <- cpt()
      
      if (cpt_type() != "with_parents") return(rep(NA_real_, nrow(df)))
      
      w  <- current_parent_weights()
      ps <- current_parent_scores()
      req(!is.null(w), !is.null(ps))
      
      compute_total_parent_score(df, parent_scores = ps, parent_weights = w)
    })
    
    # ---- Monotone UI: weights + scores
    output$mono_weights_ui <- renderUI({
      req(cpt(), cpt_type())
      if (cpt_type() != "with_parents") return(NULL)
      st <- infer_cpt_structure(cpt())
      parent_names <- names(cpt())[st$parent_cols]
      tagList(
        tags$h5("Parent weights"),
        lapply(parent_names, function(pn) numericInput(paste0("w_", pn), label = pn, value = 1, step = 0.1))
      )
    })
    
    output$mono_scores_ui <- renderUI({
      req(cpt(), cpt_type())
      if (cpt_type() != "with_parents") return(NULL)
      
      df <- cpt()
      st <- infer_cpt_structure(df)
      parent_names <- names(df)[st$parent_cols]
      
      tagList(
        tags$h5("Parent state scores"),
        lapply(parent_names, function(pn) {
          v <- df[[pn]]
          levs <- if (is.factor(v)) levels(v) else unique(as.character(v))
          levs <- as.character(levs)
          default_sc <- .default_scores(levs)
          
          tags$details(
            tags$summary(paste0("Scores for ", pn)),
            lapply(seq_along(levs), function(i) {
              numericInput(
                inputId = paste0("score_", pn, "_", i),
                label   = levs[i],
                value   = default_sc[i],
                step    = 0.5
              )
            })
          )
        })
      )
    })
    
    observeEvent(input$apply_monotone, {
      req(cpt(), cpt_type())
      if (cpt_type() != "with_parents") {
        showNotification("Monotone generator applies to CPTs with parents (P_* columns).", type = "warning")
        return()
      }
      
      df <- cpt()
      st <- infer_cpt_structure(df)
      parent_names <- names(df)[st$parent_cols]
      
      w <- setNames(rep(1, length(parent_names)), parent_names)
      for (pn in parent_names) {
        val <- input[[paste0("w_", pn)]]
        if (!is.null(val)) w[pn] <- as.numeric(val)
      }
      
      parent_scores <- list()
      for (pn in parent_names) {
        v <- df[[pn]]
        levs <- if (is.factor(v)) levels(v) else unique(as.character(v))
        levs <- as.character(levs)
        
        sc <- numeric(length(levs))
        for (i in seq_along(levs)) {
          val <- input[[paste0("score_", pn, "_", i)]]
          sc[i] <- if (is.null(val)) .default_scores(levs)[i] else as.numeric(val)
        }
        parent_scores[[pn]] <- setNames(sc, levs)
      }
      
      new_df <- build_monotone_cpt(
        df,
        parent_scores  = parent_scores,
        parent_weights = w,
        method         = input$mono_method,
        temperature    = as.numeric(input$mono_center_temp),  # used only for ordinal_logit
        center_temp    = as.numeric(input$mono_center_temp),  # used for gaussian_center
        sigma          = as.numeric(input$mono_sigma),        # used for gaussian_center
        floor          = as.numeric(input$mono_floor)
      )
      
      cpt(new_df)
      state_cols(get_state_columns(new_df))
      lock_store(matrix(FALSE, nrow = nrow(new_df), ncol = length(st$child_cols)))
      showNotification("Monotone CPT applied. You can now tweak row-wise.", type = "message")
      
      
      # ---- Save applied parameters (one line per CPT file)
      params_row <- data.frame(
        cpt_file = input$cpt_file,
        method   = input$mono_method,
        center_temp = as.numeric(input$mono_center_temp),
        sigma       = as.numeric(input$mono_sigma),
        floor    = as.numeric(input$mono_floor),
        weights  = serialize_weights(w),
        scores   = serialize_scores(parent_scores),
        saved_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        stringsAsFactors = FALSE
      )
      
      if (!file.exists(params_file)) {
        write.csv(params_row, params_file, row.names = FALSE)
      } else {
        old <- try(read.csv(params_file, stringsAsFactors = FALSE), silent = TRUE)
        if (inherits(old, "try-error") || is.null(old) || nrow(old) == 0) {
          write.csv(params_row, params_file, row.names = FALSE)
        } else {
          # replace existing row for this CPT if present, else append
          keep <- old$cpt_file != input$cpt_file
          out <- rbind(old[keep, , drop = FALSE], params_row)
          write.csv(out, params_file, row.names = FALSE)
        }
      }
      
      
      destroy_observers()
      make_slider_observers()
    })
    
    # ---- Load CPT ----
    observeEvent(input$cpt_file, {
      req(input$cpt_file)
      
      df <- read.csv(file.path(cpt_dir, input$cpt_file), stringsAsFactors = FALSE)
      
      # IMPORTANT: snap imported CPTs to grid+min so "1/3" becomes e.g. 0.34/0.33/0.33
      ct <- detect_cpt_type(df)
      sc <- get_state_columns(df)
      if (ct == "with_parents" && length(sc) > 0) {
        for (r in seq_len(nrow(df))) {
          df[r, sc] <- snap_probs_to_grid(df[r, sc], step = PROB_STEP, min_prob = PROB_MIN)
        }
      } else if (ct == "no_parents" && "Probability" %in% names(df)) {
        df$Probability <- snap_probs_to_grid(df$Probability, step = PROB_STEP, min_prob = PROB_MIN)
      }
      
      cpt(df)
      cpt_type(ct)
      state_cols(sc)
      
      if (cpt_type() == "no_parents") {
        lock_store(rep(FALSE, nrow(df)))
      } else {
        n_states <- length(state_cols())
        lock_store(matrix(FALSE, nrow = nrow(df), ncol = n_states))
      }
      
      destroy_observers()
      current_row(1)
      make_slider_observers()
    })
    
    # ---- Row selector ----
    output$row_selector <- renderUI({
      req(cpt(), current_row())
      if (cpt_type() == "with_parents") {
        selectInput("row_id", "Parent configuration (line)",
                    choices = seq_len(nrow(cpt())),
                    selected = current_row())
      }
    })
    
    observeEvent(input$row_id, {
      current_row(as.integer(input$row_id))
    })
    
    observeEvent(list(current_row(), cpt(), state_cols(), lock_store()), {
      req(cpt(), !is.null(lock_store()))
      
      df <- cpt()
      row <- current_row()
      sc  <- state_cols()
      
      n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(sc)
      req(n_controls > 0)
      
      probs <- if (cpt_type() == "no_parents") as.numeric(df$Probability) else as.numeric(df[row, sc])
      locks <- if (cpt_type() == "no_parents") lock_store() else lock_store()[row, ]
      
      updating(TRUE)
      on.exit(updating(FALSE), add = TRUE)
      
      for (i in seq_along(probs)) {
        updateSliderInput(session, paste0("p_", i), value = probs[i])
        updateCheckboxInput(session, paste0("lock_", i), value = isTRUE(locks[i]))
      }
    }, ignoreInit = TRUE)
    
    observeEvent(input$table_dblclick_row, {
      req(cpt_type())
      if (cpt_type() != "with_parents") return()
      i <- as.integer(input$table_dblclick_row)
      if (is.na(i)) return()
      current_row(i)
      updateSelectInput(session, "row_id", selected = i)
    })
    
    # Parent question (unchanged)
    output$parent_question <- renderUI({
      req(cpt(), cpt_type(), input$cpt_file)
      
      df <- cpt()
      clean <- function(x) trimws(gsub("_+", " ", as.character(x)))
      
      child_node <- gsub("\\.csv$", "", input$cpt_file)
      child_node <- gsub("^CPT_+", "", child_node)
      child_node <- clean(child_node)
      if (nchar(child_node) == 0) child_node <- "the node"
      
      if (cpt_type() == "no_parents") {
        return(tags$div(
          style = "margin-top: 6px; margin-bottom: 6px;",
          tags$span("What is the probability that "),
          tags$span(style = "font-weight: bold;", paste0("'", child_node, "'")),
          tags$span(" is:")
        ))
      }
      
      sc <- state_cols()
      row <- current_row()
      parent_cols <- setdiff(seq_along(df), sc)
      req(length(parent_cols) > 0)
      
      parent_names <- clean(names(df)[parent_cols])
      parent_names <- gsub(".", "-", parent_names, fixed = TRUE)
      parent_vals  <- clean(df[row, parent_cols, drop = TRUE])
      
      pal <- c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00")
      
      parent_color_maps <- lapply(parent_cols, function(j) {
        vals <- df[[j]]
        u_ord <- if (is.factor(vals)) levels(vals) else unique(clean(vals))
        u_ord <- clean(u_ord)
        cols <- pal[round(seq(1, length(pal), length.out = length(u_ord)))]
        stats::setNames(cols, u_ord)
      })
      
      cond_fragments <- lapply(seq_along(parent_cols), function(k) {
        p_name <- parent_names[k]
        p_val  <- parent_vals[k]
        cmap   <- parent_color_maps[[k]]
        colhex <- unname(cmap[p_val])
        
        tags$span(
          style = "font-weight: bold;",
          tags$span(paste0(p_name, " is ")),
          tags$span(style = paste0("color:", colhex, ";"), paste0(p_val))
        )
      })
      
      cond_with_and <- list()
      for (k in seq_along(cond_fragments)) {
        cond_with_and <- c(cond_with_and, list(cond_fragments[[k]]))
        if (k < length(cond_fragments)) cond_with_and <- c(cond_with_and, list(tags$span(" and ")))
      }
      
      tags$div(
        style = "margin-top: 6px; margin-bottom: 6px;",
        tags$span("When "), cond_with_and,
        tags$span(", what is the probability that "),
        tags$span(style = "font-style: italic;", paste0("'", child_node, "'")),
        tags$span(" is:")
      )
    })
    
    
    output$mono_params_text <- renderUI({
      req(input$cpt_file)
      
      if (!file.exists(params_file)) {
        return(tags$div(
          style = "font-size: 10px; margin-top: 6px; padding: 8px; border: 1px solid #ddd; border-radius: 6px; max-height: 140px; overflow-y: auto;",
          "No saved monotone parameters for any CPT yet."
        ))
      }
      
      tab <- try(read.csv(params_file, stringsAsFactors = FALSE), silent = TRUE)
      if (inherits(tab, "try-error") || is.null(tab) || nrow(tab) == 0) {
        return(tags$div(
          style = "font-size: 12px; color: #666; margin-top: 6px;",
          "Parameter registry is empty or unreadable."
        ))
      }
      
      row <- tab[tab$cpt_file == input$cpt_file, , drop = FALSE]
      if (nrow(row) == 0) {
        return(tags$div(
          style = "font-size: 12px; color: #666; margin-top: 6px;",
          "No saved monotone parameters for this CPT."
        ))
      }
      
      r <- row[nrow(row), , drop = FALSE]
      
      tags$div(
        style = "font-size: 12px; margin-top: 6px; padding: 8px; border: 1px solid #ddd; border-radius: 6px;",
        tags$div(tags$b("Last applied monotone parameters for this CPT")),
        tags$div(sprintf("Saved at: %s", r$saved_at)),
        tags$div(sprintf("Method: %s | CenterTemp: %s | Sigma: %s | Floor: %s",
                         r$method, r$center_temp, r$sigma, r$floor)),
        tags$div(style="margin-top:4px;", tags$b("Weights:"), tags$div(style="color:#444;", r$weights)),
        tags$div(style="margin-top:4px;", tags$b("Scores:"), tags$div(style="color:#444; word-break: break-word;", r$scores))
      )
    })
    
    
    # ---- Sliders ----
    output$sliders_ui <- renderUI({
      req(cpt(), cpt_type(), !is.null(lock_store()))
      df <- cpt()
      
      if (cpt_type() == "no_parents") {
        probs <- as.numeric(df$Probability)
        locks <- as.logical(lock_store())
        state_names <- as.character(df$State)
        
        lapply(seq_along(probs), function(i) {
          tags$div(
            style = "margin-bottom: 10px; padding-bottom: 6px; border-bottom: 1px solid #eee;",
            tags$div(style = "font-weight: bold; margin-bottom: 4px; text-align: center;",
                     gsub("_"," ",state_names[i])),
            checkboxInput(paste0("lock_", i), label = "Lock", value = isTRUE(locks[i])),
            sliderInput(paste0("p_", i), label = NULL, min = PROB_MIN, max = 1, step = PROB_STEP, value = probs[i])
          )
        })
      } else {
        n_controls <- length(state_cols())
        req(n_controls > 0)
        
        row <- current_row()
        sc <- state_cols()
        
        probs <- as.numeric(df[row, sc])
        locks <- as.logical(lock_store()[row, ])
        state_names <- sub("^P_", "", names(df)[sc])
        
        lapply(seq_along(probs), function(i) {
          tags$div(
            style = "margin-bottom: 10px; padding-bottom: 6px; border-bottom: 1px solid #eee;",
            tags$div(style = "font-weight: bold; margin-bottom: 4px; text-align: center;",
                     gsub("_"," ",state_names[i])),
            checkboxInput(paste0("lock_", i), label = "Lock", value = isTRUE(locks[i])),
            sliderInput(paste0("p_", i), label = NULL, min = PROB_MIN, max = 1, step = PROB_STEP, value = probs[i])
          )
        })
      }
    })
    
    output$row_sum_indicator <- renderUI({
      req(cpt())
      
      df <- cpt()
      sc <- state_cols()
      
      probs <- if (cpt_type() == "no_parents") as.numeric(df$Probability) else as.numeric(df[current_row(), sc])
      U <- round(probs / PROB_STEP)
      
      ok <- (sum(U) == TOTAL_UNITS) && all(U >= 1)
      s <- sum(probs)
      
      col <- if (ok) "#2E7D32" else "#C62828"
      txt <- if (ok) sprintf("Row sum: %.2f ✓ (grid OK)", s) else sprintf("Row sum: %.2f ✗ (must be grid-sum 1.00 and each ≥ 0.01)", s)
      
      tags$div(
        style = paste0("padding:8px 10px; border-radius:6px; font-weight:600; color:white; ",
                       "background:", col, "; margin-top:6px; margin-bottom:6px;"),
        txt
      )
    })
    
    
    output$cpt_curve <- renderPlot({
      req(cpt(), cpt_type())
      if (cpt_type() != "with_parents") {
        ggplot() + theme_void() + annotate("text", x = 0, y = 0, label = "Heatmap view is for CPTs with parents (P_* columns).")
      } else {
        
        df <- cpt()
        sc <- state_cols()
        req(length(sc) > 0)
        
        # Current axis definition (must match what users tune)
        w  <- current_parent_weights()
        ps <- current_parent_scores()
        req(!is.null(w), !is.null(ps))
        
        # Parent-score per row (x)
        x_rows <- compute_total_parent_score(df, parent_scores = ps, parent_weights = w)
        
        # CPT probabilities wide -> long
        Pmat <- as.matrix(df[, sc, drop = FALSE])
        storage.mode(Pmat) <- "numeric"
        child_states <- sub("^P_", "", names(df)[sc])
        K <- length(child_states)
        
        # Build long table: one rectangle per (row, child_state)
        long <- as.data.frame(Pmat)
        names(long) <- child_states
        long$.row <- seq_len(nrow(long))
        long$score <- x_rows
        
        long <- long |>
          tidyr::pivot_longer(cols = all_of(child_states),
                              names_to = "child_state",
                              values_to = "prob") |>
          mutate(
            child_state = factor(child_state, levels = child_states),
            y = as.integer(child_state)
          ) |>
          filter(is.finite(score), is.finite(prob))
        
        validate(need(nrow(long) > 0, "No valid CPT rows to plot."))
        
        # ---- IMPORTANT: build per-row xmin/xmax (no bins)
        # Use unique score per row; if duplicates exist, we will slightly jitter them for display
        # (otherwise rectangles collapse onto the same x)
        row_scores <- long |>
          distinct(.row, score) |>
          arrange(score) |>
          mutate(
            score_plot = score
          )
        
        # handle duplicated scores: tiny deterministic jitter by row index
        dup <- duplicated(row_scores$score_plot) | duplicated(row_scores$score_plot, fromLast = TRUE)
        if (any(dup)) {
          # jitter magnitude is tiny relative to axis range
          rng <- diff(range(row_scores$score_plot))
          eps <- if (is.finite(rng) && rng > 0) rng * 1e-4 else 1e-4
          row_scores <- row_scores |>
            group_by(score_plot) |>
            mutate(score_plot = score_plot + (row_number() - 1) * eps) |>
            ungroup()
        }
        
        # Compute boundaries by midpoints between neighbors
        s <- row_scores$score_plot
        n <- length(s)
        if (n == 1) {
          # arbitrary width if only one row
          dx <- 1
          xmin <- s - dx/2
          xmax <- s + dx/2
        } else {
          mids <- (s[-1] + s[-n]) / 2
          xmin <- c(s[1] - (mids[1] - s[1]), mids)
          xmax <- c(mids, s[n] + (s[n] - mids[n-1]))
        }
        
        row_scores$xmin <- xmin
        row_scores$xmax <- xmax
        
        # join xmin/xmax back to long
        long <- long |>
          left_join(row_scores |> select(.row, xmin, xmax, score_plot), by = ".row")
        
        # ---- OPTIONAL: overlay contour of the current parameterized monotone curve surface
        # Build dense x-grid and compute p(state|x) using your method
        method     <- input$mono_method %||% "gaussian_center"
        center_tmp <- input$mono_center_temp %||% 1
        sigma_val  <- input$mono_sigma %||% 1
        floor      <- input$mono_floor %||% 0.01
        
        # grid for preview (smooth)
        xg <- seq(min(row_scores$score_plot), max(row_scores$score_plot), length.out = 220)
        
        Pfit <- vapply(xg, function(xval) {
          p <- if (method == "ordinal_logit") {
            .probs_ordinal_logit(xval, K, temperature = center_tmp)
          } else {
            .probs_gaussian_center_cs(xval, K, center_temp = center_tmp, sigma = sigma_val)
          }
          p <- pmax(p, 0)
          p <- p + floor
          p <- p / sum(p)
          p
        }, numeric(K))
        Pfit <- t(Pfit)  # length(xg) x K
        
        fit_long <- data.frame(
          x = rep(xg, times = K),
          y = rep(seq_len(K), each = length(xg)),
          prob = as.vector(Pfit)
        )
        
        # ---- Plot
        ggplot() +
          geom_rect(
            data = long,
            aes(xmin = xmin, xmax = xmax, ymin = y - 0.5, ymax = y + 0.5, fill = prob),
            color = NA
          ) +
          # contour overlay from parameterized curve (thin)
          { if (isTRUE(input$show_contours)) {
            geom_contour(
              data = fit_long,
              aes(x = x, y = y, z = prob),
              bins = 7,
              linewidth = 0.3,
              alpha = 0.6
            )
          }
          } +
          scale_fill_gradient(low = "white", high = "red", limits = c(0, 1), name = "P") +
          scale_y_continuous(
            breaks = seq_len(K),
            labels = child_states,
            expand = expansion(add = 0.02)
          ) +
          labs(
            x = "Parent score (current weights + scores)",
            y = "Child state",
            title = "CPT probability surface (one tile per CPT row) + monotone preview contours"
          ) +
          theme_minimal(base_size = 10) +
          theme(
            legend.position = "right",
            plot.title = element_text(size = 11),
            axis.text.y = element_text(size = 8)
          )
      }
    }, res = 120)
    
    
    # ---- bar plot ----
    output$prob_bar <- renderPlot({
      req(cpt())
      df <- cpt()
      
      if (is.null(df)) return()
      
      if (cpt_type() == "no_parents") {
        req("Probability" %in% names(df), "State" %in% names(df))
        probs <- as.numeric(df$Probability)
        state_names <- as.character(df$State)
      } else {
        sc <- state_cols()
        req(length(sc) > 0)
        probs <- as.numeric(df[current_row(), sc])
        state_names <- sub("^P_", "", names(df)[sc])
      }
      
      # Defensive cleanup (should already be grid-valid, but keep safe for plotting)
      probs[!is.finite(probs)] <- 0
      probs <- pmax(probs, 0)
      s <- sum(probs)
      if (s <= 0) return()
      probs <- probs / s
      
      # Colors (same palette as before)
      ramp5 <- c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00")
      n_states <- length(probs)
      cols <- ramp5[round(seq(1, 5, length.out = n_states))]
      
      lefts  <- c(0, cumsum(probs)[-length(probs)])
      rights <- cumsum(probs)
      
      op <- par(mar = c(2.2, 1.0, 0.5, 0.5), xaxs = "i", yaxs = "i")
      on.exit(par(op), add = TRUE)
      
      plot(NA,
           xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Probability mass", ylab = "",
           yaxt = "n", xaxt = "n")
      
      min_width_label <- 0.05
      cex_max <- 0.85
      cex_min <- 0.55
      
      for (i in seq_along(probs)) {
        rect(lefts[i], 0.25, rights[i], 0.75, col = cols[i], border = "black")
        
        width <- rights[i] - lefts[i]
        if (width >= min_width_label) {
          scaled <- cex_min + (pmin(width, 0.20) - 0.05) * (cex_max - cex_min) / (0.20 - 0.05)
          cex <- max(cex_min, min(cex_max, scaled))
          mid <- (lefts[i] + rights[i]) / 2
          text(mid, 0.50,
               labels = sprintf("%s\n%.2f", gsub("_", " ", state_names[i]), probs[i]),
               cex = cex)
        }
      }
      
      ticks_minor <- seq(0, 1, by = PROB_STEP)
      ticks_major <- seq(0, 1, by = 0.1)
      axis(1, at = ticks_minor, labels = FALSE, tck = -0.03)
      axis(1, at = ticks_major, labels = sprintf("%.1f", ticks_major), cex.axis = 0.8)
    }, res = 120)
    
    # ---- Table ----
    output$cpt_table <- renderDT({
      req(cpt())
      
      df <- cpt()
      sc <- state_cols()
      
      prob_names <- character(0)
      if (cpt_type() == "with_parents" && length(sc) > 0) prob_names <- names(df)[sc]
      if (cpt_type() == "no_parents" && "Probability" %in% names(df)) prob_names <- "Probability"
      
      df_show <- df
      # ---- add computed parent score sum column (between parents and P_* columns)
      if (cpt_type() == "with_parents") {
        sc_idx <- state_cols()
        if (length(sc_idx) > 0) {
          insert_pos <- min(sc_idx)  # first P_ column position in df
          score_vec <- row_score_sum()
          
          # Build new df_show with insertion before first P_ col
          df_show <- cbind(
            df_show[, seq_len(insert_pos - 1), drop = FALSE],
            ScoreSum = round(score_vec, 3),
            df_show[, insert_pos:ncol(df_show), drop = FALSE]
          )
        }
      }
      if (length(prob_names) > 0) df_show[, prob_names] <- round(df_show[, prob_names], 2)
      
      df_show <- sanitize_for_dt(df_show)
      names(df_show) <- make.unique(names(df_show))
      
      bad <- row_bad_flags(df = df, cpt_type = cpt_type(), state_cols = state_cols())
      
      df_show <- cbind(.row_id = seq_len(nrow(df_show)),
                       .row_bad = as.integer(bad),
                       df_show)
      
      datatable(
        df_show,
        selection = "none",
        rownames = FALSE,
        options = list(
          pageLength = 100,
          columnDefs = list(
            list(targets = c(0, 1), visible = FALSE, searchable = FALSE)
          ),
          rowCallback = JS("
            function(row, data, index) {

              if (parseInt(data[1]) === 1) $(row).addClass('bad-row');
              else $(row).removeClass('bad-row');

              $(row).off('dblclick');
              $(row).on('dblclick', function() {
                Shiny.setInputValue('table_dblclick_row', parseInt(data[0]), {priority: 'event'});
              });
            }
          ")
        )
      )
    }, server = TRUE)
    
    observeEvent(list(
      cpt(),
      input$cpt_table_rows_current,
      current_parent_weights(),
      current_parent_scores()
    ), {
      
      req(cpt(), !is.null(input$cpt_table_rows_current))
      
      df <- cpt()
      sc <- state_cols()
      
      prob_names <- character(0)
      if (cpt_type() == "with_parents" && length(sc) > 0) prob_names <- names(df)[sc]
      if (cpt_type() == "no_parents" && "Probability" %in% names(df)) prob_names <- "Probability"
      
      df_show <- df
      # ---- add computed parent score sum column (between parents and P_* columns)
      if (cpt_type() == "with_parents") {
        sc_idx <- state_cols()
        if (length(sc_idx) > 0) {
          insert_pos <- min(sc_idx)
          score_vec <- row_score_sum()
          
          df_show <- cbind(
            df_show[, seq_len(insert_pos - 1), drop = FALSE],
            ScoreSum = round(score_vec, 3),
            df_show[, insert_pos:ncol(df_show), drop = FALSE]
          )
        }
      }
      if (length(prob_names) > 0) df_show[, prob_names] <- round(df_show[, prob_names], 2)
      
      df_show <- sanitize_for_dt(df_show)
      names(df_show) <- make.unique(names(df_show))
      
      bad <- row_bad_flags(df = df, cpt_type = cpt_type(), state_cols = state_cols())
      
      df_show <- cbind(.row_id = seq_len(nrow(df_show)),
                       .row_bad = as.integer(bad),
                       df_show)
      
      DT::replaceData(cpt_proxy, df_show, resetPaging = FALSE, rownames = FALSE)
    }, ignoreInit = TRUE)
    
    # ---- Save ----
    observeEvent(input$save, {
      req(cpt())
      df <- cpt()
      sc <- state_cols()
      
      if (length(sc) > 0) df[, sc] <- round(df[, sc], 2)
      write.csv(df, file.path(pathProCpt, folder, input$cpt_file), row.names = FALSE)
      
      cpt(df)
      showNotification("CPT saved", type = "message")
    })
    
    observeEvent(input$quit, {
      stopApp()
    })
  }
  
  shinyApp(ui, server)
}

if (interactive()) {
  run_cpt_app(pathProCpt = pathProCpt, folder = "01_CPT")
}