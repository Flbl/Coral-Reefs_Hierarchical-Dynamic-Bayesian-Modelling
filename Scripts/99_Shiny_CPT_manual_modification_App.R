#############################################################################################
#                                                                                           #
#                                                                                           #
#####             Simple Shiny interface to modify CPTs with mouse cursors              #####
#                                                                                           #
#                                                                                           #
#############################################################################################

library(shiny)
library(DT)

# =========================
# Utility functions
# =========================

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

.probs_gaussian_center <- function(score, K, temperature = 1) {
  temperature <- max(1e-6, temperature)
  center <- 1 + (K - 1) * plogis(score / temperature)
  sigma <- max(1e-3, temperature)
  idx <- seq_len(K)
  w <- exp(-0.5 * ((idx - center) / sigma)^2)
  w / sum(w)
}

build_monotone_cpt <- function(df,
                               parent_scores = NULL,
                               parent_weights = NULL,
                               method = c("ordinal_logit", "gaussian_center"),
                               temperature = 1,
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
    p <- if (method == "ordinal_logit") .probs_ordinal_logit(s, K, temperature = temperature) else .probs_gaussian_center(s, K, temperature = temperature)
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
        uiOutput("row_selector"),
        uiOutput("parent_question"),
        hr(),
        uiOutput("sliders_ui"),
        hr(),
        uiOutput("row_sum_indicator"),
        tags$h4("Monotone CPT generator"),
        uiOutput("mono_weights_ui"),
        uiOutput("mono_scores_ui"),
        selectInput("mono_method", "Method",
                    choices = c("Ordinal logit" = "ordinal_logit",
                                "Gaussian center" = "gaussian_center"),
                    selected = "gaussian_center"),
        sliderInput("mono_temp", "Temperature / sharpness",
                    min = 0.1, max = 5, value = 1, step = 0.1),
        numericInput("mono_floor", "Smoothness floor", value = 0.01, min = 0, step = 0.01),
        actionButton("apply_monotone", "Apply monotone CPT", class = "btn-primary"),
        hr(),
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
        plotOutput("prob_bar", height = "160px"),
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
        parent_scores = parent_scores,
        parent_weights = w,
        method = input$mono_method,
        temperature = input$mono_temp,
        floor = input$mono_floor
      )
      
      cpt(new_df)
      state_cols(get_state_columns(new_df))
      lock_store(matrix(FALSE, nrow = nrow(new_df), ncol = length(st$child_cols)))
      showNotification("Monotone CPT applied. You can now tweak row-wise.", type = "message")
      
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
    
    observeEvent(list(cpt(), input$cpt_table_rows_current), {
      req(cpt(), !is.null(input$cpt_table_rows_current))
      
      df <- cpt()
      sc <- state_cols()
      
      prob_names <- character(0)
      if (cpt_type() == "with_parents" && length(sc) > 0) prob_names <- names(df)[sc]
      if (cpt_type() == "no_parents" && "Probability" %in% names(df)) prob_names <- "Probability"
      
      df_show <- df
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