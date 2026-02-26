#############################################################################################
#                                                                                           #
#####             Simple Shiny interface to modify CPTs with mouse cursors              #####
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
  if (all(c("State", "Probability") %in% names(df))) "no_parents" else "with_parents"
}

get_state_columns <- function(df) {
  res <- grep("^P_", names(df))
  if (length(res) == 0 && "Probability" %in% names(df)) res <- which(names(df) == "Probability")
  res
}

PROB_STEP   <- 0.01
PROB_MIN    <- PROB_STEP
TOTAL_UNITS <- round(1 / PROB_STEP)

snap_probs_to_grid <- function(probs, step = PROB_STEP, min_prob = PROB_MIN) {
  probs <- as.numeric(probs)
  K <- length(probs)
  if (K == 0) return(probs)
  
  probs[!is.finite(probs)] <- min_prob
  probs <- pmax(probs, min_prob)
  probs <- pmin(probs, 1)
  
  U <- round(probs / step)
  U[U < 1] <- 1L
  
  while (sum(U) > TOTAL_UNITS) {
    j <- which.max(U)
    if (U[j] <= 1) break
    U[j] <- U[j] - 1L
  }
  while (sum(U) < TOTAL_UNITS) {
    j <- which.max(U)
    U[j] <- U[j] + 1L
  }
  U * step
}

renormalise <- function(probs, locked, changed, new_value, step = PROB_STEP) {
  n <- length(probs)
  locked <- as.logical(locked)
  total_units <- round(1 / step)
  
  probs <- as.numeric(probs)
  probs[!is.finite(probs)] <- step
  
  units <- round(probs / step)
  
  idx_all <- seq_len(n)
  locked_others <- locked & idx_all != changed
  free_idx <- which(!locked & idx_all != changed)
  
  min_units <- integer(n)
  min_units[!locked] <- 1L
  
  used_units_locked <- sum(units[locked_others])
  min_required_except_changed <- sum(min_units[idx_all != changed])
  
  max_units_changed <- total_units - used_units_locked - min_required_except_changed
  max_units_changed <- max(max_units_changed, min_units[changed])
  
  target_units_changed <- round(new_value / step)
  target_units_changed <- max(target_units_changed, min_units[changed])
  target_units_changed <- min(target_units_changed, max_units_changed)
  
  units[changed] <- target_units_changed
  units[free_idx] <- min_units[free_idx]
  
  remaining_units <- total_units - used_units_locked - units[changed] - sum(units[free_idx])
  
  if (length(free_idx) > 0 && remaining_units > 0) {
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
  
  diff <- total_units - sum(units)
  if (diff != 0) {
    if (length(free_idx) > 0) {
      j <- free_idx[1]
      units[j] <- units[j] + diff
      if (units[j] < min_units[j]) {
        d2 <- min_units[j] - units[j]
        units[j] <- min_units[j]
        units[changed] <- units[changed] - d2
      }
    } else {
      units[changed] <- units[changed] + diff
    }
  }
  
  units[!locked] <- pmax(units[!locked], 1L)
  diff2 <- total_units - sum(units)
  if (diff2 != 0) {
    j <- which(!locked)[1]
    units[j] <- units[j] + diff2
  }
  
  units * step
}

# =========================
# Monotone CPT generation
# =========================

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

.probs_gaussian_center_cs <- function(score, K, center_temp = 1, sigma = 1) {
  center_temp <- max(1e-6, center_temp)
  sigma <- max(1e-6, sigma)
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
    p <- snap_probs_to_grid(p, step = PROB_STEP, min_prob = PROB_MIN)
    P[i, ] <- p
  }
  
  out <- df
  out[, child_cols] <- P
  out
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

# =========================
# App
# =========================

run_cpt_app <- function(pathProCpt,
                        folder = "01_CPT",
                        cpt_patterns = NULL) {
  
  ui <- fluidPage(
    titlePanel("Bayesian Network CPT Editor"),
    
    sidebarLayout(
      sidebarPanel(
        selectInput("cpt_file", "Select CPT", choices = NULL),
        
        checkboxInput("mono_on", "Monotone CPT generator", value = FALSE),
        
        conditionalPanel(
          condition = "input.mono_on == true",
          tags$div(
            style = "margin-top: 6px; padding: 8px; border: 1px solid #ddd; border-radius: 6px; background: #fafafa;",
            checkboxInput("show_contours", "Show contours", value = TRUE),
            uiOutput("mono_weights_ui"),
            uiOutput("mono_scores_ui"),
            selectInput("mono_method", "Method",
                        choices = c("Ordinal logit" = "ordinal_logit",
                                    "Gaussian center" = "gaussian_center"),
                        selected = "gaussian_center"
            ),
            sliderInput("mono_center_temp", "Center saturation (lower = stronger extremes)",
                        min = 0.1, max = 5, value = 1, step = 0.1
            ),
            sliderInput("mono_sigma", "Smoothing width (higher = smoother transitions)",
                        min = 0.2, max = 4, value = 1, step = 0.1
            ),
            numericInput("mono_floor", "Smoothness floor", value = 0.01, min = 0, step = 0.01),
            uiOutput("mono_status"),
            tags$div(
              style = "margin-top: 6px; font-size: 12px; color: #8a1f1f;",
              tags$b("Warning: "),
              "Changing generator parameters recomputes the CPT from the loaded baseline and ",
              tags$u("resets any manual edits."),
              tags$br(),
              "If you start editing row-wise sliders, don't tweak the generator parameters any further !"
            )
          ),
          hr()
        ),
        
        uiOutput("row_selector"),
        uiOutput("parent_question"),
        hr(),
        uiOutput("sliders_ui"),
        hr(),
        uiOutput("row_sum_indicator"),
        actionButton("save", "Save CPT"),
        actionButton("quit", "Quit"),
        
        tags$style(HTML("
          tr.bad-row td { background-color: #f8d7da !important; }
          tr.active-row td { background-color: #5485b0 !important; color: white !important; }
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
        conditionalPanel(
          condition = "input.mono_on == true",
          fluidRow(
            column(
              width = 12,
              plotOutput("cpt_curve", height = "420px"),
              uiOutput("mono_params_text")
            )
          ),
          hr()
        ),
        
        fluidRow(
          column(width = 12, br(), plotOutput("prob_bar", height = "170px"))
        ),
        hr(),
        DTOutput("cpt_table")
      )
    )
  )
  
  server <- function(input, output, session) {
    
    cpt_dir <- file.path(pathProCpt, folder)
    if (!dir.exists(cpt_dir)) stop(paste("CPT directory does not exist:", cpt_dir))
    
    updating   <- reactiveVal(FALSE)
    mono_busy  <- reactiveVal(FALSE)
    
    cpt       <- reactiveVal(NULL)
    cpt_base  <- reactiveVal(NULL)
    cpt_type  <- reactiveVal("with_parents")
    state_cols <- reactiveVal(integer(0))
    
    lock_store <- reactiveVal(NULL)
    current_row <- reactiveVal(1)
    
    slider_obs <- reactiveVal(list())
    lock_obs   <- reactiveVal(NULL)
    
    # Parent metadata
    parent_meta <- reactiveVal(NULL)   # list(parent_names, levels_by_parent)
    
    cpt_proxy <- DT::dataTableProxy("cpt_table")
    
    # -------------------------
    # Helper holding updating until flush
    # -------------------------
    hold_updating_until_flush <- function(expr) {
      updating(TRUE)
      force(expr)
      
      session$onFlushed(function() {
        # Critical: frozen inputs "thaw" and trigger observers right after flush.
        # Keep updating==TRUE for one more event-loop tick so slider observers ignore them.
        later::later(function() updating(FALSE), 0)
      }, once = TRUE)
    }
    
    safe_update_slider <- function(id, value) {
      freezeReactiveValue(input, id)
      updateSliderInput(session, id, value = value)
    }
    
    safe_update_checkbox <- function(id, value) {
      freezeReactiveValue(input, id)
      updateCheckboxInput(session, id, value = value)
    }
    
    # -------------------------
    # File list
    # -------------------------
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
    
    # -------------------------
    # Row highlight message
    # -------------------------
    observeEvent(current_row(), {
      session$sendCustomMessage("activeRow", list(row_id = current_row()))
    }, ignoreInit = FALSE)
    
    # -------------------------
    # Bad row flags (grid validity)
    # -------------------------
    row_bad_flags <- function(df, ct, sc, step = PROB_STEP) {
      if (ct != "with_parents") return(rep(FALSE, nrow(df)))
      if (length(sc) == 0) return(rep(FALSE, nrow(df)))
      
      P <- as.matrix(df[, sc, drop = FALSE])
      storage.mode(P) <- "numeric"
      
      bad_num   <- apply(P, 1, function(x) any(!is.finite(x)))
      bad_range <- apply(P, 1, function(x) any(x < step | x > 1))
      U <- round(P / step)
      bad_grid  <- apply(U, 1, function(u) any(u < 1) || sum(u) != TOTAL_UNITS)
      
      bad_num | bad_range | bad_grid
    }
    
    # -------------------------
    # Slider observer lifecycle
    # -------------------------
    destroy_observers <- function() {
      old <- slider_obs()
      if (length(old) > 0) lapply(old, function(h) try(h$destroy(), silent = TRUE))
      slider_obs(list())
      
      lo <- lock_obs()
      if (!is.null(lo)) try(lo$destroy(), silent = TRUE)
      lock_obs(NULL)
    }
    
    # -------------------------
    # Slider observers
    # -------------------------
    make_slider_observers <- function() {
      req(cpt(), cpt_type(), !is.null(lock_store()))
      
      destroy_observers()
      
      df0 <- cpt()
      sc0 <- state_cols()
      n_controls <- if (cpt_type() == "no_parents") nrow(df0) else length(sc0)
      req(n_controls > 0)
      
      handles <- vector("list", n_controls)
      
      for (i in seq_len(n_controls)) {
        local({
          idx <- i
          pid <- paste0("p_", idx)
          
          handles[[idx]] <<- observeEvent(input[[pid]], {
            # If we are in a programmatic UI update, do nothing (prevents feedback loops)
            if (isTRUE(updating())) return()
            
            req(cpt(), !is.null(lock_store()))
            df  <- cpt()
            row <- isolate(current_row())
            sc  <- isolate(state_cols())
            n_controls2 <- if (cpt_type() == "no_parents") nrow(df) else length(sc)
            
            # Collect lock state (use isolate so locks don't retrigger mid-handler)
            lock_vals <- vapply(
              seq_len(n_controls2),
              function(k) isTRUE(isolate(input[[paste0("lock_", k)]])),
              logical(1)
            )
            
            # Persist lock store for this row (or global for no-parents)
            if (cpt_type() == "no_parents") {
              lock_store(lock_vals)
            } else {
              ls <- lock_store()
              ls[row, ] <- lock_vals
              lock_store(ls)
            }
            
            old_probs <- if (cpt_type() == "no_parents") {
              as.numeric(df$Probability)
            } else {
              as.numeric(df[row, sc])
            }
            
            # compute bounded new_value
            sum_locked_other <- sum(old_probs[lock_vals & seq_along(old_probs) != idx])
            min_mass_free_others <- sum((!lock_vals) & seq_along(old_probs) != idx) * PROB_STEP
            max_allowed <- max(0, 1 - sum_locked_other - min_mass_free_others)
            max_allowed <- max(max_allowed, PROB_STEP)
            
            new_value_raw <- as.numeric(input[[pid]])
            new_value <- min(max(new_value_raw, PROB_STEP), max_allowed)
            
            # Everything inside a "programmatic update" guard.
            # This prevents the sync observer + slider observers from chasing each other.
            hold_updating_until_flush({
              # optional clamp: OK to push *this one* back
              if (!isTRUE(all.equal(new_value, new_value_raw))) {
                safe_update_slider(pid, new_value)
              }
              
              new_probs <- renormalise(
                probs     = old_probs,
                locked    = lock_vals,
                changed   = idx,
                new_value = new_value
              )
              
              # Update CPT only. Do NOT push other sliders here.
              # The sync observer will push the full consistent vector once.
              if (cpt_type() == "no_parents") df$Probability <- new_probs else df[row, sc] <- new_probs
              cpt(df)
            })
          }, ignoreInit = TRUE)
        })
      }
      
      slider_obs(handles)
      
      # Lock observer (store lock toggles, but do NOT touch sliders here)
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
        row <- isolate(current_row())
        
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
    
    # -------------------------
    # Monotone parameters from UI
    # -------------------------
    current_parent_weights <- reactive({
      pm <- parent_meta()
      if (is.null(pm)) return(NULL)
      
      parent_names <- pm$parent_names
      w <- setNames(rep(1, length(parent_names)), parent_names)
      
      for (pn in parent_names) {
        id <- paste0("w_", pn)
        val <- input[[id]]
        if (!is.null(val) && is.finite(val)) w[pn] <- as.numeric(val)
      }
      w
    })
    
    current_parent_scores <- reactive({
      pm <- parent_meta()
      if (is.null(pm)) return(NULL)
      
      parent_names <- pm$parent_names
      ps <- list()
      
      for (pn in parent_names) {
        levs <- pm$levels[[pn]]
        default_sc <- .default_scores(levs)
        
        scv <- numeric(length(levs))
        for (i in seq_along(levs)) {
          id <- paste0("score_", pn, "_", i)
          val <- input[[id]]
          scv[i] <- if (is.null(val) || !is.finite(val)) default_sc[i] else as.numeric(val)
        }
        ps[[pn]] <- setNames(scv, levs)
      }
      ps
    })
    
    mono_params <- reactive({
      req(cpt(), cpt_type())
      if (cpt_type() != "with_parents") return(NULL)
      
      list(
        method      = input$mono_method %||% "gaussian_center",
        center_temp = as.numeric(input$mono_center_temp %||% 1),
        sigma       = as.numeric(input$mono_sigma %||% 1),
        floor       = as.numeric(input$mono_floor %||% 0.01),
        weights     = current_parent_weights(),
        scores      = current_parent_scores()
      )
    })
    
    # -------------------------
    # Load CPT (runs initially)
    # -------------------------
    observeEvent(input$cpt_file, {
      req(input$cpt_file)
      
      df <- read.csv(file.path(cpt_dir, input$cpt_file), stringsAsFactors = FALSE)
      
      ct <- detect_cpt_type(df)
      sc <- get_state_columns(df)
      
      if (ct == "with_parents" && length(sc) > 0) {
        for (r in seq_len(nrow(df))) df[r, sc] <- snap_probs_to_grid(df[r, sc], step = PROB_STEP, min_prob = PROB_MIN)
      } else if (ct == "no_parents" && "Probability" %in% names(df)) {
        df$Probability <- snap_probs_to_grid(df$Probability, step = PROB_STEP, min_prob = PROB_MIN)
      }
      
      cpt_type(ct)
      state_cols(sc)
      
      cpt_base(df)
      if (ct == "with_parents") {
        st <- infer_cpt_structure(df)
        parent_names <- names(df)[st$parent_cols]
        levels_by_parent <- lapply(parent_names, function(pn) {
          v <- df[[pn]]
          as.character(if (is.factor(v)) levels(v) else unique(v))
        })
        names(levels_by_parent) <- parent_names
        parent_meta(list(parent_names = parent_names, levels = levels_by_parent))
      } else {
        parent_meta(NULL)
      }
      cpt(df)
      
      if (ct == "no_parents") {
        lock_store(rep(FALSE, nrow(df)))
      } else {
        lock_store(matrix(FALSE, nrow = nrow(df), ncol = length(sc)))
      }
      
      current_row(1)
      destroy_observers()
      make_slider_observers()
      
      safe_update_checkbox("mono_on", FALSE)
    }, ignoreInit = FALSE)
    
    # -------------------------
    # Generator recompute (live)
    # -------------------------
    observeEvent(list(
      input$mono_on,
      input$mono_method, input$mono_center_temp, input$mono_sigma, input$mono_floor,
      current_parent_weights(), current_parent_scores()
    ), {
      req(cpt_base(), cpt_type())
      if (!isTRUE(input$mono_on)) return()
      if (cpt_type() != "with_parents") return()
      
      w  <- current_parent_weights()
      ps <- current_parent_scores()
      req(!is.null(w), !is.null(ps))
      
      mono_busy(TRUE)
      on.exit(mono_busy(FALSE), add = TRUE)
      
      base <- cpt_base()
      
      new_df <- build_monotone_cpt(
        base,
        parent_scores  = ps,
        parent_weights = w,
        method         = input$mono_method %||% "gaussian_center",
        temperature    = input$mono_center_temp %||% 1,
        center_temp    = input$mono_center_temp %||% 1,
        sigma          = input$mono_sigma %||% 1,
        floor          = input$mono_floor %||% 0.01
      )
      
      hold_updating_until_flush({
        cpt(new_df)
        state_cols(get_state_columns(new_df))
        
        st <- infer_cpt_structure(new_df)
        lock_store(matrix(FALSE, nrow = nrow(new_df), ncol = length(st$child_cols)))
        current_row(1)
        
        destroy_observers()
        make_slider_observers()
      })
    }, ignoreInit = TRUE)
    
    # -------------------------
    # Row selector
    # -------------------------
    output$row_selector <- renderUI({
      req(cpt(), current_row())
      if (cpt_type() == "with_parents") {
        selectInput("row_id", "Parent configuration (line)",
                    choices = seq_len(nrow(cpt())),
                    selected = current_row()
        )
      }
    })
    
    observeEvent(input$row_id, {
      # user-driven selectInput change
      i <- as.integer(input$row_id)
      if (is.na(i)) return()
      if (identical(i, current_row())) return()
      current_row(i)
    }, ignoreInit = TRUE)
    
    observeEvent(input$table_dblclick_row, {
      # user double-clicks DT row -> update current_row and the selectInput
      if (cpt_type() != "with_parents") return()
      i <- as.integer(input$table_dblclick_row)
      if (is.na(i)) return()
      if (identical(i, current_row())) return()
      
      # Critical: avoid bouncing back through input$row_id observer
      freezeReactiveValue(input, "row_id")
      current_row(i)
      updateSelectInput(session, "row_id", selected = i)
    }, ignoreInit = TRUE)
    
    # keep UI sliders synced when row/cpt changes
    observeEvent(list(current_row(), cpt(), state_cols(), lock_store()), {
      req(cpt(), !is.null(lock_store()))
      
      # If we're already in a controlled update, do nothing.
      # This is key to avoid loops when double-clicking rows or when CPT changes programmatically.
      if (isTRUE(updating())) return()
      
      df <- cpt()
      row <- isolate(current_row())
      sc  <- isolate(state_cols())
      n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(sc)
      req(n_controls > 0)
      
      probs <- if (cpt_type() == "no_parents") as.numeric(df$Probability) else as.numeric(df[row, sc])
      locks <- if (cpt_type() == "no_parents") lock_store() else lock_store()[row, ]
      
      hold_updating_until_flush({
        for (i in seq_along(probs)) {
          safe_update_slider(paste0("p_", i), probs[i])
          safe_update_checkbox(paste0("lock_", i), isTRUE(locks[i]))
        }
      })
    }, ignoreInit = TRUE)
    
    # -------------------------
    # Monotone UI builders
    # -------------------------
    output$mono_weights_ui <- renderUI({
      pm <- parent_meta()
      if (is.null(pm)) return(NULL)
      
      tagList(
        tags$h5("Parent weights"),
        lapply(pm$parent_names, function(pn) {
          id <- paste0("w_", pn)
          numericInput(
            inputId = id,
            label   = pn,
            value   = isolate(input[[id]] %||% 1),
            step    = 0.1
          )
        })
      )
    })
    
    output$mono_scores_ui <- renderUI({
      pm <- parent_meta()
      if (is.null(pm)) return(NULL)
      
      tagList(
        tags$h5("Parent state scores"),
        lapply(pm$parent_names, function(pn) {
          levs <- pm$levels[[pn]]
          default_sc <- .default_scores(levs)
          
          tags$details(
            tags$summary(paste0("Scores for ", pn)),
            lapply(seq_along(levs), function(i) {
              id <- paste0("score_", pn, "_", i)
              numericInput(
                inputId = id,
                label   = levs[i],
                value   = isolate(input[[id]] %||% default_sc[i]),
                step    = 0.5
              )
            })
          )
        })
      )
    })
    
    output$mono_status <- renderUI({
      if (!isTRUE(input$mono_on)) return(NULL)
      if (isTRUE(mono_busy())) {
        tags$div(style = "margin-top: 6px; font-size: 12px; color: #444;", "Recomputing CPT…")
      } else {
        tags$div(style = "margin-top: 6px; font-size: 12px; color: #2E7D32;", "Up to date.")
      }
    })
    
    output$mono_params_text <- renderUI({
      if (!isTRUE(input$mono_on)) return(NULL)
      if (cpt_type() != "with_parents") return(NULL)
      p <- mono_params()
      if (is.null(p)) return(NULL)
      
      tags$div(
        style = "font-size: 12px; margin-top: 6px; padding: 8px; border: 1px solid #ddd; border-radius: 6px;",
        tags$div(tags$b("Current monotone parameters (live)")),
        tags$div(sprintf("Method: %s", p$method)),
        tags$div(sprintf("CenterTemp: %.3g | Sigma: %.3g | Floor: %.3g", p$center_temp, p$sigma, p$floor))
      )
    })
    
    # -------------------------
    # Parent question
    # -------------------------
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
    
    # -------------------------
    # Sliders UI
    # -------------------------
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
        sc <- state_cols()
        req(length(sc) > 0)
        
        row <- current_row()
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
        style = paste0("padding:8px 10px; border-radius:6px; font-weight:600; color:white; background:", col, "; margin-top:6px; margin-bottom:6px;"),
        txt
      )
    })
    
    # -------------------------
    # Heatmap plot (generator ON only)
    # -------------------------
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
    
    output$cpt_curve <- renderPlot({
      if (!isTRUE(input$mono_on)) return(NULL)
      req(cpt(), cpt_type())
      if (cpt_type() != "with_parents") {
        ggplot() + theme_void() + annotate("text", x = 0, y = 0, label = "Heatmap view is for CPTs with parents (P_* columns).")
      } else {
        df <- cpt()
        sc <- state_cols()
        req(length(sc) > 0)
        
        w  <- current_parent_weights()
        ps <- current_parent_scores()
        req(!is.null(w), !is.null(ps))
        
        x_rows <- compute_total_parent_score(df, parent_scores = ps, parent_weights = w)
        
        Pmat <- as.matrix(df[, sc, drop = FALSE])
        storage.mode(Pmat) <- "numeric"
        child_states <- sub("^P_", "", names(df)[sc])
        K <- length(child_states)
        
        long <- as.data.frame(Pmat)
        names(long) <- child_states
        long$.row <- seq_len(nrow(long))
        long$score <- x_rows
        
        long <- long |>
          tidyr::pivot_longer(cols = all_of(child_states), names_to = "child_state", values_to = "prob") |>
          mutate(child_state = factor(child_state, levels = child_states),
                 y = as.integer(child_state)) |>
          filter(is.finite(score), is.finite(prob))
        
        validate(need(nrow(long) > 0, "No valid CPT rows to plot."))
        
        row_scores <- long |>
          distinct(.row, score) |>
          arrange(score) |>
          mutate(score_plot = score)
        
        dup <- duplicated(row_scores$score_plot) | duplicated(row_scores$score_plot, fromLast = TRUE)
        if (any(dup)) {
          rng <- diff(range(row_scores$score_plot))
          eps <- if (is.finite(rng) && rng > 0) rng * 1e-4 else 1e-4
          row_scores <- row_scores |>
            group_by(score_plot) |>
            mutate(score_plot = score_plot + (row_number() - 1) * eps) |>
            ungroup()
        }
        
        s <- row_scores$score_plot
        n <- length(s)
        if (n == 1) {
          xmin <- s - 0.5
          xmax <- s + 0.5
        } else {
          mids <- (s[-1] + s[-n]) / 2
          xmin <- c(s[1] - (mids[1] - s[1]), mids)
          xmax <- c(mids, s[n] + (s[n] - mids[n - 1]))
        }
        
        row_scores$xmin <- xmin
        row_scores$xmax <- xmax
        
        long <- long |>
          left_join(row_scores |> select(.row, xmin, xmax, score_plot), by = ".row")
        
        method     <- input$mono_method %||% "gaussian_center"
        center_tmp <- input$mono_center_temp %||% 1
        sigma_val  <- input$mono_sigma %||% 1
        floor      <- input$mono_floor %||% 0.01
        
        xg <- seq(min(row_scores$score_plot), max(row_scores$score_plot), length.out = 220)
        
        Pfit <- vapply(xg, function(xval) {
          p <- if (method == "ordinal_logit") {
            .probs_ordinal_logit(xval, K, temperature = center_tmp)
          } else {
            .probs_gaussian_center_cs(xval, K, center_temp = center_tmp, sigma = sigma_val)
          }
          p <- pmax(p, 0)
          p <- p + floor
          p / sum(p)
        }, numeric(K))
        Pfit <- t(Pfit)
        
        fit_long <- data.frame(
          x = rep(xg, times = K),
          y = rep(seq_len(K), each = length(xg)),
          prob = as.vector(Pfit)
        )
        
        ggplot() +
          geom_rect(
            data = long,
            aes(xmin = xmin, xmax = xmax, ymin = y - 0.5, ymax = y + 0.5, fill = prob),
            color = NA
          ) +
          { if (isTRUE(input$show_contours)) geom_contour(
            data = fit_long,
            aes(x = x, y = y, z = prob),
            bins = 7,
            linewidth = 0.3,
            alpha = 0.6
          )
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
    
    # -------------------------
    # Bar plot
    # -------------------------
    output$prob_bar <- renderPlot({
      req(cpt())
      df <- cpt()
      
      if (cpt_type() == "no_parents") {
        probs <- as.numeric(df$Probability)
        state_names <- as.character(df$State)
      } else {
        sc <- state_cols()
        req(length(sc) > 0)
        probs <- as.numeric(df[current_row(), sc])
        state_names <- sub("^P_", "", names(df)[sc])
      }
      
      probs[!is.finite(probs)] <- 0
      probs <- pmax(probs, 0)
      s <- sum(probs)
      if (s <= 0) return()
      probs <- probs / s
      
      ramp5 <- c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00")
      cols <- ramp5[round(seq(1, 5, length.out = length(probs)))]
      
      lefts  <- c(0, cumsum(probs)[-length(probs)])
      rights <- cumsum(probs)
      
      op <- par(mar = c(2.2, 1.0, 0.5, 0.5), xaxs = "i", yaxs = "i")
      on.exit(par(op), add = TRUE)
      
      plot(NA, xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Probability mass", ylab = "", yaxt = "n", xaxt = "n")
      
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
      
      axis(1, at = seq(0, 1, by = PROB_STEP), labels = FALSE, tck = -0.03)
      ticks_major <- seq(0, 1, by = 0.1)
      axis(1, at = ticks_major, labels = sprintf("%.1f", ticks_major), cex.axis = 0.8)
    }, res = 120)
    
    # -------------------------
    # ScoreSum + DT builder (single source of truth)
    # -------------------------
    row_score_sum <- reactive({
      req(cpt())
      df <- cpt()
      if (cpt_type() != "with_parents") return(rep(NA_real_, nrow(df)))
      
      w  <- current_parent_weights()
      ps <- current_parent_scores()
      req(!is.null(w), !is.null(ps))
      
      compute_total_parent_score(df, parent_scores = ps, parent_weights = w)
    })
    
    build_dt_df <- function(df) {
      sc <- state_cols()
      
      prob_names <- character(0)
      if (cpt_type() == "with_parents" && length(sc) > 0) prob_names <- names(df)[sc]
      if (cpt_type() == "no_parents" && "Probability" %in% names(df)) prob_names <- "Probability"
      
      df_show <- df
      
      if (cpt_type() == "with_parents" && length(sc) > 0) {
        insert_pos <- min(sc)
        score_vec <- row_score_sum()
        df_show <- cbind(
          df_show[, seq_len(insert_pos - 1), drop = FALSE],
          ScoreSum = round(score_vec, 3),
          df_show[, insert_pos:ncol(df_show), drop = FALSE]
        )
      }
      
      if (length(prob_names) > 0) df_show[, prob_names] <- round(df_show[, prob_names], 2)
      
      df_show <- sanitize_for_dt(df_show)
      names(df_show) <- make.unique(names(df_show))
      
      bad <- row_bad_flags(df = df, ct = cpt_type(), sc = sc)
      
      cbind(
        .row_id = seq_len(nrow(df_show)),
        .row_bad = as.integer(bad),
        df_show
      )
    }
    
    output$cpt_table <- renderDT({
      req(cpt())
      df_show <- build_dt_df(cpt())
      
      datatable(
        df_show,
        selection = "none",
        rownames = FALSE,
        options = list(
          pageLength = 100,
          columnDefs = list(list(targets = c(0, 1), visible = FALSE, searchable = FALSE)),
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
    
    observeEvent(list(cpt(), input$cpt_table_rows_current, row_score_sum()), {
      req(cpt(), !is.null(input$cpt_table_rows_current))
      df_show <- build_dt_df(cpt())
      DT::replaceData(cpt_proxy, df_show, resetPaging = FALSE, rownames = FALSE)
    }, ignoreInit = TRUE)
    
    # -------------------------
    # Save / Quit
    # -------------------------
    observeEvent(input$save, {
      req(cpt(), input$cpt_file)
      df <- cpt()
      sc <- state_cols()
      
      if (length(sc) > 0) df[, sc] <- round(df[, sc], 2)
      
      write.csv(df, file.path(pathProCpt, folder, input$cpt_file), row.names = FALSE)
      
      cpt(df)
      cpt_base(df)
      showNotification("CPT saved", type = "message")
    })
    
    observeEvent(input$quit, stopApp())
  }
  
  shinyApp(ui, server)
}

if (interactive()) {
  run_cpt_app(pathProCpt = pathProCpt, folder = "01_CPT")
}