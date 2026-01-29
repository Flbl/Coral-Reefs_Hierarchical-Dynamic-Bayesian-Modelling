#############################################################################################
#                                                                                           #
#                                                                                           #
#####             Simple Shiny interface to modify CPTs with mouse cursors              #####
#                                                                                           #
#                                                                                           #
#############################################################################################



library(shiny)
library(DT)

# dir.create(file.path(pathProCpt,"00_ManualShinyExportsTest"), showWarnings = FALSE)

# folder = "01_CPT"


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


PROB_STEP <- 0.01
TOTAL_UNITS <- round(1 / PROB_STEP)  # 100


# The renormalise function for keeping 0-1 between states when modifying a cursor
renormalise <- function(probs, locked, changed, new_value,
                        step = PROB_STEP) {
  
  n <- length(probs)
  locked <- as.logical(locked)
  
  total_units <- round(1 / step)
  
  # Convert probabilities to integer units
  units <- round(probs / step)
  
  # Units locked elsewhere
  locked_others <- locked & seq_len(n) != changed
  used_units_locked <- sum(units[locked_others])
  
  # Max units allowed for changed index
  max_units_changed <- max(0, total_units - used_units_locked)
  
  # Clamp changed value to feasible units
  units[changed] <- min(
    max(round(new_value / step), 0),
    max_units_changed
  )
  
  # Remaining units to distribute
  remaining_units <- total_units - used_units_locked - units[changed]
  
  # Free indices
  free_idx <- which(!locked & seq_len(n) != changed)
  
  if (length(free_idx) == 0 || remaining_units <= 0) {
    return(units * step)
  }
  
  # Current weights (use previous units, but non-negative)
  w <- pmax(units[free_idx], 0)
  sw <- sum(w)
  
  if (sw == 0) {
    # Uniform integer allocation
    base <- remaining_units %/% length(free_idx)
    extra <- remaining_units %% length(free_idx)
    
    units[free_idx] <- base
    if (extra > 0) {
      units[free_idx][seq_len(extra)] <-
        units[free_idx][seq_len(extra)] + 1
    }
  } else {
    # Proportional allocation in integer units
    alloc <- floor(w / sw * remaining_units)
    remainder <- remaining_units - sum(alloc)
    
    units[free_idx] <- alloc
    
    if (remainder > 0) {
      # Distribute remaining units deterministically
      order_idx <- order(-w)
      units[free_idx][order_idx[seq_len(remainder)]] <-
        units[free_idx][order_idx[seq_len(remainder)]] + 1
    }
  }
  
  # Convert back to probabilities
  units * step
}


# Monotone CPT generation
# Safe cleaning for display / matching
.clean_txt <- function(x) trimws(gsub("_+", " ", as.character(x)))


# Default score vector for an ordered factor/character vector of levels
# Example: 5 levels -> -2,-1,0,1,2 ; 3 levels -> -1,0,1 ; 4 levels -> -1.5,-0.5,0.5,1.5
.default_scores <- function(levels_vec) {
  m <- length(levels_vec)
  seq(-(m - 1) / 2, (m - 1) / 2, length.out = m)
}


# Extract CPT structure: parent columns and child probability columns (P_*), and child state names
infer_cpt_structure <- function(df) {
  child_cols <- grep("^P_", names(df))
  if (length(child_cols) == 0) {
    stop("Expected P_* columns for a with-parents CPT.")
  }
  parent_cols <- setdiff(seq_along(df), child_cols)
  
  child_states <- sub("^P_", "", names(df)[child_cols])
  list(parent_cols = parent_cols, child_cols = child_cols, child_states = child_states)
}


# /!\ A RULE FOR BELOW ORDINAL/GAUSSIAN FITTING TO WORK:
# CPTs HAVE TO BE ORDERED IN THEIR STATE COMBINATION 
# (1: its easier for expert elicitation, 2: The function makes a curve inside the table so it has to be ordered)

# Method 1: ordinal cumulative logit (proportional odds)
# Higher score -> distribution shifts stochastically to higher categories
.probs_ordinal_logit <- function(score, K, temperature = 1, cut_span = 2.0) {
  # temperature > 0 ; smaller => sharper
  temperature <- max(1e-6, temperature)
  
  # cutpoints (K-1) equally spaced
  theta <- seq(-cut_span, cut_span, length.out = K - 1)
  
  # latent mean increases with score
  mu <- score / temperature
  
  # cumulative probabilities P(Y <= k)
  cum <- plogis(theta - mu)
  
  p <- numeric(K)
  p[1] <- cum[1]
  if (K > 2) {
    for (k in 2:(K - 1)) p[k] <- cum[k] - cum[k - 1]
  }
  p[K] <- 1 - cum[K - 1]
  pmax(p, 0)
}


# Method 2: Gaussian-like kernel around a score-dependent "center" index
# Center moves monotonically with score; sigma controlled by temperature
.probs_gaussian_center <- function(score, K, temperature = 1) {
  temperature <- max(1e-6, temperature)
  
  # Map score -> center in [1, K] via logistic
  center <- 1 + (K - 1) * plogis(score / temperature)
  
  # sigma: smaller temperature -> sharper peak
  sigma <- max(1e-3, temperature)
  
  idx <- seq_len(K)
  w <- exp(-0.5 * ((idx - center) / sigma)^2)
  w / sum(w)
}



# Main builder
# df: existing CPT data.frame (with parents, wide P_* columns)
# parent_scores: named list, one entry per parent column name. Each entry is a named numeric vector mapping state -> score
# parent_weights: named numeric vector, per parent column name
# method: "ordinal_logit" or "gaussian_center" or "coarse"
# temperature: >0 ; smaller => sharper
# floor: minimum probability mass per child state (smoothness)
# returns: df with same parent columns + regenerated P_* columns
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
  child_states <- st$child_states
  K <- length(child_cols)
  
  parent_names <- names(df)[parent_cols]
  
  # Defaults
  if (is.null(parent_weights)) {
    parent_weights <- setNames(rep(1, length(parent_names)), parent_names)
  } else {
    # ensure all parents present
    missing <- setdiff(parent_names, names(parent_weights))
    if (length(missing) > 0) parent_weights[missing] <- 1
    parent_weights <- parent_weights[parent_names]
  }
  
  if (is.null(parent_scores)) {
    # default based on factor levels if available, else unique order in data
    parent_scores <- lapply(parent_names, function(pn) {
      v <- df[[pn]]
      levs <- if (is.factor(v)) levels(v) else unique(as.character(v))
      levs <- as.character(levs)
      sc <- .default_scores(levs)
      setNames(sc, levs)
    })
    names(parent_scores) <- parent_names
  } else {
    # ensure all parents present
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
  
  
  # Build total score per row
  total_score <- numeric(nrow(df))
  for (pn in parent_names) {
    w <- as.numeric(parent_weights[pn])
    mp <- parent_scores[[pn]]
    vals <- as.character(df[[pn]])
    sc <- mp[vals]
    sc[is.na(sc)] <- 0
    total_score <- total_score + w * as.numeric(sc)
  }
  
  # Generate probs row-wise
  P <- matrix(0, nrow(df), K)
  for (i in seq_len(nrow(df))) {
    s <- total_score[i]
    if (method == "ordinal_logit") {
      p <- .probs_ordinal_logit(s, K, temperature = temperature)
    } else {
      p <- .probs_gaussian_center(s, K, temperature = temperature)
    }
    
    # smoothness floor, renormalize
    p <- pmax(p, 0)
    p <- p + floor
    p <- p / sum(p)
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
  # cpt_patterns: NULL = no filter
  # or character vector of regex patterns (OR-ed)
  # example: cpt_patterns = c("CPT__Site_", "CPT__Env_")

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
        tags$h4("Monotone CPT generator"),
        uiOutput("mono_weights_ui"),
        uiOutput("mono_scores_ui"),
        selectInput("mono_method", "Method",
                    choices = c("Ordinal logit" = "ordinal_logit",
                                "Gaussian center" = "gaussian_center"),
                    selected = "gaussian_center"),
        sliderInput("mono_temp", "Temperature / sharpness",
                    min = 0.1, max = 5, value = 1, step = 0.1),
        # sliderInput(
        #   "mono_cutspan", "Cut-span (ordinal logit separation)",
        #   min = 0.5, max = 6, value = 3, step = 0.1
        # ),
        numericInput("mono_floor", "Smoothness floor", value = 0.01, min = 0, step = 0.01),
        actionButton("apply_monotone", "Apply monotone CPT", class = "btn-primary"),
        hr(),
        actionButton("save", "Save CPT"),
        actionButton("quit", "Quit"),
        tags$style(HTML("
          tr.active-row td {
            background-color: #5485b0 !important;
          }
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
    
    cpt_dir <- file.path(pathProCpt, folder) #"Data/01_Processed/CPT" #pathProCpt # TO CHANGE TO MAKE PROMPT FROM OTHER SCRIPTS
    if (!dir.exists(cpt_dir)) {
      stop(paste("CPT directory does not exist:", cpt_dir))
    }
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
      
      # keep current selection if still available
      sel <- isolate(input$cpt_file)
      if (is.null(sel) || !(sel %in% files)) sel <- files[1]
      
      updateSelectInput(session, "cpt_file", choices = files, selected = sel)
    })
    
    
    cpt <- reactiveVal(NULL)
    cpt_type <- reactiveVal("with_parents") #safe default
    state_cols <- reactiveVal(integer(0))
    # locked <- reactiveVal(logical(0))
    lock_store <- reactiveVal(NULL)   # will hold a matrix (rows x states) or a vector
    obs_inited <- reactiveVal(FALSE)  # ensures we don't create duplicate observers
    current_row <- reactiveVal(1)
    
    
    # Making the active row available to javascript
    observe({
      session$sendCustomMessage("activeRow", list(row = current_row()))
    })
    
    
    
    # ---- Monotone UI: weights + scores (only when CPT has parents)
    output$mono_weights_ui <- renderUI({
      req(cpt())
      req(cpt_type())
      if (cpt_type() != "with_parents") return(NULL)
      
      st <- infer_cpt_structure(cpt())
      parent_cols <- st$parent_cols
      parent_names <- names(cpt())[parent_cols]
      
      tagList(
        tags$h5("Parent weights"),
        lapply(parent_names, function(pn) {
          numericInput(paste0("w_", pn), label = pn, value = 1, step = 0.1)
        })
      )
    })
    
    output$mono_scores_ui <- renderUI({
      req(cpt())
      req(cpt_type())
      if (cpt_type() != "with_parents") return(NULL)
      
      df <- cpt()
      st <- infer_cpt_structure(df)
      parent_cols <- st$parent_cols
      parent_names <- names(df)[parent_cols]
      
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
      req(cpt())
      req(cpt_type())
      if (cpt_type() != "with_parents") {
        showNotification("Monotone generator applies to CPTs with parents (P_* columns).", type = "warning")
        return()
      }
      
      df <- cpt()
      st <- infer_cpt_structure(df)
      parent_cols <- st$parent_cols
      parent_names <- names(df)[parent_cols]
      
      # collect weights
      w <- setNames(rep(1, length(parent_names)), parent_names)
      for (pn in parent_names) {
        val <- input[[paste0("w_", pn)]]
        if (!is.null(val)) w[pn] <- as.numeric(val)
      }
      
      # collect scores
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
      
      # Replace CPT + reset locks for safety
      cpt(new_df)
      state_cols(get_state_columns(new_df)) # unchanged but safe
      lock_store(matrix(FALSE, nrow = nrow(new_df), ncol = length(st$child_cols)))
      obs_inited(FALSE)
      showNotification("Monotone CPT applied. You can now tweak row-wise.", type = "message")
    })
    
    
    # ---- Load CPT ----
    observeEvent(input$cpt_file, {
      req(input$cpt_file)
      
      df <- read.csv(
        file.path(cpt_dir, input$cpt_file),
        stringsAsFactors = FALSE
      )
      
      cpt(df)
      cpt_type(detect_cpt_type(df))
      state_cols(get_state_columns(df))
      # locked(rep(FALSE, length(state_cols())))
      if (cpt_type() == "no_parents") {
        lock_store(rep(FALSE, nrow(df)))           # one lock per row-state
      } else {
        n_states <- length(state_cols())
        lock_store(matrix(FALSE, nrow = nrow(df), ncol = n_states))
      }
      obs_inited(FALSE)
      current_row(1)
    })
    
    # ---- Row selector (only if parents exist) ----
    output$row_selector <- renderUI({
      req(cpt())
      req(current_row())
      
      if (cpt_type() == "with_parents") {
        numeric_cols <- state_cols()
        parent_cols <- setdiff(seq_along(cpt()), numeric_cols)
        
        selectInput(
          "row_id",
          "Parent configuration (line)",
          choices = seq_len(nrow(cpt())),
          selected = current_row()
        )
      }
    })
    
    observeEvent(input$row_id, {
      current_row(as.integer(input$row_id))
    })
    
    observeEvent(list(current_row(), cpt(), state_cols(), lock_store()), {
      req(cpt())
      # req(length(state_cols()) > 0)
      req(!is.null(lock_store()))
      
      df <- cpt()
      row <- current_row()
      sc  <- state_cols()
      
      n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(state_cols())
      req(n_controls > 0)
      
      # pull probabilities for the row / vector
      probs <- if (cpt_type() == "no_parents") {
        as.numeric(df$Probability)
      } else {
        as.numeric(df[row, sc])
      }
      
      # pull locks for this row
      locks <- if (cpt_type() == "no_parents") {
        lock_store()
      } else {
        lock_store()[row, ]
      }
      
      updating(TRUE)
      on.exit(updating(FALSE), add = TRUE)
      
      # update UI controls to reflect selected row
      for (i in seq_along(probs)) {
        updateSliderInput(session, paste0("p_", i), value = probs[i])
        updateCheckboxInput(session, paste0("lock_", i), value = isTRUE(locks[i]))
      }
    }, ignoreInit = TRUE)
    
    # double click update of the row for parent configuration
    observeEvent(input$table_dblclick_row, {
      req(cpt_type())
      if (cpt_type() != "with_parents") return()
      
      i <- as.integer(input$table_dblclick_row)
      if (is.na(i)) return()
      
      current_row(i)
      
      # keep dropdown selector in sync
      updateSelectInput(session, "row_id", selected = i)
    })
    
    # Parent question
    output$parent_question <- renderUI({
      req(cpt())
      req(cpt_type())
      req(input$cpt_file)
      
      df <- cpt()
      
      # Helper: clean display text
      clean <- function(x) trimws(gsub("_+", " ", as.character(x)))
      
      # Infer child node name from filename (e.g., CPT__Site_.csv -> Site)
      child_node <- gsub("\\.csv$", "", input$cpt_file)
      child_node <- gsub("^CPT_+", "", child_node)
      child_node <- clean(child_node)
      if (nchar(child_node) == 0) child_node <- "the node"
      
      # No-parents CPT: just show a simple prompt
      if (cpt_type() == "no_parents") {
        return(tags$div(
          style = "margin-top: 6px; margin-bottom: 6px;",
          tags$span("What is the probability that "),
          tags$span(style = "font-weight: bold;", paste0("'", child_node, "'")),
          tags$span(" is:")
        ))
      }
      
      # With-parents CPT
      sc <- state_cols()
      row <- current_row()
      
      # Parent columns = all non-probability columns
      parent_cols <- setdiff(seq_along(df), sc)
      req(length(parent_cols) > 0)
      
      parent_names <- clean(names(df)[parent_cols])
      # Secondary clean for parents for t-n
      parent_names <- gsub(".","-", parent_names, fixed = TRUE)
      parent_vals  <- clean(df[row, parent_cols, drop = TRUE])
      
      # --- Build a per-parent color map (unique values per parent column)
      # pal <- palette.colors(palette = "Okabe-Ito")
      pal <- c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00")
      
      # For each parent column, map its unique values to palette colors
      # parent_color_maps <- lapply(parent_cols, function(j) {
      #   vals <- clean(df[[j]])
      #   u <- unique(vals)
      #   cols <- rep(pal, length.out = length(u))
      #   stats::setNames(cols, u)
      # })
      parent_color_maps <- lapply(parent_cols, function(j) {
        vals <- df[[j]]
        # if it's a factor, use its defined level order; else use unique order
        u_ord <- if (is.factor(vals)) levels(vals) else unique(clean(vals))
        u_ord <- clean(u_ord)
        
        cols <- pal[round(seq(1, length(pal), length.out = length(u_ord)))]
        stats::setNames(cols, u_ord)
      })
      
      # Build condition fragments: "[Parent] is [Value]" with value colored
      cond_fragments <- lapply(seq_along(parent_cols), function(k) {
        p_name <- parent_names[k]
        p_val  <- parent_vals[k]
        cmap   <- parent_color_maps[[k]]
        colhex <- unname(cmap[p_val])
        
        tags$span(
          style = "font-weight: bold;",
          tags$span(paste0(p_name, " is ")),
          tags$span(
            style = paste0("color:", colhex, ";"),
            paste0(p_val)
          )
        )
      })
      
      # Interleave " and " between fragments
      cond_with_and <- list()
      for (k in seq_along(cond_fragments)) {
        cond_with_and <- c(cond_with_and, list(cond_fragments[[k]]))
        if (k < length(cond_fragments)) {
          cond_with_and <- c(cond_with_and, list(tags$span(" and ")))
        }
      }
      
      tags$div(
        style = "margin-top: 6px; margin-bottom: 6px;",
        tags$span("When "),
        cond_with_and,
        tags$span(", what is the probability that "),
        tags$span(style = "font-style: italic;", paste0("'", child_node, "'")),
        tags$span(" is:")
      )
    })
    
    
    # ---- Sliders ----
    output$sliders_ui <- renderUI({
      req(cpt())
      req(cpt_type())
      req(!is.null(lock_store()))
      
      df <- cpt()
      
      if (cpt_type() == "no_parents") {
        probs <- as.numeric(df$Probability)
        locks <- as.logical(lock_store())
        state_names <- as.character(df$State)
        
        lapply(seq_along(probs), function(i) {
          tags$div(
            style = "margin-bottom: 10px; padding-bottom: 6px; border-bottom: 1px solid #eee;",
            
            # State title
            tags$div(
              style = "font-weight: bold; margin-bottom: 4px; text-align: center;",
              gsub("_"," ",state_names[i])
            ),
            
            # Lock checkbox (label is always "Lock")
            checkboxInput(
              paste0("lock_", i),
              label = "Lock",
              value = isTRUE(locks[i])
            ),
            
            # Slider
            sliderInput(
              paste0("p_", i),
              label = NULL,
              min = 0, max = 1, step = PROB_STEP,
              value = probs[i]
            )
          )
        })
        
      } else {
        # req(length(state_cols()) > 0)
        
        n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(state_cols())
        req(n_controls > 0)
        
        row <- current_row()
        sc <- state_cols()
        
        probs <- as.numeric(df[row, sc])
        locks <- as.logical(lock_store()[row, ])
        state_names <- sub("^P_", "", names(df)[sc])
        
        lapply(seq_along(probs), function(i) {
          tags$div(
            style = "margin-bottom: 10px; padding-bottom: 6px; border-bottom: 1px solid #eee;",
            
            # State title
            tags$div(
              style = "font-weight: bold; margin-bottom: 4px; text-align: center;",
              gsub("_"," ",state_names[i])
            ),
            
            # Lock checkbox (label is always "Lock")
            checkboxInput(
              paste0("lock_", i),
              label = "Lock",
              value = isTRUE(locks[i])
            ),
            
            # Slider
            sliderInput(
              paste0("p_", i),
              label = NULL,
              min = 0, max = 1, step = PROB_STEP,
              value = probs[i]
            )
          )
        })
      }
    })
    
    
    
    # ---- Slider observer ----
    observeEvent(list(cpt(), state_cols()), {
      req(cpt())
      req(cpt_type())
      # req(length(state_cols()) > 0)
      req(!is.null(lock_store()))
      
      df <- cpt()
      n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(state_cols())
      req(n_controls > 0)
      
      # Do not re-register observers repeatedly
      if (obs_inited()) return()
      obs_inited(TRUE)
      
      n_controls <- if (cpt_type() == "no_parents") nrow(cpt()) else length(state_cols())
      for (i in seq_len(n_controls)) {
        local({
          idx <- i
          
          observeEvent(input[[paste0("p_", idx)]], {
            if (updating()) return()
            
            df <- cpt()
            row <- current_row()
            sc  <- state_cols()
            
            # read current lock state from UI
            n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(sc)
            lock_vals <- vapply(
              seq_len(n_controls),
              function(k) isTRUE(input[[paste0("lock_", k)]]),
              logical(1)
            )
            
            # persist locks for this row
            ls <- lock_store()
            if (cpt_type() == "no_parents") {
              lock_store(lock_vals)
            } else {
              ls[row, ] <- lock_vals
              lock_store(ls)
            }
            
            # current probs for this row
            old_probs <- if (cpt_type() == "no_parents") {
              as.numeric(df$Probability)
            } else {
              as.numeric(df[row, sc])
            }
            
            sum_locked_other <- sum(old_probs[lock_vals & seq_along(old_probs) != idx])
            max_allowed <- max(0, 1 - sum_locked_other)
            
            new_value <- min(max(as.numeric(input[[paste0("p_", idx)]]), 0), max_allowed)
            
            # If user dragged beyond feasible region, snap ONLY this slider to max_allowed
            if (!isTRUE(all.equal(new_value, as.numeric(input[[paste0("p_", idx)]])))) {
              updating(TRUE)
              on.exit(updating(FALSE), add = TRUE)
              updateSliderInput(session, paste0("p_", idx), value = new_value)
            }
            
            # new_value <- as.numeric(input[[paste0("p_", idx)]])
            
            updating(TRUE)
            on.exit(updating(FALSE), add = TRUE)
            
            new_probs <- renormalise(
              probs     = old_probs,
              locked    = lock_vals,
              changed   = idx,
              new_value = new_value
            )
            
            # write back only to the current row
            if (cpt_type() == "no_parents") {
              df$Probability <- new_probs
            } else {
              df[row, sc] <- new_probs
            }
            cpt(df)
            
            # update other sliders only
            for (j in seq_along(new_probs)) {
              if (j != idx) {
                updateSliderInput(session, paste0("p_", j), value = new_probs[j])
              }
            }
          }, ignoreInit = TRUE)
        })
      }
    }, ignoreInit = TRUE)
    
    
    observeEvent({
      req(cpt_type())
      df <- cpt()
      n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(state_cols())
      lapply(seq_len(n_controls), function(i) input[[paste0("lock_", i)]])
    }, {
      req(cpt_type())
      req(!is.null(lock_store()))
      df <- cpt()
      n_controls <- if (cpt_type() == "no_parents") nrow(df) else length(state_cols())
      row <- current_row()
      
      lock_vals <- vapply(
        seq_len(n_controls),
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
    
    
    observeEvent(current_row(), {
      req(!is.null(lock_store()))
      if (cpt_type() == "with_parents") {
        ls <- lock_store()
        ls[current_row(), ] <- FALSE
        lock_store(ls)
      } else {
        df <- cpt()
        lock_store(rep(FALSE, nrow(df)))
      }
    }, ignoreInit = TRUE)
    
    
    # ---- bar plot ----
    output$prob_bar <- renderPlot({
      req(cpt())
      req(cpt_type())
      req(length(state_cols()) > 0)
      req(!is.null(lock_store()))
      
      df <- cpt()
      sc <- state_cols()
      row <- current_row()
      
      probs <- if (cpt_type() == "no_parents") {
        as.numeric(df$Probability)
      } else {
        as.numeric(df[row, sc])
      }
      
      locks <- if (cpt_type() == "no_parents") {
        as.logical(lock_store())
      } else {
        as.logical(lock_store()[row, ])
      }
      
      # State names (strip P_)
      if (cpt_type() == "no_parents") {
        state_names <- as.character(df$State)
      } else {
        state_names <- names(df)[sc]
        state_names <- sub("^P_", "", state_names)
      }
      
      # Defensive normalization for plotting only
      probs <- pmax(probs, 0)
      s <- sum(probs)
      if (s > 0) probs <- probs / s
      
      # Accessible sequential ramp (good -> bad), Okabe-Ito inspired anchors
      ramp5 <- c("#0072B2", "#56B4E9", "#F0E442", "#E69F00", "#D55E00")
      
      n_states <- length(probs)
      ramp <- ramp5[round(seq(1, 5, length.out = n_states))]
      cols <- ramp  # <- use column order as-is (good -> bad)
      
      # Compute cumulative positions for stacked bar
      lefts  <- c(0, cumsum(probs)[-length(probs)])
      rights <- cumsum(probs)
      
      # Small margins for small height
      op <- par(mar = c(2.2, 1.0, 0.5, 0.5), xaxs = "i", yaxs = "i")
      on.exit(par(op), add = TRUE)
      
      plot(NA,
           xlim = c(0, 1), ylim = c(0, 1),
           xlab = "Probability mass", ylab = "",
           yaxt = "n", xaxt = "n")
      
      # Label scaling parameters
      min_width_label <- 0.05
      cex_max <- 0.85
      cex_min <- 0.55
      
      for (i in seq_along(probs)) {
        rect(
          xleft = lefts[i], ybottom = 0.25,
          xright = rights[i], ytop = 0.75,
          col = cols[i],
          border = "black",
          lwd = if (isTRUE(locks[i])) 3 else 1
        )
        
        width <- rights[i] - lefts[i]
        if (width >= min_width_label) {
          # width=0.05 -> cex_min ; width=0.20+ -> cex_max
          scaled <- cex_min + (pmin(width, 0.20) - 0.05) * (cex_max - cex_min) / (0.20 - 0.05)
          cex <- max(cex_min, min(cex_max, scaled))
          
          mid <- (lefts[i] + rights[i]) / 2
          text(mid, 0.50, labels = sprintf("%s\n%.2f", gsub("_"," ",state_names[i]), probs[i]), cex = cex)
        }
      }
      
      # Ticks at elicitation grid
      ticks_minor <- seq(0, 1, by = PROB_STEP)
      ticks_major <- seq(0, 1, by = 0.1)
      axis(1, at = ticks_minor, labels = FALSE, tck = -0.03)
      axis(1, at = ticks_major, labels = sprintf("%.1f", ticks_major), cex.axis = 0.8)
      
    }, res = 120)
    
    
    
    # ---- Table ----
    output$cpt_table <- renderDT({
      req(cpt())
      df <- cpt()
      
      # columns to round
      sc <- state_cols()
      prob_names <- character(0)
      if (cpt_type() == "with_parents" && length(sc) > 0) prob_names <- names(df)[sc]
      if (cpt_type() == "no_parents" && "Probability" %in% names(df)) prob_names <- "Probability"
      
      # make a rounded DISPLAY copy (doesn't touch your cpt())
      df_show <- df
      if (length(prob_names) > 0) {
        df_show[, prob_names] <- round(df_show[, prob_names], 2)
      }
      
      active <- current_row()  # << this is the key
      
      dt <- datatable(
        df_show,
        selection = "none",
        options = list(
          pageLength = 100,
          
          # dblclick row -> send row to Shiny
          rowCallback = JS("
        function(row, data, index) {
          $(row).off('dblclick');
          $(row).on('dblclick', function() {
            Shiny.setInputValue('table_dblclick_row', index + 1, {priority: 'event'});
          });
        }
      "),
          
          # highlight active row; embed active row value from R
          drawCallback = JS(sprintf("
        function(settings) {
          var api = new $.fn.dataTable.Api(settings);
          var active = %d;

          api.rows().every(function(rowIdx) {
            var tr = this.node();
            if ((rowIdx + 1) === active) $(tr).addClass('active-row');
            else $(tr).removeClass('active-row');
          });
        }
      ", active))
        )
      )
      
      dt
    })
    
    
    # ---- Save ----
    observeEvent(input$save, {
      req(cpt())
      df <- cpt()
      sc <- state_cols()
      
      if (length(sc) > 0) {
        df[, sc] <- round(df[, sc], 2)
      }
      
      write.csv(df, file.path(pathProCpt,folder, input$cpt_file), row.names = FALSE)
      # write.csv(df, file.path(pathProCpt,"00_ManualShinyExportsTest", input$cpt_file), row.names = FALSE)
      cpt(df)  # keep app state consistent with what you saved
      
      showNotification("CPT saved", type = "message")
    })
    
    # ---- Quit ----
    observeEvent(input$quit, {
      stopApp()
    })
  }
  
  shinyApp(ui, server)

}


# Direct runner

# If you run this file interactively, it starts the app with no filter
if (interactive()) {
  run_cpt_app(pathProCpt = pathProCpt, folder = "01_CPT")
}

