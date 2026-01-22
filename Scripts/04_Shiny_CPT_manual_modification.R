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

PROB_STEP <- 0.05
TOTAL_UNITS <- round(1 / PROB_STEP)  # 20

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


# =========================
# UI
# =========================

ui <- fluidPage(
  
  titlePanel("Bayesian Network CPT Editor"),
  
  sidebarLayout(
    sidebarPanel(
      selectInput("cpt_file", "Select CPT", choices = NULL),
      uiOutput("row_selector"),
      hr(),
      uiOutput("sliders_ui"),
      hr(),
      actionButton("save", "Save CPT"),
      actionButton("quit", "Quit")
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
  
  cpt_dir <- "Data/01_Processed/CPT" #pathProCpt
  if (!dir.exists(cpt_dir)) {
    stop(paste("CPT directory does not exist:", cpt_dir))
  }
  updating <- reactiveVal(FALSE)
  
  observe({
    files <- list.files(cpt_dir, pattern = "\\.csv$")
    if (length(files) == 0) return()
    
    updateSelectInput(
      session,
      "cpt_file",
      choices = files,
      selected = files[1]
    )
  })
  
  cpt <- reactiveVal(NULL)
  cpt_type <- reactiveVal("with_parents") #safe default
  state_cols <- reactiveVal(integer(0))
  # locked <- reactiveVal(logical(0))
  lock_store <- reactiveVal(NULL)   # will hold a matrix (rows x states) or a vector
  obs_inited <- reactiveVal(FALSE)  # ensures we don't create duplicate observers
  current_row <- reactiveVal(1)
  
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
        "Parent configuration",
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
        tagList(
          checkboxInput(
            paste0("lock_", i),
            label = state_names[i],
            value = isTRUE(locks[i])
          ),
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
        tagList(
          checkboxInput(
            paste0("lock_", i),
            label = state_names[i],
            value = isTRUE(locks[i])
          ),
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
        text(mid, 0.50, labels = sprintf("%s\n%.2f", state_names[i], probs[i]), cex = cex)
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
    
    sc <- state_cols()
    if (length(sc) == 0) {
      return(datatable(df, options = list(pageLength = 10)))
    }
    
    prob_names <- names(df)[sc]
    
    dt <- datatable(df, options = list(pageLength = 10))
    DT::formatRound(dt, columns = prob_names, digits = 2)
  })
  
  # ---- Save ----
  observeEvent(input$save, {
    req(cpt())
    df <- cpt()
    sc <- state_cols()
    
    if (length(sc) > 0) {
      df[, sc] <- round(df[, sc], 2)
    }
    
    write.csv(df, file.path(cpt_dir, input$cpt_file), row.names = FALSE)
    cpt(df)  # keep app state consistent with what you saved
    
    showNotification("CPT saved", type = "message")
  })
  
  # ---- Quit ----
  observeEvent(input$quit, {
    stopApp()
  })
}

shinyApp(ui, server)


df1 <- read.csv(file.path(pathProCpt, "CPT__Site_.csv"), header = TRUE)
df1

df2 <- read.csv(file.path(pathProCpt, "CPT_Coral_Reef_Ecosystem_Prey_Availability.csv"), header = TRUE)
df2

df3 <- read.csv(file.path(pathProCpt, "CPT_Env_Temperature_Site_XX.csv"), header = TRUE)
df3


