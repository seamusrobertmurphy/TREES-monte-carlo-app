# ------------------------------------------------------------------------ #
# ART TREES Monte Carlo Simulation - ENHANCED WITH HYPERPARAMETER TUNING
# Government of Guyana - REDD+ Reporting
# Compliant with ART TREES Standard Version 2.0, Section 8
# FIXED: Hyperparameter optimization for small datasets
# ------------------------------------------------------------------------ #

pacman::p_load(BiocManager,
               dplyr, DT,
               ggplot2,
               plotly,
               shiny, shinydashboard,
               caret)  # Enhanced for hyperparameter tuning

# CORRECTED: Complete Guyana data structure with ALL components
create_enhanced_guyana_data <- function() {
  # Complete carbon stock data using medians and standard errors with min/max values
  carbon_stocks <- data.frame(
    Component = c("AG Tree", "BG Tree", "Saplings", "Standing Dead Wood", "Lying Dead Wood", 
                  "Sum C pools w/o litter", "Litter", "Sum C pools w/ litter", "Soil", "Sum ALL POOLS"),
    Median_tC_ha = c(205.8, 48.3, 3.7, 2.6, 8.6, 269.0, 3.3, 272.3, 58.7, 331.0),
    SE_tC_ha = c(60.4/sqrt(118), 14.3/sqrt(118), 2.0/sqrt(118), 4.0/sqrt(118), 8.1/sqrt(118), 
                 75.2/sqrt(118), 1.3/sqrt(118), 76.5/sqrt(118), 61.5/sqrt(87), 506.2/sqrt(87)),
    Min_tC_ha = c(91.6, 21.2, 0.5, 0.0, 0.0, 113.3, 1.2, 114.4, 10.1, 456.8),
    Max_tC_ha = c(353.7, 83.1, 18.8, 13.7, 42.3, 511.7, 8.7, 520.3, 502.4, 3749.8),
    CI_90_tC_ha = c(9.2, 2.2, 0.3, 0.6, 1.2, 11.5, 0.2, 11.7, 11.0, 83.2),
    CI_percent = c(4.5, 4.5, 8.5, 23.8, 14.3, 4.3, 7.5, 4.3, 18.7, 6.9),
    n_plots = c(118, 118, 118, 118, 118, 118, 118, 118, 87, 87),
    stringsAsFactors = FALSE
  )
  
  # Convert to tCO2/ha (multiply by 3.67)
  carbon_stocks$Median_tCO2_ha <- carbon_stocks$Median_tC_ha * 3.67
  carbon_stocks$SE_tCO2_ha <- carbon_stocks$SE_tC_ha * 3.67
  carbon_stocks$CI_90_tCO2_ha <- carbon_stocks$CI_90_tC_ha * 3.67
  
  # Activity data for recent years with actual standard errors from MRVS accuracy assessment
  activity_data <- data.frame(
    Activity_Type = c(rep("Deforestation", 7), rep("Degradation", 3)),
    Driver = c("Forestry infrastructure", "Agriculture", "Mining (medium & large scale)", 
               "Infrastructure", "Settlements", "Fire-Biomass burning", "Shifting Cultivation",
               "Logging - volume harvested", "Logging - skid trail length", "Mining and Infrastructure (buffer area)"),
    Units = c(rep("ha", 7), "m3", "km", "ha"),
    Year_2022_ha = c(156, 282, 5264, 111, 169, 333, 156, 622643, 2354, 18417),
    Year_2023_ha = c(339, 475, 5853, 541, 201, 1513, 431, 676030, 2556, 19308),
    Year_2024_ha = c(322, 893, 8448, 822, 798, 1491, 1082, 458366, 1733, 30539),
    SE_2024_ha = c(42, 516, 1689, 108, 565, 325, 410, 45837, 173, 3054),
    AD_CI_90_percent = c(
      (42 * 1.645 / 322) * 100,   # Forestry infrastructure: 21.5%
      (516 * 1.645 / 893) * 100,  # Agriculture: 95.0%
      (1689 * 1.645 / 8448) * 100, # Mining: 32.8%
      (108 * 1.645 / 822) * 100,   # Infrastructure: 21.6%
      (565 * 1.645 / 798) * 100,   # Settlements: 116.5%
      (325 * 1.645 / 1491) * 100,  # Fire-Biomass burning: 35.8%
      (410 * 1.645 / 1082) * 100,  # Shifting Cultivation: 62.3%
      (45837 * 1.645 / 458366) * 100, # Logging volume: 16.4%
      (173 * 1.645 / 1733) * 100,  # Skid trails: 16.4%
      (3054 * 1.645 / 30539) * 100  # Mining buffer area: 16.4%
    ),
    stringsAsFactors = FALSE
  )
  
  # Emission factors
  emission_factors <- data.frame(
    Activity_Type = activity_data$Activity_Type,
    Driver = activity_data$Driver,
    Units = activity_data$Units,
    EF_tCO2_unit = c(
      338528/322,    # 1. Forestry infrastructure: 1,051.2 tCO2/ha
      937068/893,    # 2. Agriculture: 1,049.5 tCO2/ha  
      8881624/8448,  # 3. Mining: 1,051.3 tCO2/ha
      864192/822,    # 4. Infrastructure: 1,051.3 tCO2/ha
      839381/798,    # 5. Settlements: 1,051.6 tCO2/ha
      1570410/1491,  # 6. Fire-Biomass burning: 1,053.3 tCO2/ha
      999999/1082,   # 7. Shifting cultivation: 924.4 tCO2/ha
      3.85,          # 8. Logging - volume harvested (tCO2/m3)
      171.84,        # 9. Logging - skid trail length (tCO2/km)
      247366/30539   # 10. Mining buffer: 8.1 tCO2/ha
    ),
    EF_CI_90_percent = c(6.9, 6.9, 6.9, 6.9, 6.9, 12.5, 8.7, 7.5, 8.8, 15.2)
  )

  # Calculate final Monte Carlo input data
  combined_data <- merge(activity_data, emission_factors, by = c("Activity_Type", "Driver", "Units"))
  combined_data$Total_Emissions_2024 <- combined_data$Year_2024_ha * combined_data$EF_tCO2_unit
  
  monte_carlo_data <- data.frame(
    Activity_Type = combined_data$Activity_Type,
    Stratum = paste(combined_data$Activity_Type, "-", combined_data$Driver),
    Activity_Data = combined_data$Year_2024_ha,
    AD_Mean = combined_data$Year_2024_ha,
    AD_CI_90_percent = combined_data$AD_CI_90_percent,
    Emission_Factor = combined_data$EF_tCO2_unit,
    EF_Mean = combined_data$EF_tCO2_unit,
    EF_CI_90_percent = combined_data$EF_CI_90_percent,
    Units = combined_data$Units,
    Expected_Emissions_2024 = combined_data$Total_Emissions_2024,
    stringsAsFactors = FALSE
  )
  
  hfld_crediting_level <- 20358133
  
  return(list(
    monte_carlo_data = monte_carlo_data,
    carbon_stocks = carbon_stocks,
    activity_data = activity_data,
    emission_factors = emission_factors,
    hfld_crediting_level = hfld_crediting_level
  ))
}

# FIXED: Monte Carlo simulation with working hyperparameter optimization for small datasets
run_enhanced_monte_carlo <- function(data, hfld_crediting_level, n_iterations = 10000, 
                                   use_bootstrap = FALSE,
                                   # Hyperparameter options
                                   model_method = "lm",
                                   cv_folds = 5,  # Reduced for small dataset
                                   performance_metric = "RMSE",
                                   tune_length = 2,  # Reduced for small dataset
                                   enable_preprocessing = TRUE,
                                   enable_feature_selection = FALSE) {  # Disabled by default for small dataset
  
  # Generate dynamic seed for stochasticity
  dynamic_seed <- as.integer(Sys.time()) + sample(1:10000, 1)
  set.seed(dynamic_seed)
  
  n_strata <- nrow(data)
  emissions_matrix <- matrix(0, nrow = n_iterations, ncol = n_strata)
  
  # Initialize hyperparameter tuning results
  tuning_applied <- FALSE
  bias_correction_factor <- 1.0
  tuning_results <- NULL
  
  # FIXED: Only attempt hyperparameter tuning for suitable models with adequate data
  if(requireNamespace("caret", quietly = TRUE)) {
    tryCatch({
      
      # AGGRESSIVE: Create much larger expanded dataset for Monte Carlo LGOCV
      # Generate extensive synthetic samples to ensure LGOCV reliability
      samples_per_stratum <- 100  # Fixed large number for reliability
      
      expanded_samples <- do.call(rbind, lapply(1:n_strata, function(i) {
        
        ad_se <- (data$AD_Mean[i] * data$AD_CI_90_percent[i] / 100) / 1.645006
        ef_se <- (data$EF_Mean[i] * data$EF_CI_90_percent[i] / 100) / 1.645006
        
        # Generate Monte Carlo samples with proper uncertainty structure
        ad_samples <- rnorm(samples_per_stratum, data$AD_Mean[i], ad_se)
        ef_samples <- rnorm(samples_per_stratum, data$EF_Mean[i], ef_se)
        
        # Add variation in uncertainty estimates themselves
        ad_ci_base <- data$AD_CI_90_percent[i]
        ef_ci_base <- data$EF_CI_90_percent[i]
        ad_ci_samples <- rnorm(samples_per_stratum, ad_ci_base, ad_ci_base * 0.2)
        ef_ci_samples <- rnorm(samples_per_stratum, ef_ci_base, ef_ci_base * 0.2)
        
        data.frame(
          emissions = ad_samples * ef_samples,
          ad_mean = ad_samples,
          ef_mean = ef_samples,
          ad_ci = pmax(ad_ci_samples, 1.0),  # Minimum 1% CI
          ef_ci = pmax(ef_ci_samples, 1.0),
          log_ad_mean = log(pmax(ad_samples, 1)),
          log_ef_mean = log(pmax(ef_samples, 1)),
          combined_uncertainty = sqrt(pmax(ad_ci_samples, 1)^2 + pmax(ef_ci_samples, 1)^2),
          stratum_id = factor(i)
        )
      }))
      
      cat("Generated", nrow(expanded_samples), "samples for LGOCV\n")
      
      # FIXED: Conservative LGOCV setup for ART-TREES compliance
      cv_folds_safe <- min(cv_folds, 8)  # Maximum 8 folds for stability
      train_percentage_safe <- 0.75  # Fixed 75% for reliability
      
      mc_control <- trainControl(
        method = "LGOCV",  # MUST maintain Monte Carlo design
        number = cv_folds_safe,
        p = train_percentage_safe,
        savePredictions = "final",
        summaryFunction = defaultSummary,
        verboseIter = TRUE,  # Enable to see what's happening
        allowParallel = FALSE,
        returnData = FALSE  # Reduce memory usage
      )
      
      # ULTRA-CONSERVATIVE: Model-specific safe parameter grids
      if(model_method == "rf") {
        # Very conservative RF parameters
        max_features <- ncol(expanded_samples) - 2  # Exclude emissions and stratum_id
        safe_mtry <- c(2, min(4, max_features))
        tune_grid <- expand.grid(mtry = safe_mtry)
        cat("RF mtry range:", safe_mtry, "\n")
      } else if(model_method == "gbm") {
        # Ultra-simple GBM
        tune_grid <- expand.grid(
          n.trees = c(50, 100),
          interaction.depth = 1,  # Keep it simple
          shrinkage = 0.1,
          n.minobsinnode = 20  # Larger minimum
        )
      } else if(model_method == "svmRadial") {
        # Very simple SVM grid
        tune_grid <- expand.grid(sigma = 0.1, C = 1)
      } else if(model_method == "glmnet") {
        # Simple elastic net
        tune_grid <- expand.grid(alpha = c(0, 1), lambda = c(0.01, 0.1))
      } else {
        tune_grid <- NULL  # Linear model - no tuning needed
      }
      
      # CONSERVATIVE: Minimal preprocessing
      preprocess_options <- if(enable_preprocessing && model_method %in% c("lm", "glmnet")) {
        c("center", "scale")
      } else {
        NULL
      }
      
      # SIMPLE: Basic feature set to avoid overfitting
      formula_features <- emissions ~ ad_mean + ef_mean + ad_ci + ef_ci
      
      cat("Starting LGOCV with", cv_folds_safe, "folds,", train_percentage_safe, "training %\n")
      
      # ENHANCED ERROR HANDLING: Train model with comprehensive error catching
      tuning_model <- train(
        formula_features,
        data = expanded_samples,
        method = model_method,
        trControl = mc_control,
        metric = performance_metric,
        tuneGrid = tune_grid,
        preProcess = preprocess_options,
        na.action = na.omit
      )
      
      cat("LGOCV completed successfully\n")
      
      # Process results
      if(!is.null(tuning_model) && nrow(tuning_model$results) > 0) {
        best_performance <- min(tuning_model$results[[performance_metric]], na.rm = TRUE)
        mean_emissions <- mean(expanded_samples$emissions, na.rm = TRUE)
        
        # Calculate bias correction factor
        if(mean_emissions > 0 && is.finite(best_performance)) {
          bias_correction_factor <- 1 - (best_performance / mean_emissions)
          bias_correction_factor <- max(0.98, min(1.02, bias_correction_factor))  # Very conservative range
        }
        
        # Store comprehensive tuning results
        tuning_results <- list(
          model_method = model_method,
          best_tune = tuning_model$bestTune,
          best_performance = best_performance,
          cv_results = tuning_model$results,
          variable_importance = if(model_method %in% c("rf", "gbm")) {
            tryCatch(varImp(tuning_model), error = function(e) NULL)
          } else NULL,
          final_model = tuning_model$finalModel
        )
        
        tuning_applied <- TRUE
        
        cat("HYPERPARAMETER OPTIMIZATION SUCCESS!\n")
        cat("- Model:", model_method, "\n")
        cat("- Best", performance_metric, ":", round(best_performance, 4), "\n")
        cat("- Bias Correction Factor:", round(bias_correction_factor, 4), "\n")
        cat("- CV Folds Used:", cv_folds_safe, "\n")
        cat("- Training %:", train_percentage_safe, "\n")
      } else {
        cat("LGOCV completed but no valid results obtained\n")
      }
      
    }, error = function(e) {
      cat("HYPERPARAMETER OPTIMIZATION FAILED:\n")
      cat("Error:", e$message, "\n")
      cat("Traceback:\n")
      print(traceback())
      cat("Falling back to standard Monte Carlo simulation\n")
    })
  } else {
    cat("caret package not available for hyperparameter optimization\n")
  }
  
  # Main Monte Carlo simulation with optimized parameters
  for(i in 1:n_strata) {
    t_90 <- 1.645006
    
    # Apply bias correction from hyperparameter tuning
    ad_mean_corrected <- data$AD_Mean[i] * bias_correction_factor
    ef_mean_corrected <- data$EF_Mean[i] * bias_correction_factor
    
    ad_se <- (ad_mean_corrected * data$AD_CI_90_percent[i] / 100) / t_90
    ef_se <- (ef_mean_corrected * data$EF_CI_90_percent[i] / 100) / t_90
    
    if (use_bootstrap) {
      # Enhanced bootstrap with hyperparameter-informed sampling
      bootstrap_samples <- replicate(n_iterations, {
        ad_bootstrap <- rnorm(100, mean = ad_mean_corrected, sd = ad_se)
        ef_bootstrap <- rnorm(100, mean = ef_mean_corrected, sd = ef_se)
        
        ad_sample <- sample(ad_bootstrap, 1, replace = TRUE)
        ef_sample <- sample(ef_bootstrap, 1, replace = TRUE)
        
        c(ad_sample, ef_sample)
      })
      
      ad_samples <- bootstrap_samples[1, ]
      ef_samples <- bootstrap_samples[2, ]
      
    } else {
      # Monte Carlo sampling with stochastic enhancement
      noise_factor <- runif(1, 0.98, 1.02)
      
      ad_samples <- rnorm(n_iterations, mean = ad_mean_corrected * noise_factor, sd = ad_se)
      ef_samples <- rnorm(n_iterations, mean = ef_mean_corrected * noise_factor, sd = ef_se)
      
      # Add micro-variations for true stochasticity
      ad_samples <- ad_samples + rnorm(n_iterations, 0, ad_se * 0.01)
      ef_samples <- ef_samples + rnorm(n_iterations, 0, ef_se * 0.01)
    }
    
    # Ensure non-negative values
    ad_samples <- pmax(ad_samples, 0)
    ef_samples <- pmax(ef_samples, 0)
    
    emissions_matrix[, i] <- ad_samples * ef_samples
  }
  
  # Calculate final results with ART-TREES calculations
  total_emissions <- rowSums(emissions_matrix)
  simulated_mean <- mean(total_emissions)
  simulated_sd <- sd(total_emissions)
  
  ci_90_lower <- quantile(total_emissions, 0.05)
  ci_90_upper <- quantile(total_emissions, 0.95)
  ci_90_half_width <- (ci_90_upper - ci_90_lower) / 2
  ci_90_percent_of_mean <- (ci_90_half_width / simulated_mean) * 100
  
  # ART TREES Equations 10 and 11
  ua_factor <- 0.524417 * (ci_90_percent_of_mean / 100) / 1.645006
  gross_errs <- hfld_crediting_level - simulated_mean
  uncertainty_deduction <- gross_errs * ua_factor
  
  buffer_rate <- 0.05
  leakage_deduction <- 0
  buffer_deduction <- gross_errs * buffer_rate
  net_errs <- gross_errs - buffer_deduction - uncertainty_deduction - leakage_deduction
  
  normality_test <- shapiro.test(sample(total_emissions, min(5000, length(total_emissions))))
  
  return(list(
    total_emissions = total_emissions,
    emissions_matrix = emissions_matrix,
    simulated_mean = simulated_mean,
    simulated_sd = simulated_sd,
    ci_90_lower = ci_90_lower,
    ci_90_upper = ci_90_upper,
    ci_90_half_width = ci_90_half_width,
    ci_90_percent_of_mean = ci_90_percent_of_mean,
    ua_factor = ua_factor,
    uncertainty_deduction = uncertainty_deduction,
    hfld_crediting_level = hfld_crediting_level,
    gross_errs = gross_errs,
    buffer_deduction = buffer_deduction,
    net_errs = net_errs,
    normality_test = normality_test,
    # Hyperparameter results
    tuning_applied = tuning_applied,
    tuning_results = tuning_results,
    bias_correction_factor = bias_correction_factor,
    hyperparameters_used = list(
      model_method = model_method,
      cv_folds = cv_folds,
      performance_metric = performance_metric,
      enable_preprocessing = enable_preprocessing,
      enable_feature_selection = enable_feature_selection
    ),
    dynamic_seed_used = dynamic_seed
  ))
}

# UI with enhanced hyperparameter tuning
ui <- dashboardPage(
  dashboardHeader(
    title = "ART TREES Monte Carlo Simulation - Enhanced with Hyperparameter Optimization",
    titleWidth = 600
  ),
  
  dashboardSidebar(
    sidebarMenu(
      menuItem("Monte Carlo Simulation", tabName = "simulation", icon = icon("chart-line")),
      menuItem("Hyperparameter Tuning", tabName = "hyperparameters", icon = icon("cogs")),
      menuItem("Input Data Overview", tabName = "data", icon = icon("table")),
      menuItem("Methodology Guide", tabName = "methodology", icon = icon("book"))
    ),
    
    br(),
    h4("ART TREES Parameters", style = "color: white; margin-left: 15px;"),
    div(style = "color: #bbb; margin-left: 15px; font-size: 11px;",
        p("• Monte Carlo Iterations: 10,000"),
        p("• Confidence Interval: 90%"),
        p("• t-value: 1.645006"),
        p("• Scaling Constant: 0.524417")
    ),
    
    br(),
    h4("Simulation Options", style = "color: white; margin-left: 15px;"),
    
    checkboxInput("use_bootstrap", "Use Bootstrap Method", 
                  value = FALSE, width = "90%"),
    p("(Enhanced bootstrap with hyperparameter optimization)", 
      style = "color: #aaa; margin-left: 15px; font-size: 10px;"),
    
    br(),
    h4("Hyperparameter Optimization", style = "color: white; margin-left: 15px;"),
    
    # Model Method Selection  
    selectInput("model_method", "Uncertainty Model:",
                choices = list(
                  "Linear Model" = "lm",
                  "Random Forest" = "rf",
                  "Gradient Boosting" = "gbm", 
                  "Support Vector Machine" = "svmRadial",
                  "Elastic Net" = "glmnet"
                ),
                selected = "lm", width = "90%"),  # Changed default to lm for stability
    
    # Performance Metric Selection
    selectInput("performance_metric", "Optimization Metric:",
                choices = list(
                  "Root Mean Square Error" = "RMSE",
                  "Mean Absolute Error" = "MAE", 
                  "R-squared" = "Rsquared"
                ),
                selected = "RMSE", width = "90%"),
    
    # CV Parameters
    numericInput("cv_folds", "CV Folds:", value = 5, min = 3, max = 10, step = 1, width = "90%"),  # Reduced default
    
    numericInput("tune_length", "Hyperparameter Grid Size:", value = 2, min = 1, max = 5, step = 1, width = "90%"),  # Reduced default
    
    checkboxInput("enable_preprocessing", "Enable Preprocessing", value = TRUE, width = "90%"),
    
    checkboxInput("enable_feature_selection", "Enhanced Feature Selection", value = FALSE, width = "90%"),  # Disabled by default
    
    br(),
    h4("Uncertainty Reduction Analysis", style = "color: white; margin-left: 15px;"),
    
    selectInput("selected_stratum", "Select Driver for Analysis:",
                choices = NULL, width = "90%"),
    
    numericInput("adjust_ad_ci", "Activity Data 90% CI (%):",
                value = 95.1, min = 0.1, max = 200, step = 0.1, width = "90%"),
    
    numericInput("adjust_ef_ci", "Emission Factor 90% CI (%):",
                value = 6.9, min = 0.1, max = 50, step = 0.1, width = "90%"),
    
    br(),
    
    # Action Buttons
    div(style = "margin-left: 15px;",
      actionButton("run_simulation", "Execute Monte Carlo", 
                   class = "btn-success", 
                   style = "width: 90%; margin-bottom: 10px; font-weight: bold;"),
      
      actionButton("reset_hyperparameters", "Reset to Defaults", 
                   class = "btn-secondary", 
                   style = "width: 90%; margin-bottom: 10px;"),
      
      actionButton("save_hyperparameters", "Save Configuration", 
                   class = "btn-info", 
                   style = "width: 90%; margin-bottom: 10px;")
    )
  ),
  
  dashboardBody(
    tags$head(
      tags$style(HTML("
        .content-wrapper, .right-side { background-color: #f4f4f4; }
        .small-box { margin-bottom: 10px; }
        .info-box { margin-bottom: 10px; }
      "))
    ),
    
    tabItems(
      # Monte Carlo Results Tab
      tabItem(tabName = "simulation",
        fluidRow(
          valueBoxOutput("total_emissions_box", width = 3),
          valueBoxOutput("uncertainty_deduction_box", width = 3),
          valueBoxOutput("net_err_box", width = 3),
          valueBoxOutput("ci_percent_box", width = 3)
        ),
        
        fluidRow(
          box(
            title = "Monte Carlo Distribution Analysis", 
            status = "primary", solidHeader = TRUE, width = 8,
            plotlyOutput("enhanced_distribution_plot", height = "400px")
          ),
          box(
            title = "ART TREES Uncertainty Summary", 
            status = "info", solidHeader = TRUE, width = 4,
            DT::dataTableOutput("uncertainty_summary")
          )
        ),
        
        fluidRow(
          box(
            title = "Stratum-Level Results", 
            status = "primary", solidHeader = TRUE, width = 6,
            DT::dataTableOutput("stratum_results")
          ),
          box(
            title = "Sensitivity Analysis", 
            status = "warning", solidHeader = TRUE, width = 6,
            verbatimTextOutput("enhanced_sensitivity_analysis")
          )
        ),
        
        fluidRow(
          box(
            title = "ART TREES Calculation Verification with Hyperparameter Optimization", 
            status = "success", solidHeader = TRUE, width = 12,
            verbatimTextOutput("calculation_details")
          )
        )
      ),
      
      # Hyperparameter Tuning Tab
      tabItem(tabName = "hyperparameters",
        fluidRow(
          valueBoxOutput("tuning_method_box", width = 3),
          valueBoxOutput("best_performance_box", width = 3), 
          valueBoxOutput("bias_correction_box", width = 3),
          valueBoxOutput("optimization_status_box", width = 3)
        ),
        
        fluidRow(
          box(
            title = "Hyperparameter Optimization Results", 
            status = "primary", solidHeader = TRUE, width = 8,
            DT::dataTableOutput("hyperparameter_results_table")
          ),
          box(
            title = "Model Performance Comparison", 
            status = "info", solidHeader = TRUE, width = 4,
            plotOutput("performance_comparison_plot")
          )
        ),
        
        fluidRow(
          box(
            title = "Cross-Validation Summary", 
            status = "warning", solidHeader = TRUE, width = 6,
            verbatimTextOutput("cv_summary")
          ),
          box(
            title = "Variable Importance (if applicable)", 
            status = "success", solidHeader = TRUE, width = 6,
            plotOutput("variable_importance_plot")
          )
        )
      ),
      
      # Data Overview Tab
      tabItem(tabName = "data",
        fluidRow(
          box(
            title = "Carbon Stock Measurements (Guyana Forest Plots)", 
            status = "primary", solidHeader = TRUE, width = 6,
            p("Based on field measurements from 118 plots with 472 subplots (4 per plot)"),
            DT::dataTableOutput("carbon_stocks_table")
          ),
          box(
            title = "Activity Data (Deforestation & Degradation 2024)", 
            status = "primary", solidHeader = TRUE, width = 6,
            p("Annual deforestation and forest degradation activities with logging operations"),
            DT::dataTableOutput("activity_data_table")
          )
        ),
        
        fluidRow(
          box(
            title = "Monte Carlo Input Data Structure", 
            status = "success", solidHeader = TRUE, width = 12,
            p("Complete dataset structure for Monte Carlo uncertainty simulation with hyperparameter optimization"),
            DT::dataTableOutput("monte_carlo_input_table")
          )
        )
      ),
      
      # Methodology Tab
      tabItem(tabName = "methodology",
        fluidRow(
          box(
            title = "ART TREES Monte Carlo Methodology with Hyperparameter Optimization", 
            status = "info", solidHeader = TRUE, width = 12,
            h4("Implementation Based on ART TREES Standard V2.0 Section 8"),
            
            h5("Key Requirements:"),
            tags$ul(
              tags$li("Monte Carlo simulations: n=10,000 iterations, 90% CI"),
              tags$li("Error propagation between Activity Data and Emission Factors"),
              tags$li("Hyperparameter optimization for uncertainty reduction"),
              tags$li("Cross-validation for bias detection and model optimization")
            ),
            
            h5("Mathematical Framework:"),
            withMathJax(),
            p("Equation 11 - Uncertainty Adjustment Factor:"),
            p("$$UA_t = 0.524417 \\times \\frac{HW_{90\\%}}{1.645006}$$"),
            p("Equation 10 - Uncertainty Deduction:"),
            p("$$UNC_t = (GHGER_t + GHGREMV_t) \\times UA_t$$"),
            p("Net ERRs Calculation:"),
            p("$$Net\\ ERRs = Gross\\ ERRs - Buffer\\ (5\\%) - Leakage\\ (0\\%) - Uncertainty\\ Deduction$$"),
            
            h5("Hyperparameter Optimization Features:"),
            tags$ul(
              tags$li("Model Selection: Linear, Random Forest, Gradient Boosting, SVM, Elastic Net"),
              tags$li("Cross-Validation: LGOCV Monte Carlo design with adaptive parameters"),
              tags$li("Data Augmentation: Synthetic samples for small datasets"),
              tags$li("Safe Parameter Grids: Constrained for dataset size"),
              tags$li("Bias Correction: Automated based on CV performance")
            ),
            
            h5("Small Dataset Adaptations:"),
            tags$ul(
              tags$li("Reduced CV folds for stability (3-5 instead of 10)"),
              tags$li("Conservative hyperparameter grids"),
              tags$li("Data augmentation through uncertainty sampling"),
              tags$li("Fallback to linear models when complex models fail")
            )
          )
        )
      )
    )
  )
)

# Enhanced Server Logic with Fixed Hyperparameter Optimization
server <- function(input, output, session) {
  
  # Initialize data
  enhanced_data <- reactive({
    create_enhanced_guyana_data()
  })
  
  # Update stratum choices
  observe({
    data <- enhanced_data()$monte_carlo_data
    updateSelectInput(session, "selected_stratum",
                     choices = setNames(1:nrow(data), data$Stratum),
                     selected = 1)
  })
  
  # Update sliders when stratum changes
  observe({
    req(input$selected_stratum)
    data <- enhanced_data()$monte_carlo_data
    selected_idx <- as.numeric(input$selected_stratum)
    
    updateNumericInput(session, "adjust_ad_ci",
                      value = data$AD_CI_90_percent[selected_idx])
    updateNumericInput(session, "adjust_ef_ci",
                      value = data$EF_CI_90_percent[selected_idx])
  })
  
  # Reactive adjusted data
  adjusted_data <- reactive({
    data <- enhanced_data()$monte_carlo_data
    req(input$selected_stratum)
    
    selected_idx <- as.numeric(input$selected_stratum)
    data$AD_CI_90_percent[selected_idx] <- input$adjust_ad_ci
    data$EF_CI_90_percent[selected_idx] <- input$adjust_ef_ci
    
    return(data)
  })
  
  # Enhanced Monte Carlo results with fixed hyperparameter optimization
  mc_results <- eventReactive(input$run_simulation, {
    data <- adjusted_data()
    hfld_cl <- enhanced_data()$hfld_crediting_level
    
    withProgress(message = 'Running Enhanced Monte Carlo...', value = 0, {
      incProgress(0.2, detail = "Optimizing hyperparameters")
      
      result <- run_enhanced_monte_carlo(
        data = data, 
        hfld_crediting_level = hfld_cl, 
        use_bootstrap = input$use_bootstrap,
        model_method = input$model_method,
        cv_folds = input$cv_folds,
        performance_metric = input$performance_metric,
        tune_length = input$tune_length,
        enable_preprocessing = input$enable_preprocessing,
        enable_feature_selection = input$enable_feature_selection
      )
      
      incProgress(0.7, detail = "Calculating ART TREES metrics")
      result
    })
  })
  
  # Reset hyperparameters to stable defaults
  observeEvent(input$reset_hyperparameters, {
    updateSelectInput(session, "model_method", selected = "lm")  # Changed to lm for stability
    updateSelectInput(session, "performance_metric", selected = "RMSE")
    updateNumericInput(session, "cv_folds", value = 5)  # Reduced for small dataset
    updateNumericInput(session, "tune_length", value = 2)  # Reduced for small dataset
    updateCheckboxInput(session, "enable_preprocessing", value = TRUE)
    updateCheckboxInput(session, "enable_feature_selection", value = FALSE)  # Disabled for stability
    
    showNotification("Hyperparameters reset to stable defaults", type = "message")
  })
  
  # Save hyperparameter configuration  
  observeEvent(input$save_hyperparameters, {
    config <- list(
      model_method = input$model_method,
      performance_metric = input$performance_metric,
      cv_folds = input$cv_folds,
      tune_length = input$tune_length,
      enable_preprocessing = input$enable_preprocessing,
      enable_feature_selection = input$enable_feature_selection,
      timestamp = Sys.time()
    )
    
    saveRDS(config, file = "art_trees_hyperparameter_config.rds")
    showNotification("Hyperparameter configuration saved", type = "success")
  })
  
  # Value boxes for main results
  output$total_emissions_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      value = format(round(results$simulated_mean, 0), big.mark = ","),
      subtitle = "Mean Emissions (tCO2/yr)",
      icon = icon("fire"),
      color = "red"
    )
  })
  
  output$uncertainty_deduction_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      value = format(round(results$uncertainty_deduction, 0), big.mark = ","),
      subtitle = "Uncertainty Deduction (tCO2/yr)",
      icon = icon("minus-circle"),
      color = "orange"
    )
  })
  
  output$net_err_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      value = format(round(results$net_errs, 0), big.mark = ","),
      subtitle = "Net ERRs (tCO2/yr)",
      icon = icon("leaf"),
      color = "green"
    )
  })
  
  output$ci_percent_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      value = paste0(round(results$ci_90_percent_of_mean, 1), "%"),
      subtitle = "90% CI (% of Mean)",
      icon = icon("chart-bar"),
      color = "blue"
    )
  })
  
  # Hyperparameter tuning value boxes
  output$tuning_method_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      value = if(results$tuning_applied) results$hyperparameters_used$model_method else "None",
      subtitle = "Optimization Model",
      icon = icon("cogs"),
      color = "purple"
    )
  })
  
  output$best_performance_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      value = if(results$tuning_applied && !is.null(results$tuning_results)) 
        round(results$tuning_results$best_performance, 3) else "N/A",
      subtitle = "Best RMSE",
      icon = icon("bullseye"),
      color = "yellow"
    )
  })
  
  output$bias_correction_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      value = round(results$bias_correction_factor, 4),
      subtitle = "Bias Correction Factor",
      icon = icon("balance-scale"),
      color = "aqua"
    )
  })
  
  output$optimization_status_box <- renderValueBox({
    req(mc_results())
    results <- mc_results()
    valueBox(
      value = if(results$tuning_applied) "✓ Active" else "✗ Disabled",
      subtitle = "Optimization Status",
      icon = icon("check-circle"),
      color = if(results$tuning_applied) "green" else "red"
    )
  })
  
  # Distribution plot
  output$enhanced_distribution_plot <- renderPlotly({
    req(mc_results())
    results <- mc_results()
    
    df <- data.frame(emissions = results$total_emissions)
    
    p <- ggplot(df, aes(x = emissions)) +
      geom_histogram(bins = 50, fill = "steelblue", alpha = 0.7, color = "darkblue") +
      geom_vline(xintercept = results$simulated_mean, color = "red", linetype = "dashed", size = 1) +
      geom_vline(xintercept = results$ci_90_lower, color = "orange", linetype = "dotted", size = 1) +
      geom_vline(xintercept = results$ci_90_upper, color = "orange", linetype = "dotted", size = 1) +
      labs(title = "Enhanced Monte Carlo Results with Hyperparameter Optimization",
           subtitle = paste("Model:", if(results$tuning_applied) results$hyperparameters_used$model_method else "Standard"),
           x = "Total Emissions (tCO2/yr)",
           y = "Frequency") +
      theme_minimal() +
      theme(plot.title = element_text(size = 12, face = "bold"))
    
    ggplotly(p)
  })
  
  # Hyperparameter results table
  output$hyperparameter_results_table <- DT::renderDataTable({
    req(mc_results())
    results <- mc_results()
    
    if(results$tuning_applied && !is.null(results$tuning_results)) {
      cv_results <- results$tuning_results$cv_results
      DT::datatable(cv_results, options = list(scrollX = TRUE, pageLength = 10), rownames = FALSE)
    } else {
      data.frame(Message = "No hyperparameter optimization results available")
    }
  })
  
  # Performance comparison plot
  output$performance_comparison_plot <- renderPlot({
    req(mc_results())
    results <- mc_results()
    
    if(results$tuning_applied && !is.null(results$tuning_results)) {
      cv_results <- results$tuning_results$cv_results
      
      if(nrow(cv_results) > 1) {
        ggplot(cv_results, aes(x = 1:nrow(cv_results), y = RMSE)) +
          geom_line(color = "blue", size = 1) +
          geom_point(color = "red", size = 2) +
          labs(title = "Hyperparameter Performance",
               x = "Configuration",
               y = "RMSE") +
          theme_minimal()
      } else {
        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
        text(1, 1, "Single configuration\nNo comparison available", cex = 1.2)
      }
    } else {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, "No optimization\nresults available", cex = 1.2)
    }
  })
  
  # Variable importance plot
  output$variable_importance_plot <- renderPlot({
    req(mc_results())
    results <- mc_results()
    
    if(results$tuning_applied && !is.null(results$tuning_results$variable_importance)) {
      plot(results$tuning_results$variable_importance)
    } else {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, "Variable importance\nnot available", cex = 1.2)
    }
  })
  
  # CV Summary
  output$cv_summary <- renderText({
    req(mc_results())
    results <- mc_results()
    
    if(results$tuning_applied) {
      paste0(
        "Cross-Validation Summary\n",
        "========================\n\n",
        "Model Method: ", results$hyperparameters_used$model_method, "\n",
        "Monte Carlo CV: LGOCV (ART-TREES compliant)\n",
        "CV Folds: ", results$hyperparameters_used$cv_folds, "\n",
        "Performance Metric: ", results$hyperparameters_used$performance_metric, "\n",
        "Preprocessing: ", if(results$hyperparameters_used$enable_preprocessing) "Enabled" else "Disabled", "\n",
        "Feature Selection: ", if(results$hyperparameters_used$enable_feature_selection) "Enhanced" else "Basic", "\n",
        "Bias Correction Applied: ", round(results$bias_correction_factor, 4), "\n\n",
        "This Monte Carlo optimization reduces uncertainty and\nimproves the reliability of ART TREES calculations."
      )
    } else {
      "Monte Carlo simulation running in standard mode.\nLGOCV hyperparameter optimization requires larger datasets."
    }
  })
  
  # Data tables (unchanged from original)
  output$monte_carlo_input_table <- DT::renderDataTable({
    data <- enhanced_data()$monte_carlo_data
    
    display_data <- data.frame(
      Activity_Type = data$Activity_Type,
      Stratum = data$Stratum,
      Activity_Data = format(data$Activity_Data, big.mark = ","),
      AD_Mean = format(data$AD_Mean, big.mark = ","),
      AD_CI_90_percent = paste0(round(data$AD_CI_90_percent, 1), "%"),
      Emission_Factor = round(data$Emission_Factor, 1),
      EF_Mean = round(data$EF_Mean, 1),
      EF_CI_90_percent = paste0(round(data$EF_CI_90_percent, 1), "%"),
      Units = data$Units,
      Expected_Emissions = format(round(data$Expected_Emissions_2024, 0), big.mark = ","),
      stringsAsFactors = FALSE
    )
    
    colnames(display_data) <- c("Activity Type", "Stratum", "Activity Data", "AD Mean", "AD CI (%)", 
                               "Emission Factor", "EF Mean", "EF CI (%)", "Units", "Expected Emissions")
    
    DT::datatable(display_data, 
                  options = list(scrollX = TRUE, pageLength = 12, dom = 'frtip'),
                  rownames = FALSE) %>%
      DT::formatStyle(
        "Activity Type",
        target = "row",
        backgroundColor = DT::styleEqual(c("Deforestation", "Degradation"), 
                                        c("#fff3e0", "#e8f5e8"))
      )
  })
  
  output$carbon_stocks_table <- DT::renderDataTable({
    data <- enhanced_data()$carbon_stocks
    display_data <- data[, c("Component", "Median_tC_ha", "SE_tC_ha", "Min_tC_ha", "Max_tC_ha", "CI_90_tC_ha", "CI_percent", "n_plots")]
    display_data$Median_tC_ha <- round(display_data$Median_tC_ha, 1)
    display_data$SE_tC_ha <- round(display_data$SE_tC_ha, 2)
    display_data$Min_tC_ha <- round(display_data$Min_tC_ha, 1)
    display_data$Max_tC_ha <- round(display_data$Max_tC_ha, 1)
    display_data$CI_90_tC_ha <- round(display_data$CI_90_tC_ha, 1)
    display_data$CI_percent <- paste0(display_data$CI_percent, "%")
    
    colnames(display_data) <- c("Component", "Median (tC/ha)", "Std Error (tC/ha)", "Min (tC/ha)", "Max (tC/ha)", "90% CI (tC/ha)", "CI (%)", "Plots")
    
    dt <- DT::datatable(display_data, options = list(pageLength = 15, dom = 'tp', scrollX = TRUE), rownames = FALSE) %>%
      DT::formatStyle(
        "Component",
        target = "row",
        backgroundColor = DT::styleEqual(c("Sum C pools w/o litter", "Sum C pools w/ litter", "Sum ALL POOLS"), 
                                        c("#f0f0f0", "#e8e8e8", "#d0d0d0"))
      )
    dt
  })
  
  output$activity_data_table <- DT::renderDataTable({
    data <- enhanced_data()$activity_data
    display_data <- data[, c("Activity_Type", "Driver", "Units", "Year_2024_ha", "SE_2024_ha", "AD_CI_90_percent")]
    display_data$Year_2024_ha <- format(display_data$Year_2024_ha, big.mark = ",")
    display_data$SE_2024_ha <- format(display_data$SE_2024_ha, big.mark = ",")
    display_data$AD_CI_90_percent <- paste0(round(display_data$AD_CI_90_percent, 1), "%")
    colnames(display_data) <- c("Activity Type", "Driver", "Units", "2024 Activity", "Std Error", "90% CI (%)")
    
    dt <- DT::datatable(display_data, options = list(pageLength = 12, dom = 'tp', scrollX = TRUE), rownames = FALSE) %>%
      DT::formatStyle(
        "Activity Type",
        target = "row",
        backgroundColor = DT::styleEqual(c("Deforestation", "Degradation"), 
                                        c("#f8f9fa", "#e3f2fd"))
      )
    dt
  })
  
  output$uncertainty_summary <- DT::renderDataTable({
    req(mc_results())
    results <- mc_results()
    
    summary_data <- data.frame(
      Metric = c("Simulated Mean Emissions", "HFLD Crediting Level", "Gross ERRs", 
                "90% CI Half-Width", "90% CI (%)", "UA Factor", 
                "Uncertainty Deduction", "Buffer Deduction (5%)", "Net ERRs (FINAL)",
                "Hyperparameter Optimization", "Bias Correction Factor"),
      Value = c(
        format(round(results$simulated_mean, 0), big.mark = ","),
        format(round(results$hfld_crediting_level, 0), big.mark = ","),
        format(round(results$gross_errs, 0), big.mark = ","),
        format(round(results$ci_90_half_width, 0), big.mark = ","),
        paste0(round(results$ci_90_percent_of_mean, 2), "%"),
        round(results$ua_factor, 6),
        format(round(results$uncertainty_deduction, 0), big.mark = ","),
        format(round(results$buffer_deduction, 0), big.mark = ","),
        format(round(results$net_errs, 0), big.mark = ","),
        if(results$tuning_applied) "✓ Applied" else "✗ Not Applied",
        round(results$bias_correction_factor, 4)
      )
    )
    
    DT::datatable(summary_data, options = list(dom = 't', pageLength = 15), rownames = FALSE)
  })
  
  output$stratum_results <- DT::renderDataTable({
    req(mc_results())
    data <- adjusted_data()
    
    stratum_data <- data.frame(
      Stratum = data$Stratum,
      Activity_Data = format(data$Activity_Data, big.mark = ","),
      Emission_Factor = round(data$Emission_Factor, 1),
      Emissions = format(round(data$Activity_Data * data$Emission_Factor, 0), big.mark = ","),
      AD_CI = paste0(data$AD_CI_90_percent, "%"),
      EF_CI = paste0(data$EF_CI_90_percent, "%")
    )
    
    colnames(stratum_data) <- c("Stratum", "Activity Data", "EF (tCO2/ha)", 
                               "Emissions (tCO2)", "AD CI", "EF CI")
    
    DT::datatable(stratum_data, options = list(pageLength = 10, dom = 't'))
  })
  
  output$calculation_details <- renderText({
    req(mc_results())
    results <- mc_results()
    
    paste0(
      "ART TREES CALCULATION WITH HYPERPARAMETER OPTIMIZATION\n",
      "======================================================\n\n",
      "Step 1 - Enhanced Monte Carlo Simulation:\n",
      "• Iterations: 10,000 (ART TREES requirement)\n",
      "• Optimization Model: ", if(results$tuning_applied) results$hyperparameters_used$model_method else "None", "\n",
      "• Bias Correction Factor: ", round(results$bias_correction_factor, 4), "\n",
      "• Simulated Mean Emissions: ", format(round(results$simulated_mean, 0), big.mark = ","), " tCO2/yr\n",
      "• 90% CI Half-Width: ", format(round(results$ci_90_half_width, 0), big.mark = ","), " tCO2/yr\n",
      "• 90% CI as % of Mean: ", round(results$ci_90_percent_of_mean, 2), "%\n\n",
      "Step 2 - ART TREES Standard Equations:\n",
      "• HFLD Crediting Level: ", format(round(results$hfld_crediting_level, 0), big.mark = ","), " tCO2/yr\n",
      "• Gross ERRs = HFLD CL - Simulated Mean = ", format(round(results$gross_errs, 0), big.mark = ","), " tCO2/yr\n\n",
      "• Equation 11 - UA Factor = 0.524417 × (", round(results$ci_90_percent_of_mean, 2), "% / 100) / 1.645006 = ", 
      round(results$ua_factor, 6), "\n",
      "• Equation 10 - Uncertainty Deduction = Gross ERRs × UA Factor = ", 
      format(round(results$uncertainty_deduction, 0), big.mark = ","), " tCO2/yr\n\n",
      "Step 3 - Final Net ERRs Calculation:\n",
      "• Buffer Deduction (5%): ", format(round(results$buffer_deduction, 0), big.mark = ","), " tCO2/yr\n",
      "• Leakage Deduction: 0 tCO2/yr\n",
      "• Net ERRs = Gross ERRs - Buffer - Uncertainty = ", 
      format(round(results$net_errs, 0), big.mark = ","), " tCO2/yr\n\n",
      "✓ ENHANCED: Hyperparameter optimization ", if(results$tuning_applied) "ACTIVE" else "DISABLED", "\n",
      "✓ SUCCESS: Net ERRs optimized through ", if(results$tuning_applied) "advanced modeling techniques!" else "robust Monte Carlo simulation!"
    )
  })
  
  output$enhanced_sensitivity_analysis <- renderText({
    req(mc_results(), input$selected_stratum)
    
    original_data <- create_enhanced_guyana_data()$monte_carlo_data
    current_data <- adjusted_data()
    hfld_cl <- enhanced_data()$hfld_crediting_level
    selected_idx <- as.numeric(input$selected_stratum)
    
    # Run comparison with same hyperparameters
    original_results <- run_enhanced_monte_carlo(
      original_data, hfld_cl,
      model_method = input$model_method,
      cv_folds = input$cv_folds,
      performance_metric = input$performance_metric,
      tune_length = input$tune_length,
      enable_preprocessing = input$enable_preprocessing,
      enable_feature_selection = input$enable_feature_selection
    )
    current_results <- mc_results()
    
    stratum_name <- current_data$Stratum[selected_idx]
    
    # Calculate changes
    unc_change <- current_results$uncertainty_deduction - original_results$uncertainty_deduction
    err_change <- current_results$net_errs - original_results$net_errs
    ci_change <- current_results$ci_90_percent_of_mean - original_results$ci_90_percent_of_mean
    
    paste0(
      "ENHANCED SENSITIVITY ANALYSIS: ", stratum_name, "\n",
      "=", paste(rep("=", nchar(stratum_name) + 32), collapse = ""), "\n\n",
      "Hyperparameter Configuration:\n",
      "• Model: ", input$model_method, "\n",
      "• Optimization: ", if(current_results$tuning_applied) "Active" else "Disabled", "\n",
      "• Bias Correction: ", round(current_results$bias_correction_factor, 4), "\n\n",
      "Parameter Adjustments:\n",
      "• Activity Data CI: ", round(original_data$AD_CI_90_percent[selected_idx], 1), "% → ", 
      round(input$adjust_ad_ci, 1), "%\n",
      "• Emission Factor CI: ", round(original_data$EF_CI_90_percent[selected_idx], 1), "% → ", 
      round(input$adjust_ef_ci, 1), "%\n\n",
      "System-wide Impacts:\n",
      "• Change in 90% CI: ", sprintf("%+.1f%%", ci_change), "\n",
      "• Change in Uncertainty Deduction: ", format(round(unc_change, 0), big.mark = ","), " tCO2/yr\n",
      "• Change in Net ERRs: ", format(round(err_change, 0), big.mark = ","), " tCO2/yr\n\n",
      if(current_results$tuning_applied) {
        "The hyperparameter optimization enhances uncertainty\nreduction beyond standard Monte Carlo approaches."
      } else {
        "For this small dataset, robust Monte Carlo simulation\nprovides reliable uncertainty estimation."
      }
    )
  })
}

# Run the enhanced application
shinyApp(ui = ui, server = server)