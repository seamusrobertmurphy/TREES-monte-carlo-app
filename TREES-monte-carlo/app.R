# ------------------------------------------------------------------------ #
# ART TREES Monte Carlo Simulation
# Government of Guyana
# Compliant with ART TREES Standard Version 2.0, Section 8
# ------------------------------------------------------------------------ #

pacman::p_load(
	BiocManager,
  caret,
	dplyr, DT,
  ggplot2,
  MASS,
	plotly,
  shiny, shinydashboard,shinyjs
	)  

# Guyana data structure with ALL components
create_enhanced_guyana_data <- function() {
  # Carbon stock data using medians and standard errors with min/max values
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

# Monte Carlo simulation with working hyperparameter optimization for small datasets
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
  
  # Only attempt hyperparameter tuning for suitable models with adequate data
  if(requireNamespace("caret", quietly = TRUE)) {
    tryCatch({
      
      # Generate extensive synthetic samples to ensure LGOCV reliability
      samples_per_stratum <- 100
      
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
      
      # Conservative LGOCV setup for ART-TREES compliance
      cv_folds_safe <- min(cv_folds, 8)  # Maximum 8 folds for stability
      train_percentage_safe <- 0.75  # 75% found most reliable
      
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

# ============================================================================
# NEW DISTRIBUTION ANALYSIS FUNCTIONS
# ============================================================================

# Load pre-existing CSV data files
load_preloaded_data <- function() {
  # Initialize list to store data
  preloaded_data <- list()
  
  # Try to load each CSV file
  tryCatch({
    # Load LIF data
    if(file.exists("LIF_Prepared.csv")) {
      lif_data <- read.csv("LIF_Prepared.csv", stringsAsFactors = FALSE)
      # Clean column names - remove spaces and special characters
      names(lif_data) <- gsub("\\s+", "_", names(lif_data))
      names(lif_data) <- gsub("[()]", "", names(lif_data))
      preloaded_data$LIF <- lif_data
    }
    
    # Load AllBiomassData
    if(file.exists("AllBiomassData_Prepared.csv")) {
      biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
      # Clean column names
      names(biomass_data) <- gsub("\\s+", "_", names(biomass_data))
      names(biomass_data) <- gsub("[()]", "", names(biomass_data))
      preloaded_data$Biomass <- biomass_data
    }
    
    # Load LDF data
    if(file.exists("LDF_Prepared.csv")) {
      ldf_data <- read.csv("LDF_Prepared.csv", stringsAsFactors = FALSE)
      # Clean column names
      names(ldf_data) <- gsub("\\s+", "_", names(ldf_data))
      names(ldf_data) <- gsub("[()]", "", names(ldf_data))
      preloaded_data$LDF <- ldf_data
    }
  }, error = function(e) {
    cat("Error loading CSV files:", e$message, "\n")
  })
  
  return(preloaded_data)
}

# Extract data from preloaded datasets
extract_from_preloaded <- function(preloaded_data, variable_name) {
  if(is.null(preloaded_data) || length(preloaded_data) == 0) {
    return(NULL)
  }
  
  # Map variable names to CSV columns
  column_mapping <- list(
    "LIF" = list(file = "LIF", column = "LIF_tC/km"),
    "AGB" = list(file = "Biomass", column = "AGB_tC/ha"),
    "BGB" = list(file = "Biomass", column = "BGB_tC/ha"),
    "Saplings" = list(file = "Biomass", column = "Saplings_tC/ha"),
    "Litter" = list(file = "Biomass", column = "Litter_tC/ha"),
    "Standing_Dead" = list(file = "Biomass", column = "Standing_Dead_tC/ha"),
    "Lying_Dead" = list(file = "Biomass", column = "Lying_Dead_tC/ha"),
    "Soil" = list(file = "Biomass", column = "Soil_tC/ha"),
    "LDF" = list(file = "LDF", column = "LDF")
  )
  
  if(variable_name %in% names(column_mapping)) {
    mapping <- column_mapping[[variable_name]]
    
    if(!is.null(preloaded_data[[mapping$file]])) {
      data <- preloaded_data[[mapping$file]]
      
      # Try to find the column (with flexibility for naming variations)
      possible_columns <- c(
        mapping$column,
        gsub("_", " ", mapping$column),
        gsub("tC/ha", "(tC/ha)", mapping$column),
        gsub("tC/km", "(tC/km)", mapping$column)
      )
      
      for(col in possible_columns) {
        if(col %in% names(data)) {
          return(as.numeric(data[, col]))
        }
      }
      
      # If exact match not found, try partial match
      col_idx <- grep(strsplit(mapping$column, "_")[[1]][1], names(data), ignore.case = TRUE)
      if(length(col_idx) > 0) {
        return(as.numeric(data[, col_idx[1]]))
      }
    }
  }
  
  return(NULL)
}

# Core distribution analysis function
perform_distribution_analysis <- function(data_vector, variable_name, units) {
  # Remove NA values
  data_clean <- na.omit(data_vector)
  n <- length(data_clean)
  
  if(n < 3) {
    return(list(error = "Insufficient data for analysis (n < 3)"))
  }
  
  # Descriptive statistics
  desc_stats <- list(
    n = n,
    mean = mean(data_clean),
    median = median(data_clean),
    sd = sd(data_clean),
    min = min(data_clean),
    max = max(data_clean),
    q25 = quantile(data_clean, 0.25),
    q75 = quantile(data_clean, 0.75),
    iqr = IQR(data_clean),
    cv = if(mean(data_clean) != 0) sd(data_clean) / abs(mean(data_clean)) * 100 else NA
  )
  
  # Simple skewness calculation (avoid external packages)
  skewness <- function(x) {
    n <- length(x)
    if(n < 3) return(NA)
    m <- mean(x)
    s <- sd(x)
    if(s == 0) return(0)
    (sum((x - m)^3) / n) / (s^3)
  }
  
  # Simple kurtosis calculation
  kurtosis <- function(x) {
    n <- length(x)
    if(n < 4) return(NA)
    m <- mean(x)
    s <- sd(x)
    if(s == 0) return(0)
    (sum((x - m)^4) / n) / (s^4) - 3
  }
  
  desc_stats$skewness <- skewness(data_clean)
  desc_stats$kurtosis <- kurtosis(data_clean)
  
  # Shapiro-Wilk test (limit to 5000 samples for performance)
  test_sample <- if(n > 5000) sample(data_clean, 5000) else data_clean
  shapiro_test <- if(n >= 3 && n <= 5000) {
    shapiro.test(test_sample)
  } else {
    list(statistic = NA, p.value = NA)
  }
  
  # Determine normality interpretation
  normality_interpretation <- if(!is.na(shapiro_test$p.value)) {
    if(shapiro_test$p.value > 0.05) {
      "Data appears normally distributed (p > 0.05)"
    } else if(shapiro_test$p.value > 0.01) {
      "Mild deviation from normality (0.01 < p < 0.05)"
    } else {
      "Significant deviation from normality (p < 0.01)"
    }
  } else {
    "Unable to test normality"
  }
  
  # Transformation recommendation based on skewness
  transformation_rec <- if(!is.na(desc_stats$skewness)) {
    if(abs(desc_stats$skewness) < 0.5) {
      "No transformation needed"
    } else if(desc_stats$skewness > 0.5) {
      "Consider log or square root transformation (right-skewed)"
    } else {
      "Consider square or cubic transformation (left-skewed)"
    }
  } else {
    "Unable to assess skewness"
  }
  
  return(list(
    variable_name = variable_name,
    units = units,
    desc_stats = desc_stats,
    shapiro_test = shapiro_test,
    normality_interpretation = normality_interpretation,
    transformation_rec = transformation_rec,
    data = data_clean
  ))
}

# Create distribution visualization
create_distribution_plot <- function(analysis_result) {
  data_df <- data.frame(value = analysis_result$data)
  
  # Create the plot
  p <- ggplot(data_df, aes(x = value)) +
    geom_histogram(aes(y = ..density..), bins = 30, 
                   fill = "steelblue", alpha = 0.7, color = "darkblue") +
    geom_density(color = "red", size = 1.2) +
    geom_vline(xintercept = analysis_result$desc_stats$mean, 
               color = "darkgreen", linetype = "dashed", size = 1) +
    geom_vline(xintercept = analysis_result$desc_stats$median, 
               color = "orange", linetype = "dotted", size = 1) +
    labs(
      title = paste("Distribution Analysis:", analysis_result$variable_name),
      subtitle = paste("Shapiro-Wilk p-value:", 
                      format(analysis_result$shapiro_test$p.value, digits = 4),
                      " | Skewness:", round(analysis_result$desc_stats$skewness, 2),
                      " | n =", analysis_result$desc_stats$n),
      x = paste(analysis_result$variable_name, "(", analysis_result$units, ")"),
      y = "Density"
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      plot.subtitle = element_text(size = 11)
    )
  
  return(p)
}

# Generate synthetic variable data for demonstration
generate_synthetic_variable_data <- function(variable_name, units = "tC/ha") {
  set.seed(123)  # For reproducibility
  
  n <- 500  # Sample size
  
  data <- switch(variable_name,
    "LIF" = rlnorm(n, meanlog = 3, sdlog = 0.5),  # Log-normal for LIF
    "AGB" = rgamma(n, shape = 20, rate = 0.1),    # Gamma for AGB
    "BGB" = rgamma(n, shape = 10, rate = 0.2),    # Gamma for BGB
    "Saplings" = rgamma(n, shape = 2, rate = 0.5), # Low values
    "Litter" = rnorm(n, mean = 3.3, sd = 0.8),    # Normal distribution
    "Standing_Dead" = c(rep(0, n*0.3), rgamma(n*0.7, shape = 1, rate = 0.3)), # Many zeros
    "Lying_Dead" = c(rep(0, n*0.2), rgamma(n*0.8, shape = 2, rate = 0.2)),   # Some zeros
    "Soil" = rnorm(n, mean = 60, sd = 20) + abs(rnorm(n, 0, 10)), # Variable with outliers
    "LDF" = rnorm(n, mean = 8, sd = 2),           # Normal distribution
    rnorm(n, mean = 100, sd = 30)  # Default
  )
  
  # Ensure positive values
  data <- abs(data)
  
  return(data)
}

# Extract variable from CSV based on column mapping
extract_variable_from_csv <- function(csv_data, variable_name, file_path) {
  # Try to find the column in various formats
  possible_names <- c(
    variable_name,
    tolower(variable_name),
    toupper(variable_name),
    gsub("_", " ", variable_name),
    gsub("_", "", variable_name)
  )
  
  # Find matching column
  col_idx <- which(tolower(names(csv_data)) %in% tolower(possible_names))
  
  if(length(col_idx) > 0) {
    return(as.numeric(csv_data[, col_idx[1]]))
  } else {
    # Try to find partial matches
    col_idx <- grep(variable_name, names(csv_data), ignore.case = TRUE)
    if(length(col_idx) > 0) {
      return(as.numeric(csv_data[, col_idx[1]]))
    }
  }
  
  return(NULL)
}

# Read and analyze CSV data
read_and_analyze_csv_data <- function(file_path, variable_name, units) {
  tryCatch({
    # Read CSV
    csv_data <- read.csv(file_path, stringsAsFactors = FALSE)
    
    # Extract variable data
    var_data <- extract_variable_from_csv(csv_data, variable_name, file_path)
    
    if(!is.null(var_data)) {
      # Perform analysis
      analysis <- perform_distribution_analysis(var_data, variable_name, units)
      return(analysis)
    } else {
      return(list(error = paste("Variable", variable_name, "not found in CSV")))
    }
  }, error = function(e) {
    return(list(error = paste("Error reading CSV:", e$message)))
  })
}

# Generate hyperparameter recommendations based on distribution
generate_hyperparameter_recommendations <- function(analysis_result) {
  if(is.null(analysis_result$shapiro_test)) {
    return("No analysis available")
  }
  
  p_value <- analysis_result$shapiro_test$p.value
  skewness <- analysis_result$desc_stats$skewness
  kurtosis <- analysis_result$desc_stats$kurtosis
  cv <- analysis_result$desc_stats$cv
  
  recommendations <- character()
  
  # Model selection based on distribution
  if(p_value > 0.05) {
    recommendations <- c(recommendations, 
      "• Linear models (lm) are appropriate for normally distributed data",
      "• Consider elastic net (glmnet) for regularization if needed")
  } else {
    recommendations <- c(recommendations,
      "• Non-linear models recommended due to non-normal distribution",
      "• Random Forest (rf) can handle skewed distributions well",
      "• Gradient Boosting (gbm) effective for complex patterns")
  }
  
  # Transformation recommendations
  if(abs(skewness) > 1) {
    recommendations <- c(recommendations,
      paste("• High skewness (", round(skewness, 2), ") - consider data transformation"))
  }
  
  # CV strategy based on variability
  if(cv > 50) {
    recommendations <- c(recommendations,
      "• High variability - use more CV folds (8-10) for robust estimates",
      "• Enable preprocessing to normalize scales")
  } else {
    recommendations <- c(recommendations,
      "• Moderate variability - standard CV folds (5) should suffice")
  }
  
  # Feature engineering
  if(kurtosis > 3 || kurtosis < -3) {
    recommendations <- c(recommendations,
      "• Extreme kurtosis detected - consider feature engineering",
      "• Add polynomial features or interaction terms")
  }
  
  return(paste(recommendations, collapse = "\n"))
}

# ============================================================================
# ENHANCED UI WITH DISTRIBUTION ANALYSIS TAB
# ============================================================================

ui <- dashboardPage(
  dashboardHeader(
    title = "ART-TREES Monte Carlo Tool: Normality + Hyperparamter Tuning Enhancements",
    titleWidth = 900,
    tags$li(class = "dropdown",
      tags$a(
        img(src = "https://winrock.org/wp-content/uploads/2021/12/Winrock-logo-R.png", height = "40px", style = "padding: 5px;"),
        href = "https://winrock.org/wp-content/uploads/2021/12/Winrock-logo-R.png",
        target = "_blank"
      )
    )
  ),
  
  dashboardSidebar(
    sidebarMenu(
    	id = "tabs",  
      menuItem("Monte Carlo Simulation", tabName = "simulation", icon = icon("chart-line")),
      menuItem("Input Data Overview", tabName = "data", icon = icon("table")),
      menuItem("Distribution Analysis", tabName = "distribution", icon = icon("chart-bar")),
      menuItem("Hyperparameter Tuning", tabName = "hyperparameters", icon = icon("cogs")),
      menuItem("Methodology Guide", tabName = "methodology", icon = icon("book")),
    	menuItem("Source Code", tabName = "sourcecode", icon = icon("code"))
    ),
    
    br(),
    
    # Distribution Analysis Controls (NEW)
    conditionalPanel(
      condition = "input.tabs == 'distribution'",
      h4("Distribution Analysis", style = "color: white; margin-left: 15px;"),
      
      selectInput("analysis_variable", "Select Variable:",
                  choices = list(
                    "Logging Impact Factor" = "LIF",
                    "Above Ground Biomass" = "AGB",
                    "Below Ground Biomass" = "BGB",
                    "Saplings Carbon" = "Saplings",
                    "Litter Carbon" = "Litter",
                    "Standing Dead Wood" = "Standing_Dead",
                    "Lying Dead Wood" = "Lying_Dead",
                    "Soil Carbon" = "Soil"
                  ),
                  selected = "AGB", width = "90%"),
      
      checkboxInput("use_preloaded_data", "Use Pre-loaded Data",
                    value = TRUE, width = "90%"),
      
      fileInput(
      	"csv_upload", 
      	"Upload CSV File:",
      	accept = c("text/csv", ".csv"),
      	width = "90%"),
      
      actionButton("analyze_distribution", "Analyze Distribution",
                   class = "btn-primary",
                   style = "width: 90%; margin-bottom: 10px;"),
      
      br()
    ),
    
    # ART TREES Parameters
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
                selected = "lm", width = "90%"),
    
    # Performance Metric Selection
    selectInput("performance_metric", "Optimization Metric:",
                choices = list(
                  "Root Mean Square Error" = "RMSE",
                  "Mean Absolute Error" = "MAE", 
                  "R-squared" = "Rsquared"
                ),
                selected = "RMSE", width = "90%"),
    
    # CV Parameters
    numericInput("cv_folds", "CV Folds:", value = 5, min = 3, max = 10, step = 1, width = "90%"),
    
    numericInput("tune_length", "Hyperparameter Grid Size:", value = 2, min = 1, max = 5, step = 1, width = "90%"),
    
    checkboxInput("enable_preprocessing", "Enable Preprocessing", value = TRUE, width = "90%"),
    
    checkboxInput("enable_feature_selection", "Enhanced Feature Selection", value = FALSE, width = "90%"),
    
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
      
      # NEW: Distribution Analysis Tab
      tabItem(tabName = "distribution",
        # Value boxes for distribution status
        fluidRow(
          valueBoxOutput("normality_status_box", width = 3),
          valueBoxOutput("shapiro_p_value_box", width = 3),
          valueBoxOutput("sample_size_box", width = 3),
          valueBoxOutput("transformation_rec_box", width = 3)
        ),
        
        # Main distribution plot and summary
        fluidRow(
          box(
            title = "Variable Distribution Analysis",
            status = "primary", solidHeader = TRUE, width = 8,
            plotOutput("distribution_plot", height = "400px")
          ),
          box(
            title = "Normality Test Summary",
            status = "info", solidHeader = TRUE, width = 4,
            verbatimTextOutput("normality_summary")
          )
        ),
        
        # Statistics and recommendations
        fluidRow(
          box(
            title = "Descriptive Statistics",
            status = "warning", solidHeader = TRUE, width = 6,
            DT::dataTableOutput("descriptive_stats_table")
          ),
          box(
            title = "Hyperparameter Recommendations",
            status = "success", solidHeader = TRUE, width = 6,
            verbatimTextOutput("hyperparameter_recommendations")
          )
        ),
        
        # All variables summary
        fluidRow(
          box(
            title = "Distribution Analysis for All Variables",
            status = "primary", solidHeader = TRUE, width = 12,
            DT::dataTableOutput("all_variables_analysis")
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
            title = "Carbon Stock Measurements", 
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
            title = "ART TREES Monte Carlo Methodology with Distribution Analysis", 
            status = "info", solidHeader = TRUE, width = 12,
            h4("Implementation Based on ART TREES Standard V2.0 Section 8"),
            
            h5("Key Requirements:"),
            tags$ul(
              tags$li("Monte Carlo simulations: n=10,000 iterations, 90% CI"),
              tags$li("Error propagation between Activity Data and Emission Factors"),
              tags$li("Distribution analysis for input variables"),
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
            
            h5("Distribution Analysis Features:"),
            tags$ul(
              tags$li("Shapiro-Wilk normality tests for all carbon pools"),
              tags$li("Skewness and kurtosis assessment"),
              tags$li("CSV data upload capability"),
              tags$li("Automated hyperparameter recommendations based on distributions"),
              tags$li("Transformation suggestions for non-normal data")
            ),
            
            h5("Hyperparameter Optimization Features:"),
            tags$ul(
              tags$li("Model Selection: Linear, Random Forest, Gradient Boosting, SVM, Elastic Net"),
              tags$li("Cross-Validation: LGOCV Monte Carlo design with adaptive parameters"),
              tags$li("Data Augmentation: Synthetic samples for small datasets"),
              tags$li("Safe Parameter Grids: Constrained for dataset size"),
              tags$li("Bias Correction: Automated based on CV performance")
            )
          )
        )
      ),

      # Source Code Tab for Auditors
      tabItem(tabName = "sourcecode",
        fluidRow(
          box(
            title = "Monte Carlo Core Algorithm - ART TREES Equations", 
            status = "primary", solidHeader = TRUE, width = 12,
            h4("Uncertainty Propagation Algorithm"),
            verbatimTextOutput("monte_carlo_algorithm"),
            br(),
            h4("ART TREES Standard Calculations (Equations 10 & 11)"),
            verbatimTextOutput("art_trees_equations")
          )
        ),
        
        fluidRow(
          box(
            title = "Distribution Analysis Algorithm", 
            status = "info", solidHeader = TRUE, width = 6,
            h4("Shapiro-Wilk Normality Testing"),
            verbatimTextOutput("normality_test_code"),
            br(),
            h4("Skewness and Kurtosis Calculations"),
            verbatimTextOutput("distribution_metrics_code")
          ),
          box(
            title = "Hyperparameter Optimization Logic", 
            status = "warning", solidHeader = TRUE, width = 6,
            h4("Cross-Validation Implementation"),
            verbatimTextOutput("cv_implementation_code"),
            br(),
            h4("Bias Correction Algorithm"),
            verbatimTextOutput("bias_correction_code")
          )
        ),
        
        fluidRow(
          box(
            title = "Data Validation and Quality Checks", 
            status = "success", solidHeader = TRUE, width = 12,
            h4("Input Data Validation Rules"),
            verbatimTextOutput("validation_rules_code"),
            br(),
            h4("Computational Verification Steps"),
            verbatimTextOutput("verification_steps")
          )
        ),
        
        fluidRow(
          box(
            title = "Download Source Files", 
            status = "danger", solidHeader = TRUE, width = 12,
            p("Source code files are available for independent verification and audit purposes."),
            br(),
            downloadButton("download_source", "Download Complete Source Code", class = "btn-primary"),
            br(), br(),
            p("Repository Information:"),
            tags$ul(
              tags$li("Version: 2.0 - ART TREES Compliant"),
              tags$li("Last Updated: ", format(Sys.Date(), "%B %d, %Y")),
              tags$li("Compliance: ART TREES Standard V2.0 Section 8"),
              tags$li("Language: R version ", R.version.string)
            )
          )
        )
      )
    )
  )
)

# ============================================================================
# ENHANCED SERVER WITH DISTRIBUTION ANALYSIS
# ============================================================================

server <- function(input, output, session) {
  
  # Initialize data
  enhanced_data <- reactive({
    create_enhanced_guyana_data()
  })
  
  # Reactive values for distribution analysis
  current_analysis <- reactiveVal(NULL)
  all_analyses <- reactiveVal(NULL)  
  #all_analyses <- reactiveVal(list())
  
  # Auto-load distribution analysis on startup
	observe({
	  if(is.null(all_analyses())) {
	    if(file.exists("AllBiomassData_Prepared.csv")) {
	      biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
	      analysis_results <- data.frame()
	      
	      # Process all biomass columns
	      for(col_name in c("AGB_tC_ha", "BGB_tC_ha", "Saplings_tC_ha", 
	                        "Litter_tC_ha", "Standing_Dead_tC_ha", 
	                        "Lying_Dead_tC_ha", "Soil_tC_ha")) {
	        if(col_name %in% names(biomass_data)) {
	          data_vec <- as.numeric(biomass_data[, col_name])
	          data_clean <- na.omit(data_vec)
	          
	          if(length(data_clean) >= 3) {
	            sw_test <- shapiro.test(data_clean[1:min(5000, length(data_clean))])
	            
	            analysis_results <- rbind(analysis_results, data.frame(
	              Variable = col_name,
	              N = length(data_clean),
	              Mean = round(mean(data_clean), 2),
	              SD = round(sd(data_clean), 2),
	              Min = round(min(data_clean), 2),
	              Max = round(max(data_clean), 2),
	              W_statistic = round(sw_test$statistic, 4),
	              p_value = format(sw_test$p.value, digits = 4),
	              Normality = ifelse(sw_test$p.value > 0.05, "Normal", "Non-Normal"),
	              stringsAsFactors = FALSE
	            ))
	          }
	        }
	      }
	      
	      # Load LIF data
	      if(file.exists("LIF_Prepared.csv")) {
	        lif_data <- read.csv("LIF_Prepared.csv", stringsAsFactors = FALSE)
	        if("LIF_tC_km" %in% names(lif_data)) {
	          data_vec <- as.numeric(lif_data$LIF_tC_km)
	          data_clean <- na.omit(data_vec)
	          
	          if(length(data_clean) >= 3) {
	            sw_test <- shapiro.test(data_clean[1:min(5000, length(data_clean))])
	            
	            analysis_results <- rbind(analysis_results, data.frame(
	              Variable = "LIF_tC_km",
	              N = length(data_clean),
	              Mean = round(mean(data_clean), 2),
	              SD = round(sd(data_clean), 2),
	              Min = round(min(data_clean), 2),
	              Max = round(max(data_clean), 2),
	              W_statistic = round(sw_test$statistic, 4),
	              p_value = format(sw_test$p.value, digits = 4),
	              Normality = ifelse(sw_test$p.value > 0.05, "Normal", "Non-Normal"),
	              stringsAsFactors = FALSE
	            ))
	          }
	        }
	      }
	      
	      all_analyses(analysis_results)
	      if(nrow(analysis_results) > 0) {
	        current_analysis(analysis_results[1, , drop = FALSE])
	      }
	    }
	  }
	})

  # Variable mapping
  var_mapping <- list(
    "LIF" = list(name = "Logging Impact Factor", units = "tC/km"),
    "AGB" = list(name = "Above Ground Biomass", units = "tC/ha"),
    "BGB" = list(name = "Below Ground Biomass", units = "tC/ha"),
    "Saplings" = list(name = "Saplings Carbon", units = "tC/ha"),
    "Litter" = list(name = "Litter Carbon", units = "tC/ha"),
    "Standing_Dead" = list(name = "Standing Dead Wood", units = "tC/ha"),
    "Lying_Dead" = list(name = "Lying Dead Wood", units = "tC/ha"),
    "Soil" = list(name = "Soil Carbon", units = "tC/ha"),
    "LDF" = list(name = "Logging Degradation Factor", units = "tC/m³")
  )
  
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
  
  # Enhanced Monte Carlo results with hyperparameter optimization
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
  
  # ============================================================================
  # DISTRIBUTION ANALYSIS EVENT OBSERVERS 
  # ============================================================================

  # Distribution Analysis Event Observer
	observeEvent(input$analyze_distribution, {
	  withProgress(message = 'Analyzing distributions...', value = 0, {
	    
	    analysis_results <- data.frame()
	    
	    # Read AllBiomassData_Prepared.csv with CORRECT column names
	    if(file.exists("AllBiomassData_Prepared.csv")) {
	      biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
	      
	      # Debug: Print actual column names
	      cat("Actual column names in biomass data:", paste(names(biomass_data), collapse=", "), "\n")
	      
	      # Use EXACT column names from your CSV
	      biomass_columns <- c("AGB_tC_ha", "BGB_tC_ha", "Saplings_tC_ha", 
	                          "Litter_tC_ha", "Standing_Dead_tC_ha", 
	                          "Lying_Dead_tC_ha", "Soil_tC_ha")
	      
	      for(col_name in biomass_columns) {
	        if(col_name %in% names(biomass_data)) {
	          data_vec <- as.numeric(biomass_data[, col_name])
	          data_clean <- na.omit(data_vec)
	          
	          if(length(data_clean) >= 3) {
	            # Calculate statistics
	            sw_test <- shapiro.test(data_clean[1:min(5000, length(data_clean))])
	            skewness_val <- (sum((data_clean - mean(data_clean))^3) / length(data_clean)) / (sd(data_clean)^3)
	            
	            analysis_results <- rbind(analysis_results, data.frame(
	              Variable = col_name,
	              N = length(data_clean),
	              Mean = round(mean(data_clean), 2),
	              SD = round(sd(data_clean), 2),
	              Min = round(min(data_clean), 2),
	              Max = round(max(data_clean), 2),
	              Skewness = round(skewness_val, 3),
	              W_statistic = round(sw_test$statistic, 4),
	              p_value = format(sw_test$p.value, digits = 4),
	              Normality = ifelse(sw_test$p.value > 0.05, "Normal", "Non-Normal"),
	              stringsAsFactors = FALSE
	            ))
	          }
	        } else {
	          cat("Warning: Column", col_name, "not found in biomass data\n")
	        }
	      }
	    }
	    
	    # Read LIF_Prepared.csv with CORRECT column name
	    if(file.exists("LIF_Prepared.csv")) {
	      lif_data <- read.csv("LIF_Prepared.csv", stringsAsFactors = FALSE)
	      
	      # Debug: Print actual column names
	      cat("Actual column names in LIF data:", paste(names(lif_data), collapse=", "), "\n")
	      
	      if("LIF_tC_km" %in% names(lif_data)) {
	        data_vec <- as.numeric(lif_data$LIF_tC_km)
	        data_clean <- na.omit(data_vec)
	        
	        if(length(data_clean) >= 3) {
	          sw_test <- shapiro.test(data_clean[1:min(5000, length(data_clean))])
	          skewness_val <- (sum((data_clean - mean(data_clean))^3) / length(data_clean)) / (sd(data_clean)^3)
	          
	          analysis_results <- rbind(analysis_results, data.frame(
	            Variable = "LIF_tC_km",
	            N = length(data_clean),
	            Mean = round(mean(data_clean), 2),
	            SD = round(sd(data_clean), 2),
	            Min = round(min(data_clean), 2),
	            Max = round(max(data_clean), 2),
	            Skewness = round(skewness_val, 3),
	            W_statistic = round(sw_test$statistic, 4),
	            p_value = format(sw_test$p.value, digits = 4),
	            Normality = ifelse(sw_test$p.value > 0.05, "Normal", "Non-Normal"),
	            stringsAsFactors = FALSE
	          ))
	        }
	      }
	    }
	    
	    # Read LDF_Prepared.csv
	    if(file.exists("LDF_Prepared.csv")) {
	      ldf_data <- read.csv("LDF_Prepared.csv", stringsAsFactors = FALSE)
	      
	      # Debug: Print actual column names
	      cat("Actual column names in LDF data:", paste(names(ldf_data), collapse=", "), "\n")
	      
	      if("LDF_tC_m3" %in% names(ldf_data)) {
	        # Convert to numeric, handling any non-numeric values
	        data_vec <- suppressWarnings(as.numeric(ldf_data$LDF_tC_m3))
	        data_clean <- na.omit(data_vec)
	        
	        if(length(data_clean) >= 3) {
	          sw_test <- shapiro.test(data_clean[1:min(5000, length(data_clean))])
	          skewness_val <- (sum((data_clean - mean(data_clean))^3) / length(data_clean)) / (sd(data_clean)^3)
	          
	          analysis_results <- rbind(analysis_results, data.frame(
	            Variable = "LDF_tC_m3",
	            N = length(data_clean),
	            Mean = round(mean(data_clean), 2),
	            SD = round(sd(data_clean), 2),
	            Min = round(min(data_clean), 2),
	            Max = round(max(data_clean), 2),
	            Skewness = round(skewness_val, 3),
	            W_statistic = round(sw_test$statistic, 4),
	            p_value = format(sw_test$p.value, digits = 4),
	            Normality = ifelse(sw_test$p.value > 0.05, "Normal", "Non-Normal"),
	            stringsAsFactors = FALSE
	          ))
	        }
	      }
	    }
	    
	    # Store results
	    all_analyses(analysis_results)
	    if(nrow(analysis_results) > 0) {
	      current_analysis(analysis_results[1, , drop = FALSE])
	      showNotification(paste("Analyzed", nrow(analysis_results), "variables successfully"), type = "success")
	    } else {
	      showNotification("No data found to analyze. Check CSV files and column names.", type = "error")
	    }
	  })
	})
	
	# Auto-load on startup
	observe({
	  if(is.null(all_analyses())) {
	    # Trigger the same analysis as the button click
	    isolate({
	      analysis_results <- data.frame()
	      
	      if(file.exists("AllBiomassData_Prepared.csv")) {
	        biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
	        
	        biomass_columns <- c("AGB_tC_ha", "BGB_tC_ha", "Saplings_tC_ha", 
	                            "Litter_tC_ha", "Standing_Dead_tC_ha", 
	                            "Lying_Dead_tC_ha", "Soil_tC_ha")
	        
	        for(col_name in biomass_columns) {
	          if(col_name %in% names(biomass_data)) {
	            data_vec <- as.numeric(biomass_data[, col_name])
	            data_clean <- na.omit(data_vec)
	            
	            if(length(data_clean) >= 3) {
	              sw_test <- shapiro.test(data_clean[1:min(5000, length(data_clean))])
	              
	              analysis_results <- rbind(analysis_results, data.frame(
	                Variable = col_name,
	                N = length(data_clean),
	                Mean = round(mean(data_clean), 2),
	                SD = round(sd(data_clean), 2),
	                Min = round(min(data_clean), 2),
	                Max = round(max(data_clean), 2),
	                W_statistic = round(sw_test$statistic, 4),
	                p_value = format(sw_test$p.value, digits = 4),
	                Normality = ifelse(sw_test$p.value > 0.05, "Normal", "Non-Normal"),
	                stringsAsFactors = FALSE
	              ))
	            }
	          }
	        }
	      }
	      
	      if(file.exists("LIF_Prepared.csv")) {
	        lif_data <- read.csv("LIF_Prepared.csv", stringsAsFactors = FALSE)
	        if("LIF_tC_km" %in% names(lif_data)) {
	          data_vec <- as.numeric(lif_data$LIF_tC_km)
	          data_clean <- na.omit(data_vec)
	          
	          if(length(data_clean) >= 3) {
	            sw_test <- shapiro.test(data_clean[1:min(5000, length(data_clean))])
	            
	            analysis_results <- rbind(analysis_results, data.frame(
	              Variable = "LIF_tC_km",
	              N = length(data_clean),
	              Mean = round(mean(data_clean), 2),
	              SD = round(sd(data_clean), 2),
	              Min = round(min(data_clean), 2),
	              Max = round(max(data_clean), 2),
	              W_statistic = round(sw_test$statistic, 4),
	              p_value = format(sw_test$p.value, digits = 4),
	              Normality = ifelse(sw_test$p.value > 0.05, "Normal", "Non-Normal"),
	              stringsAsFactors = FALSE
	            ))
	          }
	        }
	      }
	      
	      all_analyses(analysis_results)
	    })
	  }
	})
  

  # Reset hyperparameters to stable defaults
  observeEvent(input$reset_hyperparameters, {
    updateSelectInput(session, "model_method", selected = "lm")
    updateSelectInput(session, "performance_metric", selected = "RMSE")
    updateNumericInput(session, "cv_folds", value = 5)
    updateNumericInput(session, "tune_length", value = 2)
    updateCheckboxInput(session, "enable_preprocessing", value = TRUE)
    updateCheckboxInput(session, "enable_feature_selection", value = FALSE)
    
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
  
  # ============================================================================
  # DISTRIBUTION ANALYSIS VALUE BOXES
  # ============================================================================
  
	output$normality_status_box <- renderValueBox({
    analysis <- current_analysis()
    if(is.null(analysis) || nrow(analysis) == 0) {
      valueBox("No Analysis", "Run analysis first", icon = icon("chart-bar"), color = "purple")
    } else {
      status <- analysis$Normality[1]
      color <- ifelse(status == "Normal", "green", "red")
      valueBox(status, "Distribution Status", icon = icon("bell-curve"), color = color)
    }
  })

  output$shapiro_p_value_box <- renderValueBox({
    analysis <- current_analysis()
    if(is.null(analysis) || nrow(analysis) == 0) {
      valueBox("--", "Shapiro-Wilk p-value", icon = icon("calculator"), color = "orange")
    } else {
      valueBox(
        analysis$p_value[1],
        "Shapiro-Wilk p-value",
        icon = icon("calculator"),
        color = "blue"
      )
    }
  })

  output$sample_size_box <- renderValueBox({
    analysis <- current_analysis()
    if(is.null(analysis) || nrow(analysis) == 0) {
      valueBox("--", "Sample Size", icon = icon("database"), color = "teal")
    } else {
      valueBox(
        analysis$N[1],
        "Sample Size",
        icon = icon("database"),
        color = "purple"
      )
    }
  })
  
  output$transformation_rec_box <- renderValueBox({
    analyses <- all_analyses()
    if(is.null(analyses) || nrow(analyses) == 0) {
      valueBox("--", "Transformation", icon = icon("exchange-alt"), color = "green")
    } else {
      # Get the first variable's results
      first_var <- analyses[1, ]
      
      # Determine recommendation based on skewness
      if("Skewness" %in% names(first_var) && !is.na(first_var$Skewness)) {
        skew_val <- as.numeric(first_var$Skewness)
        if(abs(skew_val) < 0.5) {
          rec <- "None Needed"
          color <- "green"
        } else if(skew_val > 0.5) {
          rec <- "Log/Sqrt"
          color <- "orange"
        } else {
          rec <- "Square/Cubic"
          color <- "orange"
        }
      } else {
        rec <- "None Needed"
        color <- "green"
      }
      
      valueBox(rec, "Transformation", icon = icon("exchange-alt"), color = color)
    }
  })
  
  # ============================================================================
  # DISTRIBUTION ANALYSIS OUTPUTS
  # ============================================================================
  
  output$distribution_plot <- renderPlot({
	  analyses <- all_analyses()
	  
	  if(!is.null(analyses) && nrow(analyses) > 0) {
	    # Create multi-panel plot for all variables
	    n_vars <- nrow(analyses)
	    n_cols <- ceiling(sqrt(n_vars))
	    n_rows <- ceiling(n_vars / n_cols)
	    
	    par(mfrow = c(n_rows, n_cols))
	    par(mar = c(3, 3, 3, 1))
	    
	    # Read data files once
	    biomass_data <- NULL
	    lif_data <- NULL
	    ldf_data <- NULL
	    
	    if(file.exists("AllBiomassData_Prepared.csv")) {
	      biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
	    }
	    if(file.exists("LIF_Prepared.csv")) {
	      lif_data <- read.csv("LIF_Prepared.csv", stringsAsFactors = FALSE)
	    }
	    if(file.exists("LDF_Prepared.csv")) {
	      ldf_data <- read.csv("LDF_Prepared.csv", stringsAsFactors = FALSE)
	    }
	    
	    for(i in 1:n_vars) {
	      var_name <- analyses$Variable[i]
	      
	      # Get the actual data for plotting
	      data_to_plot <- NULL
	      
	      # Check which file contains this variable
	      if(!is.null(biomass_data) && var_name %in% names(biomass_data)) {
	        data_to_plot <- as.numeric(biomass_data[, var_name])
	      } else if(!is.null(lif_data) && var_name == "LIF_tC_km") {
	        data_to_plot <- as.numeric(lif_data$LIF_tC_km)
	      } else if(!is.null(ldf_data) && var_name == "LDF_tC_m3") {
	        data_to_plot <- suppressWarnings(as.numeric(ldf_data$LDF_tC_m3))
	      }
	      
	      data_to_plot <- na.omit(data_to_plot)
	      
	      if(length(data_to_plot) > 0) {
	        # Create histogram with density overlay
	        hist(data_to_plot, 
	             main = paste(var_name, "\np =", analyses$p_value[i]),
	             xlab = "",
	             ylab = "Frequency",
	             col = ifelse(analyses$Normality[i] == "Normal", 
	                         rgb(0.2, 0.8, 0.2, 0.5), 
	                         rgb(0.8, 0.2, 0.2, 0.5)),
	             border = "darkgray",
	             breaks = 30)
	        
	        # Add normal curve overlay
	        if(length(data_to_plot) > 1) {
	          x_seq <- seq(min(data_to_plot), max(data_to_plot), length.out = 100)
	          y_norm <- dnorm(x_seq, mean = mean(data_to_plot), sd = sd(data_to_plot))
	          par(new = TRUE)
	          plot(x_seq, y_norm * length(data_to_plot) * diff(range(data_to_plot))/30, 
	               type = "l", col = "blue", lwd = 2, 
	               axes = FALSE, xlab = "", ylab = "")
	        }
	        
	        # Add mean line
	        abline(v = analyses$Mean[i], col = "red", lwd = 2, lty = 2)
	      } else {
	        # Empty plot if no data
	        plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
	        text(1, 1, paste(var_name, "\nNo data"), cex = 1)
	      }
	    }
	    
	    # Reset plotting parameters
	    par(mfrow = c(1, 1))
	    
	  } else {
	    # No analyses available
	    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
	    text(1, 1, "Click 'Analyze Distribution' to generate plots", cex = 1.5)
	  }
	})

  output$normality_summary <- renderText({
    analyses <- all_analyses()
    
    if(is.null(analyses) || nrow(analyses) == 0) {
      "No analysis performed yet. Click 'Analyze Distribution' to run tests."
    } else {
      # Get the first variable's results for display
      first_var <- analyses[1, ]
      
      paste0(
        "SHAPIRO-WILK NORMALITY TEST SUMMARY\n",
        "====================================\n\n",
        "Variable: ", first_var$Variable, "\n",
        "Sample Size: ", first_var$N, "\n",
        "Mean: ", first_var$Mean, "\n",
        "Std Dev: ", first_var$SD, "\n\n",
        "Test Statistic W: ", first_var$W_statistic, "\n",
        "p-value: ", first_var$p_value, "\n\n",
        "Interpretation:\n",
        if(as.numeric(first_var$p_value) > 0.05) {
          "Data appears normally distributed (p > 0.05)"
        } else if(as.numeric(first_var$p_value) > 0.01) {
          "Mild deviation from normality (0.01 < p < 0.05)"
        } else {
          "Significant deviation from normality (p < 0.01)"
        }, "\n\n",
        "Skewness: ", ifelse("Skewness" %in% names(first_var), first_var$Skewness, "N/A"), "\n\n",
        "Recommendation:\n",
        if(as.numeric(first_var$p_value) > 0.05) {
          "No transformation needed - data is approximately normal"
        } else if("Skewness" %in% names(first_var) && !is.na(first_var$Skewness)) {
          if(first_var$Skewness > 0.5) {
            "Consider log or square root transformation (right-skewed)"
          } else if(first_var$Skewness < -0.5) {
            "Consider square or cubic transformation (left-skewed)"
          } else {
            "Mild skewness - transformation may not be necessary"
          }
        } else {
          "Consider data transformation to achieve normality"
        }
      )
    }
  })
  
  output$descriptive_stats_table <- DT::renderDataTable({
    analyses <- all_analyses()
    
    if(!is.null(analyses) && nrow(analyses) > 0) {
      # Get the first variable for detailed stats
      first_var <- analyses[1, ]
      
      stats_df <- data.frame(
        Statistic = c("Sample Size", "Mean", "Std Dev", "Min", "Max", 
                     "Skewness", "W Statistic", "P-value", "Normality"),
        Value = c(
          as.character(first_var$N),
          as.character(first_var$Mean),
          as.character(first_var$SD),
          as.character(first_var$Min),
          as.character(first_var$Max),
          ifelse("Skewness" %in% names(first_var), as.character(first_var$Skewness), "N/A"),
          as.character(first_var$W_statistic),
          as.character(first_var$p_value),
          first_var$Normality
        ),
        stringsAsFactors = FALSE
      )
      
      DT::datatable(stats_df, 
                    options = list(dom = 't', pageLength = 15), 
                    rownames = FALSE)
    } else {
      DT::datatable(
        data.frame(Message = "No analysis available"),
        options = list(dom = 't'),
        rownames = FALSE
      )
    }
  })

  output$hyperparameter_recommendations <- renderText({
    analyses <- all_analyses()
    
    if(is.null(analyses) || nrow(analyses) == 0) {
      "Perform distribution analysis to receive hyperparameter recommendations."
    } else {
      # Count normal vs non-normal variables
      n_normal <- sum(analyses$Normality == "Normal")
      n_nonnormal <- sum(analyses$Normality == "Non-Normal")
      
      recommendations <- c(
        "HYPERPARAMETER RECOMMENDATIONS",
        "==============================\n",
        paste("Analysis Summary: ", n_normal, "normal,", n_nonnormal, "non-normal distributions\n")
      )
      
      # Model recommendations based on distribution mix
      if(n_nonnormal > n_normal) {
        recommendations <- c(recommendations,
          "• Non-linear models recommended due to non-normal distributions",
          "• Random Forest (rf) can handle skewed distributions well",
          "• Gradient Boosting (gbm) effective for complex patterns")
      } else {
        recommendations <- c(recommendations,
          "• Linear models (lm) appropriate for mostly normal data",
          "• Consider elastic net (glmnet) for regularization")
      }
      
      # CV strategy recommendations
      recommendations <- c(recommendations,
        "\n• Cross-validation strategy:",
        "  - Use 5-8 CV folds for this dataset size",
        "  - Enable preprocessing to normalize scales")
      
      # Specific variable recommendations
      if(any(analyses$Normality == "Non-Normal")) {
        non_normal_vars <- analyses$Variable[analyses$Normality == "Non-Normal"]
        recommendations <- c(recommendations,
          "\n• Variables requiring attention:",
          paste("  -", head(non_normal_vars, 3), "is non-normal"))
      }
      
      paste(recommendations, collapse = "\n")
    }
  })
  
  output$all_variables_analysis <- DT::renderDataTable({
    analyses <- all_analyses()
    
    if(is.null(analyses) || nrow(analyses) == 0) {
      DT::datatable(
        data.frame(Message = "Click 'Analyze Distribution' to run Shapiro-Wilk tests"),
        options = list(dom = 't'),
        rownames = FALSE
      )
    } else {
      DT::datatable(
        analyses,
        options = list(pageLength = 15, dom = 'frtip', scrollX = TRUE),
        rownames = FALSE
      ) %>%
      DT::formatStyle(
        "Normality",
        backgroundColor = DT::styleEqual(
          c("Normal", "Non-Normal"),
          c("#d4f1d4", "#ffd4d4")
        )
      )
    }
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
  
  # Data tables 
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
  
  # ============================================================================
  # SOURCE CODE TAB OUTPUTS FOR AUDITORS
  # ============================================================================
  
  output$monte_carlo_algorithm <- renderText({
    "# Monte Carlo Uncertainty Propagation (10,000 iterations)
# ART TREES Standard Section 8 Implementation

for(i in 1:n_strata) {
  # Extract uncertainty parameters
  t_90 <- 1.645006  # 90% confidence t-value
  
  # Calculate standard errors from 90% CI
  ad_se <- (AD_Mean[i] * AD_CI_90_percent[i] / 100) / t_90
  ef_se <- (EF_Mean[i] * EF_CI_90_percent[i] / 100) / t_90
  
  # Generate random samples (normal distribution)
  ad_samples <- rnorm(n = 10000, 
                      mean = AD_Mean[i], 
                      sd = ad_se)
  ef_samples <- rnorm(n = 10000, 
                      mean = EF_Mean[i], 
                      sd = ef_se)
  
  # Ensure non-negative values
  ad_samples <- pmax(ad_samples, 0)
  ef_samples <- pmax(ef_samples, 0)
  
  # Calculate emissions for each iteration
  emissions_matrix[, i] <- ad_samples * ef_samples
}

# Sum across all strata
total_emissions <- rowSums(emissions_matrix)

# Calculate statistics
simulated_mean <- mean(total_emissions)
ci_90_lower <- quantile(total_emissions, 0.05)
ci_90_upper <- quantile(total_emissions, 0.95)"
  })
  
  output$art_trees_equations <- renderText({
    "# ART TREES Standard Equations Implementation

# Equation 11: Uncertainty Adjustment Factor
UA_t = 0.524417 × (CI_90_percent / 100) / 1.645006

# Where:
# - 0.524417 = scaling constant from ART TREES
# - CI_90_percent = 90% confidence interval as percentage of mean
# - 1.645006 = t-value for 90% confidence

# Equation 10: Uncertainty Deduction
UNC_t = (GHGER_t + GHGREMV_t) × UA_t

# Where:
# - GHGER_t = Gross emission reductions
# - GHGREMV_t = Gross removals (if applicable)
# - UA_t = Uncertainty adjustment factor

# Net ERRs Calculation
Gross_ERRs = HFLD_Crediting_Level - Simulated_Mean_Emissions
Buffer_Deduction = Gross_ERRs × 0.05  # 5% buffer
Uncertainty_Deduction = Gross_ERRs × UA_t
Leakage_Deduction = 0  # Project-specific

Net_ERRs = Gross_ERRs - Buffer_Deduction - Uncertainty_Deduction - Leakage_Deduction"
  })
  
  output$normality_test_code <- renderText({
    "# Shapiro-Wilk Normality Test Implementation

shapiro_wilk_test <- function(data) {
  # Remove NA values
  data_clean <- na.omit(data)
  n <- length(data_clean)
  
  # Shapiro-Wilk requires 3 ≤ n ≤ 5000
  if(n < 3) {
    return(list(statistic = NA, p.value = NA))
  }
  
  # Sample if n > 5000 for computational efficiency
  if(n > 5000) {
    test_sample <- sample(data_clean, 5000)
  } else {
    test_sample <- data_clean
  }
  
  # Apply Shapiro-Wilk test
  result <- shapiro.test(test_sample)
  
  # Interpretation
  if(result$p.value > 0.05) {
    interpretation <- 'Data appears normally distributed'
  } else if(result$p.value > 0.01) {
    interpretation <- 'Mild deviation from normality'
  } else {
    interpretation <- 'Significant deviation from normality'
  }
  
  return(list(
    W = result$statistic,
    p_value = result$p.value,
    interpretation = interpretation
  ))
}"
  })
  
  output$distribution_metrics_code <- renderText({
    "# Statistical Moments Calculation

# Skewness (3rd moment)
calculate_skewness <- function(x) {
  n <- length(x)
  if(n < 3) return(NA)
  
  mean_x <- mean(x)
  sd_x <- sd(x)
  
  if(sd_x == 0) return(0)
  
  # Fisher-Pearson standardized moment coefficient
  skewness <- (sum((x - mean_x)^3) / n) / (sd_x^3)
  return(skewness)
}

# Kurtosis (4th moment)
calculate_kurtosis <- function(x) {
  n <- length(x)
  if(n < 4) return(NA)
  
  mean_x <- mean(x)
  sd_x <- sd(x)
  
  if(sd_x == 0) return(0)
  
  # Excess kurtosis (normal distribution = 0)
  kurtosis <- (sum((x - mean_x)^4) / n) / (sd_x^4) - 3
  return(kurtosis)
}

# Transformation recommendations
if(abs(skewness) < 0.5) {
  recommendation <- 'No transformation needed'
} else if(skewness > 0.5) {
  recommendation <- 'Log or square root transformation (right-skewed)'
} else {
  recommendation <- 'Square or cubic transformation (left-skewed)'
}"
  })
  
  output$cv_implementation_code <- renderText({
    "# Leave-Group-Out Cross-Validation (LGOCV) for Monte Carlo

# Generate synthetic samples for small dataset augmentation
samples_per_stratum <- 100

expanded_samples <- lapply(1:n_strata, function(i) {
  # Calculate standard errors
  ad_se <- (AD_Mean[i] * AD_CI_90[i] / 100) / 1.645006
  ef_se <- (EF_Mean[i] * EF_CI_90[i] / 100) / 1.645006
  
  # Generate Monte Carlo samples
  ad_samples <- rnorm(samples_per_stratum, AD_Mean[i], ad_se)
  ef_samples <- rnorm(samples_per_stratum, EF_Mean[i], ef_se)
  
  # Create feature matrix
  data.frame(
    emissions = ad_samples * ef_samples,
    ad_mean = ad_samples,
    ef_mean = ef_samples,
    ad_ci = AD_CI_90[i],
    ef_ci = EF_CI_90[i],
    stratum_id = factor(i)
  )
})

# LGOCV Configuration (ART TREES compliant)
cv_control <- trainControl(
  method = 'LGOCV',        # Monte Carlo CV
  number = 5,              # Number of resampling iterations
  p = 0.75,               # Training set percentage
  savePredictions = 'final',
  summaryFunction = defaultSummary
)"
  })
  
  output$bias_correction_code <- renderText({
    "# Bias Correction Factor Calculation

# Train model and evaluate performance
model <- train(
  emissions ~ ad_mean + ef_mean + ad_ci + ef_ci,
  data = expanded_samples,
  method = model_method,  # lm, rf, gbm, etc.
  trControl = cv_control,
  metric = 'RMSE'
)

# Extract best performance
best_rmse <- min(model$results$RMSE, na.rm = TRUE)
mean_emissions <- mean(expanded_samples$emissions)

# Calculate bias correction factor
if(mean_emissions > 0 && is.finite(best_rmse)) {
  bias_correction_factor <- 1 - (best_rmse / mean_emissions)
  
  # Conservative bounds [0.98, 1.02]
  bias_correction_factor <- max(0.98, min(1.02, bias_correction_factor))
} else {
  bias_correction_factor <- 1.0  # No correction
}

# Apply to Monte Carlo simulation
AD_Mean_corrected <- AD_Mean * bias_correction_factor
EF_Mean_corrected <- EF_Mean * bias_correction_factor"
  })
  
  output$validation_rules_code <- renderText({
    "# Input Data Validation Rules

validate_input_data <- function(data) {
  errors <- character()
  warnings <- character()
  
  # Rule 1: Check for required columns
  required_cols <- c('AD_Mean', 'AD_CI_90_percent', 
                    'EF_Mean', 'EF_CI_90_percent')
  missing_cols <- setdiff(required_cols, names(data))
  if(length(missing_cols) > 0) {
    errors <- c(errors, paste('Missing columns:', 
                             paste(missing_cols, collapse=', ')))
  }
  
  # Rule 2: Check for non-negative values
  if(any(data$AD_Mean < 0, na.rm = TRUE)) {
    errors <- c(errors, 'Activity Data cannot be negative')
  }
  if(any(data$EF_Mean < 0, na.rm = TRUE)) {
    errors <- c(errors, 'Emission Factors cannot be negative')
  }
  
  # Rule 3: Check CI ranges (0.1% to 200%)
  if(any(data$AD_CI_90_percent < 0.1 | 
         data$AD_CI_90_percent > 200, na.rm = TRUE)) {
    warnings <- c(warnings, 'AD CI outside typical range [0.1%, 200%]')
  }
  if(any(data$EF_CI_90_percent < 0.1 | 
         data$EF_CI_90_percent > 200, na.rm = TRUE)) {
    warnings <- c(warnings, 'EF CI outside typical range [0.1%, 200%]')
  }
  
  # Rule 4: Check for NA values
  na_count <- sum(is.na(data[, required_cols]))
  if(na_count > 0) {
    errors <- c(errors, paste(na_count, 'NA values found'))
  }
  
  return(list(errors = errors, warnings = warnings))
}"
  })
  
  output$verification_steps <- renderText({
    "# Computational Verification Steps

1. Monte Carlo Convergence Check:
   - Run with n = 1,000, 5,000, 10,000, 20,000 iterations
   - Verify mean stabilizes (< 1% change)
   - Verify CI converges (< 2% change)

2. Distribution Verification:
   - Shapiro-Wilk test on final emissions
   - Q-Q plot for normality assessment
   - Outlier detection (> 3 standard deviations)

3. Mathematical Verification:
   - Sum of strata emissions = Total emissions
   - UA factor ∈ [0, 1]
   - Net ERRs ≤ Gross ERRs
   - All deductions ≥ 0

4. ART TREES Compliance Check:
   - 90% confidence interval used
   - 10,000 Monte Carlo iterations
   - Equation 10 & 11 correctly applied
   - 5% buffer deduction included

5. Cross-Validation Integrity:
   - Training set ≠ Test set (no data leakage)
   - Performance metrics calculated on held-out data
   - Bias correction within conservative bounds

6. Reproducibility Check:
   - Set seed for random number generation
   - Document R version and package versions
   - Save intermediate results for audit trail"
  })
  
  # Download handler for source code
  output$download_source <- downloadHandler(
    filename = function() {
      paste0("ART_TREES_Source_Code_", Sys.Date(), ".R")
    },
    content = function(file) {
      # Read the current app.R file
      source_code <- readLines("app.R")
      writeLines(source_code, file)
    }
  )
}

# Run the enhanced application
shinyApp(ui = ui, server = server)