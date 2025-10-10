# ------------------------------------------------------------------------ #
# ART TREES Monte Carlo Simulation Tool
# Client: Guyana Forestry Commission
# Compliance: ART TREES Standard (V2.0) Section 8 
# Author: Murphy, S 
# Date: 2025-OCT-05
# ------------------------------------------------------------------------ #

pacman::p_load(
	BiocManager, bslib,
	caret,
	dplyr, DT,
	ellmer,
	ggplot2, gbm,
	MASS, magrittr,
	plotly,
	shiny, shinyWidgets, shinyjs, sessioninfo,
	thematic,
	devtools # Added as requested for logging functionality
) 
options(scipen = 999) 
options(digits = 6) 

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
			(42 * 1.645 / 322) * 100, 	# Forestry infrastructure: 21.5%
			(516 * 1.645 / 893) * 100, 	# Agriculture: 95.0%
			(1689 * 1.645 / 8448) * 100, # Mining: 32.8%
			(108 * 1.645 / 822) * 100, 	# Infrastructure: 21.6%
			(565 * 1.645 / 798) * 100, 	# Settlements: 116.5%
			(325 * 1.645 / 1491) * 100, 	# Fire-Biomass burning: 35.8%
			(410 * 1.645 / 1082) * 100, 	# Shifting Cultivation: 62.3%
			(45837 * 1.645 / 458366) * 100, # Logging volume: 16.4%
			(173 * 1.645 / 1733) * 100, 	# Skid trails: 16.4%
			(3054 * 1.645 / 30539) * 100 	# Mining buffer area: 16.4%
		),
		stringsAsFactors = FALSE
	)
	
	# Emission factors
	emission_factors <- data.frame(
		Activity_Type = activity_data$Activity_Type,
		Driver = activity_data$Driver,
		Units = activity_data$Units,
		EF_tCO2_unit = c(
			338528/322, 	# 1. Forestry infrastructure: 1,051.2 tCO2/ha
			937068/893, 	# 2. Agriculture: 1,049.5 tCO2/ha 	
			8881624/8448, 	# 3. Mining: 1,051.3 tCO2/ha
			864192/822, 	# 4. Infrastructure: 1,051.3 tCO2/ha
			839381/798, 	# 5. Settlements: 1,051.6 tCO2/ha
			1570410/1491, 	# 6. Fire-Biomass burning: 1,053.3 tCO2/ha
			999999/1082, 	# 7. Shifting cultivation: 924.4 tCO2/ha
			3.85, 			# 8. Logging - volume harvested (tCO2/m3)
			171.84, 		# 9. Logging - skid trail length (tCO2/km)
			247366/30539 	# 10. Mining buffer: 8.1 tCO2/ha
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
								   cv_folds = 5, 	# Reduced for small dataset
								   performance_metric = "RMSE",
								   tune_length = 2, 	# Reduced for small dataset
								   enable_preprocessing = TRUE,
								   enable_feature_selection = FALSE) { 	# Disabled by default for small dataset
	
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
					ad_ci = pmax(ad_ci_samples, 1.0), 	# Minimum 1% CI
					ef_ci = pmax(ef_ci_samples, 1.0),
					log_ad_mean = log(pmax(ad_samples, 1)),
					log_ef_mean = log(pmax(ef_samples, 1)),
					combined_uncertainty = sqrt(pmax(ad_ci_samples, 1)^2 + pmax(ef_ci_samples, 1)^2),
					stratum_id = factor(i)
				)
			}))
			
			cat("Generated", nrow(expanded_samples), "samples for LGOCV\n")
			
			# Conservative LGOCV setup for ART-TREES compliance
			cv_folds_safe <- min(cv_folds, 8) 	# Maximum 8 folds for stability
			train_percentage_safe <- 0.75 	# 75% found most reliable
			
			mc_control <- trainControl(
				method = "LGOCV", 	# MUST maintain Monte Carlo design
				number = cv_folds_safe,
								p = train_percentage_safe,
				savePredictions = "final",
				summaryFunction = defaultSummary,
				verboseIter = TRUE, 	# Enable to see what's happening
				allowParallel = FALSE,
				returnData = FALSE 	# Reduce memory usage
			)
			
			# ULTRA-CONSERVATIVE: Model-specific safe parameter grids
			if(model_method == "rf") {
				max_features <- ncol(expanded_samples) - 2 	# Exclude emissions and stratum_id
				safe_mtry <- c(2, min(4, max_features))
				tune_grid <- expand.grid(mtry = safe_mtry)
				cat("RF mtry range:", safe_mtry, "\n")
			} else if(model_method == "gbm") {
				# Ultra-simple GBM
				tune_grid <- expand.grid(
					n.trees = c(50, 100),
					interaction.depth = 1, 	# Keep it simple
					shrinkage = 0.1,
					n.minobsinnode = 20 	# Larger minimum
				)
			} else if(model_method == "svmRadial") {
				# Very simple SVM grid
				tune_grid <- expand.grid(sigma = 0.1, C = 1)
			} else if(model_method == "glmnet") {
				# Simple elastic net
				tune_grid <- expand.grid(alpha = c(0, 1), lambda = c(0.01, 0.1))
			} else {
				tune_grid <- NULL 	# Linear model - no tuning needed
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
					bias_correction_factor <- max(0.98, min(1.02, bias_correction_factor)) 	# Very conservative range
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
# NEW DISTRIBUTION ANALYSIS FUNCTIONS FROM BETA.R
# ============================================================================

# Load pre-existing CSV data files
load_preloaded_data <- function() {
	preloaded_data <- list()
	
	tryCatch({
		if(file.exists("LIF_Prepared.csv")) {
			lif_data <- read.csv("LIF_Prepared.csv", stringsAsFactors = FALSE)
			names(lif_data) <- gsub("\\s+", "_", names(lif_data))
			names(lif_data) <- gsub("[()]", "", names(lif_data))
			preloaded_data$LIF <- lif_data
		}
		
		if(file.exists("AllBiomassData_Prepared.csv")) {
			biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
			names(biomass_data) <- gsub("\\s+", "_", names(biomass_data))
			names(biomass_data) <- gsub("[()]", "", names(biomass_data))
			preloaded_data$Biomass <- biomass_data
		}
		
		if(file.exists("LDF_Prepared.csv")) {
			ldf_data <- read.csv("LDF_Prepared.csv", stringsAsFactors = FALSE)
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
	data_clean <- na.omit(data_vector)
	n <- length(data_clean)
	
	if(n < 3) {
		return(list(error = "Insufficient data for analysis (n < 3)"))
	}
	
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
	
	# Simple skewness calculation
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
	
	test_sample <- if(n > 5000) sample(data_clean, 5000) else data_clean
	shapiro_test <- if(n >= 3 && n <= 5000) {
		shapiro.test(test_sample)
	} else {
		list(statistic = NA, p.value = NA)
	}
	
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

# Generate synthetic variable data for demonstration
generate_synthetic_variable_data <- function(variable_name, units = "tC/ha") {
	set.seed(123)
	n <- 500
	
	data <- switch(variable_name,
		"LIF" = rlnorm(n, meanlog = 3, sdlog = 0.5),
		"AGB" = rgamma(n, shape = 20, rate = 0.1),
		"BGB" = rgamma(n, shape = 10, rate = 0.2),
		"Saplings" = rgamma(n, shape = 2, rate = 0.5),
		"Litter" = rnorm(n, mean = 3.3, sd = 0.8),
		"Standing_Dead" = c(rep(0, n*0.3), rgamma(n*0.7, shape = 1, rate = 0.3)),
		"Lying_Dead" = c(rep(0, n*0.2), rgamma(n*0.8, shape = 2, rate = 0.2)),
		"Soil" = rnorm(n, mean = 60, sd = 20) + abs(rnorm(n, 0, 10)),
		"LDF" = rnorm(n, mean = 8, sd = 2),
		rnorm(n, mean = 100, sd = 30)
	)
	
	return(abs(data))
}

# Extract variable from CSV based on column mapping
extract_variable_from_csv <- function(csv_data, variable_name, file_path) {
	possible_names <- c(
		variable_name,
		tolower(variable_name),
		toupper(variable_name),
		gsub("_", " ", variable_name),
		gsub("_", "", variable_name)
	)
	
	col_idx <- which(tolower(names(csv_data)) %in% tolower(possible_names))
	
	if(length(col_idx) > 0) {
		return(as.numeric(csv_data[, col_idx[1]]))
	} else {
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
		csv_data <- read.csv(file_path, stringsAsFactors = FALSE)
		var_data <- extract_variable_from_csv(csv_data, variable_name, file_path)
		
		if(!is.null(var_data)) {
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
	
	if(abs(skewness) > 1) {
		recommendations <- c(recommendations,
			paste("• High skewness (", round(skewness, 2), ") - consider data transformation"))
	}
	
	if(cv > 50) {
		recommendations <- c(recommendations,
			"• High variability - use more CV folds (8-10) for robust estimates",
			"• Enable preprocessing to normalize scales")
	} else {
		recommendations <- c(recommendations,
			"• Moderate variability - standard CV folds (5) should suffice")
	}
	
	if(kurtosis > 3 || kurtosis < -3) {
		recommendations <- c(recommendations,
			"• Extreme kurtosis detected - consider feature engineering",
			"• Add polynomial features or interaction terms")
	}
	
	return(paste(recommendations, collapse = "\n"))
}

# ============================================================================
# ENHANCED UI WITH DISTRIBUTION ANALYSIS TAB - BSLIB DESIGN
# ============================================================================

ui <- page_fillable(
	title = "ART TREES Monte Carlo Tool: Normality + Hyperparameter Tuning",
	fillable_mobile = FALSE,
	
	theme = bs_theme(
		version = 5,
		bootswatch = "darkly",
		base_font = font_google("Inter"),
		heading_font = font_google("Poppins"),
		code_font = font_google("Fira Code")
	),
	
	tags$head(
		tags$style(HTML("
			.card-header { padding: 0.5rem !important; font-size: 0.9rem !important; }
			.bslib-value-box .value-box-title { font-size: 0.8rem !important; }
			.bslib-value-box .value-box-value { font-size: 1.2rem !important; }
			.dataTables_wrapper { font-size: 0.85rem !important; }
			.nav-pills .nav-link { padding: 0.3rem 0.8rem !important; font-size: 0.85rem !important; }
			pre { white-space: pre-wrap; word-wrap: break-word; }
			#code_tools_dropdown {
				position: absolute;
				top: 10px;
				right: 10px;
				z-index: 1000; /* FIXED: Corrected CSS property separation from = to : */
			}
		"))
	),

	# CODE TOOLS WIDGET
	tags$div(id = "code_tools_dropdown",
		shinyWidgets::dropdown(
			tagList(
				tags$h5("Code Tools & Auditing", style="color:#007bff; margin-top: 5px;"),
				hr(style="margin: 5px 0;"),
				actionButton("link_to_source", "Expose/View Source Code", class = "btn-secondary btn-sm", width = "100%", 
							 onclick = "Shiny.setInputValue('main_tabs', 'Source Code', {priority: 'event'});"),
				downloadButton("download_full_source", "Download App.R File", class = "btn-info btn-sm", width = "100%"),
				hr(style="margin: 5px 0;"),
				actionButton("link_to_github", "Edit on GitHub", class = "btn-secondary btn-sm", width = "100%", 
							 onclick = "window.open('https://github.com/seamusrobertmurphy/TREES-monte-carlo-app.git', '_blank')"),
				actionButton("link_to_issues", "Report an Issue", class = "btn-danger btn-sm", width = "100%", 
							 onclick = "window.open('https://github.com/seamusrobertmurphy/TREES-monte-carlo-app/issues', '_blank')"),
				hr(style="margin: 5px 0;"),
				fileInput("upload_data_inline", "Upload New Data:", 
						  accept = c("text/csv", ".csv"), 
						  buttonLabel = "Upload CSV", 
						  placeholder = "Choose file",
						  width = "100%")
			),
			icon = icon("tools"),
			size = "xs",
			up = TRUE,
			status = "primary",
			width = "200px" 
			# FIXED: Removed 'style = "padding: 5px;"' which caused the match.arg error
		)
	),

	layout_sidebar(
		fillable = TRUE,
		
		sidebar = sidebar(
			width = 300,
			bg = "steelblue",
			tags$li(class = "dropdown",
							tags$a(
								img(src = "https://winrock.org/wp-content/uploads/2021/12/Winrock-logo-R.png", height = "40px", style = "padding: 5px;"),
								href = "https://winrock.org/program-areas/environment-and-climate/ecosystem-services/",
								target = "_blank"
								)
							),
			h6("ART TREES Parameters"),
			tags$small(
				tags$div(style = "line-height: 1.2;",
					"• MC Simulations: 10,000", br(),
					"• Confidence Interval: 90%", br(),
					"• t-value: 1.645006"
				)
			),
			actionButton("run_simulation", "Run Monte Carlo", 
						 class = "btn-success btn-sm", style = "width: 100%;"),
			actionButton("reset_hyperparameters", "Reset Defaults", 
						 class = "btn-secondary btn-sm", style = "width: 100%"),
			
			tags$hr(style="margin: 0.5rem 0;"),
			h6("Simulation Regime"),
			checkboxInput("use_bootstrap", "Use Bootstrap", value = FALSE),
			
			tags$hr(style="margin: 0.5rem 0;"),
			h6("Hyperparameter Tuning"),
			selectInput("model_method", "Model:",
						choices = list("Linear" = "lm", "RF" = "rf", "GBM" = "gbm",
									   "SVM" = "svmRadial", "Elastic Net" = "glmnet"),
						selected = "lm"),
			
			selectInput("performance_metric", "Metric:",
						choices = list("RMSE" = "RMSE", "MAE" = "MAE", "R-squared" = "Rsquared"),
						selected = "RMSE"),
			
			numericInput("cv_folds", "CV Folds:", value = 5, min = 3, max = 10),
			numericInput("tune_length", "Grid Size:", value = 2, min = 1, max = 5),
			
			checkboxInput("enable_preprocessing", "Preprocessing", value = TRUE),
			checkboxInput("enable_feature_selection", "Feature Selection", value = FALSE),
			
			tags$hr(style="margin: 0.5rem 0;"),
			h6("Uncertainty Reduction"),
			selectInput("selected_stratum", "Driver:", choices = NULL),
			numericInput("adjust_ad_ci", "AD CI (%):", value = 95.1, min = 0.1, max = 200),
			numericInput("adjust_ef_ci", "EF CI (%):", value = 6.9, min = 0.1, max = 50)
		),
		
		navset_card_pill(
			id = "main_tabs",
			selected = "2. Distribution",
			
			# TAB 1: Input Data Overview
			nav_panel(
				title = "1. Input Data",
				icon = icon("table"),
				layout_columns(
					col_widths = c(6, 6),
					card(
						card_header("Carbon Stock Measurements"), 
						card_body(
							tags$p("Based on field measurements from 118 plots with 472 subplots (4 per plot)", style = "font-size: 0.7rem;"),
							DT::dataTableOutput("carbon_stocks_table")
						)
					),
					card(
						card_header("Activity Data (Deforestation & Degradation 2024)"), 
						card_body(
							tags$p("Annual deforestation and forest degradation activities with logging operations", style = "font-size: 0.7rem;"),
							DT::dataTableOutput("activity_data_table")
						)
					)
				),
				card(
					card_header("Monte Carlo Input Data Structure"), 
					card_body(
						tags$p("Complete dataset structure for Monte Carlo uncertainty simulation with hyperparameter optimization", style = "font-size: 0.7rem;"),
						DT::dataTableOutput("monte_carlo_input_table")
					)
				)
			),
			
# TAB 2: Distribution
nav_panel(
  title = "2. Distribution",
  icon = icon("chart-bar"),
  
  # Value boxes row
  layout_columns(
    col_widths = c(3, 3, 3, 3),
    fill = FALSE,
    uiOutput("normality_box_output"),
    value_box(
      title = "Shapiro p-value", 
      value = textOutput("shapiro_p_value_text"),
      showcase = icon("calculator"), 
      theme = "info"
    ),
    value_box(
      title = "Sample Size", 
      value = textOutput("sample_size_text"),
      showcase = icon("database"), 
      theme = "secondary"
    ),
    value_box(
      title = "Transformation", 
      value = textOutput("transformation_rec_text"),
      showcase = icon("arrows-rotate"), 
      theme = "success"
    )
  ),

  # Controls + Distribution Plot + Normality Summary
  layout_columns(
    col_widths = c(3, 5, 4),  # ADJUSTED: 3 for controls, 5 for plot, 4 for summary
    fill = FALSE,
    
    card(
      card_header("Controls"),
      card_body(
        fillable = FALSE,
        selectInput("analysis_variable", "Variable:",
          choices = list(
            "LIF tC/km" = "LIF",
            "AGB tC/ha" = "AGB",
            "BGB tC/ha" = "BGB", 
            "Saplings tC/ha" = "Saplings",
            "Litter tC/ha" = "Litter",
            "Standing Dead tC/ha" = "Standing_Dead",
            "Lying Dead tC/ha" = "Lying_Dead",
            "Soil tC/ha" = "Soil",
            "LDF tC/m³" = "LDF"
          ),
          selected = "Litter"
        ),
        checkboxInput("use_preloaded_data", "Use Pre-loaded", value = TRUE)
      )
    ),
    
    card(
      full_screen = TRUE,
      card_header("Distribution Plot"),
      card_body(
        fillable = FALSE,
        plotOutput("distribution_plot", height = "400px")  # INCREASED HEIGHT
      )
    ),
    
    card(
      card_header("Normality Summary"),
      card_body(
        fillable = FALSE,
        div(style = "height: 400px; overflow-y: auto; font-size: 0.75rem;",
          verbatimTextOutput("normality_summary")
        )
      )
    )
  ),
  
  # NEW: Descriptive Statistics + All Variables Analysis SIDE-BY-SIDE
  layout_columns(
    col_widths = c(6, 6),  # 50-50 split
    fill = FALSE,
    card(
      card_header("Descriptive Statistics"),
      card_body(
        fillable = FALSE,
        div(style = "height: 450px; overflow-y: auto;",
          DT::dataTableOutput("descriptive_stats_table")
        )
      )
    ),
    card(
      card_header("All Variables Analysis"),
      card_body(
        fillable = FALSE,
        div(style = "height: 450px; overflow-y: auto;",
          DT::dataTableOutput("all_variables_analysis")
        )
      )
    )
  )
),

# TAB 3: Hyperparameter Tuning
nav_panel(
  title = "3. Hyperparameters",
  icon = icon("cogs"),
  
  # Value boxes (existing)
  layout_columns(
    col_widths = c(3, 3, 3, 3),
    fill = FALSE,
    value_box(title = "Optimization Model", value = textOutput("tuning_method_text"),
              showcase = icon("cogs"), theme = "primary"),
    value_box(title = "Best RMSE", value = textOutput("best_performance_text"),
              showcase = icon("bullseye"), theme = "warning"),
    value_box(title = "Bias Correction", value = textOutput("bias_correction_text"),
              showcase = icon("balance-scale"), theme = "secondary"),
    value_box(title = "Tuning", value = textOutput("optimization_status_text"),
              showcase = icon("check-circle"), theme = "info")
  ),
  
  # Current run results (existing)
  layout_columns(
    col_widths = c(8, 4),
    fill = FALSE,
    card(
      card_header("Hyperparameter Optimization Results"),
      card_body(
        fillable = FALSE,
        DT::dataTableOutput("hyperparameter_results_table")
      )
    ),
    card(
      card_header("Model Performance Comparison"),
      card_body(
        fillable = FALSE,
        plotOutput("performance_comparison_plot", height = "300px")
      )
    )
  ),
  
  # CV Summary + Variable Importance (existing)
  layout_columns(
    col_widths = c(6, 6),
    fill = FALSE,
    card(
      card_header("Cross-Validation Summary"),
      card_body(
        fillable = FALSE,
        div(style = "height: 300px; overflow-y: auto; font-size: 0.85rem;",
          verbatimTextOutput("cv_summary")
        )
      )
    ),
    card(
      card_header("Variable Importance"),
      card_body(
        fillable = FALSE,
        plotOutput("variable_importance_plot", height = "300px")
      )
    )
  ),
  
  # NEW: Model History Comparison Section
  card(
    card_header("Multi-Model Comparison History"),
    card_body(
      layout_columns(
        col_widths = c(4, 4, 4),
        fill = FALSE,
        actionButton("clear_history", "Clear History", 
                     class = "btn-warning btn-sm", style = "width: 100%;"),
        downloadButton("download_model_history", "Export History", 
                      class = "btn-info btn-sm", style = "width: 100%;"),
        tags$div(
          style = "padding-top: 8px; text-align: center;",
          textOutput("history_count")
        )
      )
    )
  ),
  
  # Model comparison visualizations
  layout_columns(
    col_widths = c(12),
    fill = FALSE,
    card(
      card_header("Model History Table"),
      card_body(
        fillable = FALSE,
        DT::dataTableOutput("model_history_table")
      )
    )
  ),
  
  layout_columns(
    col_widths = c(6, 6),
    fill = FALSE,
    card(
      card_header("Performance Comparison Across Runs"),
      card_body(
        fillable = FALSE,
        plotOutput("model_comparison_plot", height = "400px")
      )
    ),
    card(
      card_header("Net ERRs Comparison Across Runs"),
      card_body(
        fillable = FALSE,
        plotOutput("net_errs_comparison_plot", height = "400px")
      )
    )
  )
),
			
			# TAB 4: Monte Carlo
			nav_panel(
				title = "4. Monte Carlo",
				icon = icon("chart-line"),
				
				layout_columns(
					col_widths = c(3, 3, 3, 3),
					fill = FALSE,
					value_box(title = "Mean Emissions", value = textOutput("total_emissions_text"),
							  showcase = icon("fire"), theme = "danger"),
					value_box(title = "Uncertainty", value = textOutput("uncertainty_deduction_text"),
							  showcase = icon("minus-circle"), theme = "warning"),
					value_box(title = "Net ERRs", value = textOutput("net_err_text"),
							  showcase = icon("leaf"), theme = "success"),
					value_box(title = "90% CI", value = textOutput("ci_percent_text"),
							  showcase = icon("chart-bar"), theme = "info")
				),
				
				layout_columns(
					col_widths = c(8, 4),
					fill = FALSE,
					card(
						card_header("Monte Carlo Distribution"),
						card_body(
							fillable = FALSE, 
							plotOutput("enhanced_distribution_plot", 
													 height = "450px"))
					),
					card(
						card_header("Uncertainty Summary"),
						card_body(
							fillable = FALSE, 
							DT::dataTableOutput("uncertainty_summary")
							)
						)
					),
				layout_columns(
					col_widths = c(6, 6),
					fill = FALSE,
					card(
						card_header("Stratum-Level Results"),
						card_body(fillable = FALSE, DT::dataTableOutput("stratum_results"))
					),
					card(
						card_header("Sensitivity Analysis"),
						card_body(fillable = FALSE, verbatimTextOutput("enhanced_sensitivity_analysis"))
					)
				),
				
				card(
						card_header("ART TREES Calculation Verification with Hyperparameter Optimization"),
						card_body(
								div(style = "overflow-x: auto;",
										verbatimTextOutput("calculation_details")
								)
						)
				)
			),
			
			# TAB 5: Methodology Guide
			nav_panel(
				title = "Guide",
				icon = icon("book"),
				card(
					card_header("ART TREES Monte Carlo Methodology with Distribution Analysis"),
					card_body(
						tags$h4("Implementation Based on ART TREES Standard V2.0 Section 8"),
						
						tags$h5("Key Requirements:"),
						tags$ul(
							tags$li("Monte Carlo simulations: n=10,000 iterations, 90% CI"),
							tags$li("Error propagation between Activity Data and Emission Factors"),
							tags$li("Distribution analysis for input variables"),
							tags$li("Hyperparameter optimization for uncertainty reduction"),
							tags$li("Cross-validation for bias detection and model optimization")
						),
						
						tags$h5("Mathematical Framework:"),
						withMathJax(),
						tags$p("Equation 11 - Uncertainty Adjustment Factor:"),
						tags$p("$$UA_t = 0.524417 \\times \\frac{HW_{90\\%}}{1.645006}$$"),
						tags$p("Equation 10 - Uncertainty Deduction:"),
						tags$p("$$UNC_t = (GHGER_t + GHGREMV_t) \\times UA_t$$"),
						tags$p("Net ERRs Calculation:"),
						tags$p("$$Net\\ ERRs = Gross\\ ERRs - Buffer\\ (5\\%) - Leakage\\ (0\\%) - Uncertainty\\ Deduction$$"),
						
						tags$h5("Distribution Analysis Features:"),
						tags$ul(
							tags$li("Shapiro-Wilk normality tests for all carbon pools"),
							tags$li("Skewness and kurtosis assessment"),
							tags$li("Automated hyperparameter recommendations based on distributions"),
							tags$li("Transformation suggestions for non-normal data")
						),
						
						tags$h5("Hyperparameter Optimization Features:"),
						tags$ul(
							tags$li("Model Selection: Linear, Random Forest, Gradient Boosting, SVM, Elastic Net"),
							tags$li("Cross-Validation: LGOCV Monte Carlo design with adaptive parameters"),
							tags$li("Data Augmentation: Synthetic samples for small datasets"),
							tags$li("Bias Correction: Automated based on CV performance")
						)
					)
				)
			),
			
			# TAB 6: Source Code
			nav_panel(
				title = "Source Code",
				icon = icon("code"),
				
				layout_columns(
					col_widths=c(6,6),
					
					# Card 1 Runtime log
					card(
						card_header("Runtime Log"),
					  card_body(
						div(style = "height: 550px; overflow-y: auto; font-size: 0.5rem;",
						    verbatimTextOutput("session_info_log")
								)
						)),
					# Card 2 Render Source
					card(
						card_header("Runtime Source"),
						card_body(
							div(style = "height: 550px; overflow-y: scroll; white-space: pre-wrap; font-family: 'Fira Code', monospace; font-size: 0.6rem; background-color: #2b3238; padding: 1px; border-radius: 1px; color: #ced4da;",
								verbatimTextOutput("runtime_source_code")
							)
						)
					)
				),
				

				

				# Card 3: Code Tools and Download/Edit Links (Simplified)
				card(
					card_body(style = "padding: 0.75rem 0.5rem;",
						layout_columns(
							col_widths = c(4, 4, 4),

						# Simple Download button
						downloadButton("download_report", "Download Script & Log", class = "btn-primary btn-sm", style = "width: 100%; margin: 1px;"),
						
						# GitHub Edit
						actionButton("link_to_github_main", "Contributor Tools", class = "btn-secondary", style = "width: 100%; margin: 1px;",
									 onclick = "window.open('https://github.com/seamusrobertmurphy/TREES-monte-carlo-app', '_blank')"),
						# Issue Tracker
						actionButton("link_to_issues_main", "Report an Issue", class = "btn-danger", style = "width: 100%; margin: 1px;",
									 onclick = "window.open('https://github.com/seamusrobertmurphy/TREES-monte-carlo-app/issues', '_blank')")
						)
					)
				)
			)
		)
	)
)


# Enable thematic
thematic::thematic_shiny(font = "auto")
theme_set(theme_bw(base_size = 16))


# ============================================================================
# ENHANCED SERVER WITH DISTRIBUTION ANALYSIS
# ============================================================================

server <- function(input, output, session) {
	
	enhanced_data <- reactive({
		create_enhanced_guyana_data()
	})
	
	current_analysis <- reactiveVal(NULL)
	all_analyses <- reactiveVal(NULL)
	
  # NEW: Store history of model runs for comparison
  model_history <- reactiveVal(data.frame(
    run_id = integer(),
    timestamp = character(),
    model_method = character(),
    cv_folds = integer(),
    performance_metric = character(),
    best_performance = numeric(),
    bias_correction = numeric(),
    mean_emissions = numeric(),
    ci_percent = numeric(),
    net_errs = numeric(),
    use_bootstrap = logical(),
    stringsAsFactors = FALSE
  ))
  
	# Auto-load distribution analysis on startup
	observe({
		priority = 100

		if(is.null(all_analyses())) {
			isolate({
				if(file.exists("AllBiomassData_Prepared.csv")) {
					biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
					analysis_results <- data.frame()
					
					biomass_columns <- c("AGB_tC_ha", "BGB_tC_ha", "Saplings_tC_ha", 
										 "Litter_tC_ha", "Standing_Dead_tC_ha", 
										 "Lying_Dead_tC_ha", "Soil_tC_ha")
					
					for(col_name in biomass_columns) {
						if(col_name %in% names(biomass_data)) {
							data_vec <- as.numeric(biomass_data[, col_name])
							data_clean <- na.omit(data_vec)

							# FILTER ZEROS for variables that commonly have zeros
							if(col_name %in% c("Litter_tC_ha", "Standing_Dead_tC_ha", "Lying_Dead_tC_ha")) {
								data_clean_filtered <- data_clean[data_clean > 0]
								if(length(data_clean_filtered) >= 3) {
									data_clean <- data_clean_filtered
								} else {
									next
								}
							}
							
							if(length(data_clean) >= 3) {
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
					
					# Load LDF data
					if(file.exists("LDF_Prepared.csv")) {
						ldf_data <- read.csv("LDF_Prepared.csv", stringsAsFactors = FALSE)
						if("LDF_tC_m3" %in% names(ldf_data)) {
							data_vec <- suppressWarnings(as.numeric(ldf_data$LDF_tC_m3))
							data_clean <- na.omit(data_vec)

							# FILTER ZEROS for LDF
							data_clean_filtered <- data_clean[data_clean > 0]
							if(length(data_clean_filtered) >= 3) {
								data_clean <- data_clean_filtered
							} else {
								next
							}
							
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
					
					all_analyses(analysis_results)
					# Set current analysis to Litter by default
					litter_analysis <- analysis_results[analysis_results$Variable == "Litter_tC_ha", ]
					if(nrow(litter_analysis) > 0) {
						current_analysis(litter_analysis[1, , drop = FALSE])
					} else if(nrow(analysis_results) > 0) {
						current_analysis(analysis_results[1, , drop = FALSE])
					}
				}
			})
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
    
    # NEW: Save to model history
    current_history <- model_history()
    new_run <- data.frame(
      run_id = nrow(current_history) + 1,
      timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
      model_method = input$model_method,
      cv_folds = input$cv_folds,
      performance_metric = input$performance_metric,
      best_performance = if(result$tuning_applied && !is.null(result$tuning_results)) {
        min(result$tuning_results$cv_results[[input$performance_metric]], na.rm = TRUE)
      } else {
        result$simulated_sd  # Use SD as proxy for non-tuned models
      },
      bias_correction = result$bias_correction_factor,
      mean_emissions = result$simulated_mean,
      ci_percent = result$ci_90_percent_of_mean,
      net_errs = result$net_errs,
      use_bootstrap = input$use_bootstrap,
      stringsAsFactors = FALSE
    )
    
    model_history(rbind(current_history, new_run))
    
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
			
			# Read AllBiomassData_Prepared.csv
			if(file.exists("AllBiomassData_Prepared.csv")) {
				biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
				
				cat("Actual column names in biomass data:", paste(names(biomass_data), collapse=", "), "\n")
				
				biomass_columns <- c("AGB_tC_ha", "BGB_tC_ha", "Saplings_tC_ha", 
									 "Litter_tC_ha", "Standing_Dead_tC_ha", 
									 "Lying_Dead_tC_ha", "Soil_tC_ha")
				
				for(col_name in biomass_columns) {
					if(col_name %in% names(biomass_data)) {
						data_vec <- as.numeric(biomass_data[, col_name])
						data_clean <- na.omit(data_vec)

						# FILTER ZEROS
						if(col_name %in% c("Litter_tC_ha", "Standing_Dead_tC_ha", "Lying_Dead_tC_ha")) {
							data_clean_filtered <- data_clean[data_clean > 0]
							if(length(data_clean_filtered) >= 3) {
								data_clean <- data_clean_filtered
								cat("Filtered zeros for:", col_name, ". New N:", length(data_clean), "\n")
							} else {
								cat("Warning: Filtering zeros on", col_name, "resulted in < 3 samples.\n")
								next
							}
						}
						
						if(length(data_clean) >= 3) {
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
			
			# Read LIF_Prepared.csv
			if(file.exists("LIF_Prepared.csv")) {
				lif_data <- read.csv("LIF_Prepared.csv", stringsAsFactors = FALSE)
				
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
				
				cat("Actual column names in LDF data:", paste(names(ldf_data), collapse=", "), "\n")
				
				if("LDF_tC_m3" %in% names(ldf_data)) {
					data_vec <- suppressWarnings(as.numeric(ldf_data$LDF_tC_m3))
					data_clean <- na.omit(data_vec)

					# FILTER ZEROS for LDF
					data_clean_filtered <- data_clean[data_clean > 0]
					if(length(data_clean_filtered) >= 3) {
						data_clean <- data_clean_filtered
					} else {
						next
					}
					
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
			
			# Update current_analysis to the selected variable
			selected_var_name <- switch(input$analysis_variable,
										"LIF" = "LIF_tC_km",
										"AGB" = "AGB_tC_ha",
										"BGB" = "BGB_tC_ha",
										"Saplings" = "Saplings_tC_ha",
										"Litter" = "Litter_tC_ha",
										"Standing_Dead" = "Standing_Dead_tC_ha",
										"Lying_Dead" = "Lying_Dead_tC_ha",
										"Soil" = "Soil_tC_ha",
										"LDF" = "LDF_tC_m3",
										"Litter_tC_ha")
			
			current_row <- analysis_results[analysis_results$Variable == selected_var_name, ]
			if(nrow(current_row) > 0) {
				current_analysis(current_row[1, , drop = FALSE])
				showNotification(paste("Analyzed", nrow(analysis_results), "variables successfully"), type = "success")
			} else if(nrow(analysis_results) > 0) {
				current_analysis(analysis_results[1, , drop = FALSE])
				showNotification("Selected variable not found. Displaying first analyzed variable.", type = "warning")
			} else {
				showNotification("No data found to analyze. Check CSV files and column names.", type = "error")
			}
		})
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
	

# Clear model history
observeEvent(input$clear_history, {
  model_history(data.frame(
    run_id = integer(),
    timestamp = character(),
    model_method = character(),
    cv_folds = integer(),
    performance_metric = character(),
    best_performance = numeric(),
    bias_correction = numeric(),
    mean_emissions = numeric(),
    ci_percent = numeric(),
    net_errs = numeric(),
    use_bootstrap = logical(),
    stringsAsFactors = FALSE
  ))
  showNotification("Model history cleared", type = "message")
})

# Download model history
output$download_model_history <- downloadHandler(
  filename = function() {
    paste0("ART_TREES_Model_History_", Sys.Date(), ".csv")
  },
  content = function(file) {
    write.csv(model_history(), file, row.names = FALSE)
  }
)

# NEW: Model History Comparison Table
output$model_history_table <- DT::renderDataTable({
  history <- model_history()
  
  if(nrow(history) == 0) {
    return(DT::datatable(
      data.frame(Message = "Run simulations to build model comparison history"),
      options = list(dom = 't'),
      rownames = FALSE
    ))
  }
  
  # Format the display
  display_history <- history
  display_history$best_performance <- format(round(display_history$best_performance, 0), big.mark = ",")
  display_history$mean_emissions <- format(round(display_history$mean_emissions, 0), big.mark = ",")
  display_history$net_errs <- format(round(display_history$net_errs, 0), big.mark = ",")
  display_history$ci_percent <- paste0(round(display_history$ci_percent, 2), "%")
  display_history$bias_correction <- round(display_history$bias_correction, 4)
  
  DT::datatable(
    display_history,
    options = list(
      pageLength = 10,
      dom = 'frtip',
      scrollX = TRUE,
      order = list(list(0, 'desc'))  # Sort by run_id descending (newest first)
    ),
    rownames = FALSE
  ) %>%
    DT::formatStyle(
      'model_method',
      backgroundColor = DT::styleEqual(
        c('lm', 'rf', 'gbm', 'svmRadial', 'glmnet'),
        c('#e3f2fd', '#fff3e0', '#f3e5f5', '#e8f5e9', '#fce4ec')
      )
    )
})

# NEW: Model Comparison Plot
output$model_comparison_plot <- renderPlot({
  history <- model_history()
  
  if(nrow(history) == 0) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "Run multiple simulations\nto see model comparisons", cex = 1.0)
    return()
  }
  
  if(nrow(history) == 1) {
    par(mar = c(4, 4, 3, 1))
    barplot(history$best_performance[1],
            names.arg = history$model_method[1],
            col = rgb(0, 0.48, 1, 0.7),
            border = "darkblue",
            main = paste("Model Performance (", history$performance_metric[1], ")", sep = ""),
            ylab = history$performance_metric[1],
            las = 1)
    return()
  }
  
  par(mar = c(6, 4, 3, 1))
  
  # Create comparison plot
  x_labels <- paste0(history$model_method, "\n(Run ", history$run_id, ")")
  colors <- rainbow(nrow(history), alpha = 0.7)
  
  bp <- barplot(history$best_performance,
          names.arg = x_labels,
          col = colors,
          border = "darkblue",
          main = paste("Model Performance Comparison (", history$performance_metric[1], ")", sep = ""),
          ylab = history$performance_metric[1],
          las = 2,
          cex.names = 0.8,
          cex.main = 1.0)
  
  # Add values on bars
  text(bp, history$best_performance, 
       labels = format(round(history$best_performance, 0), big.mark = ","),
       pos = 3, cex = 0.7, srt = 90, adj = 0)
  
  # Highlight best model
  best_idx <- which.min(history$best_performance)
  points(bp[best_idx], history$best_performance[best_idx], 
         col = "red", pch = 8, cex = 3, lwd = 2)
  
  legend("topright",
         legend = c("Best Model"),
         pch = 8,
         col = "red",
         pt.cex = 2,
         bty = "n")
  
}, height = 400)

output$history_count <- renderText({
  paste("Models Run:", nrow(model_history()))
})

# NEW: Net ERRs Comparison Plot
output$net_errs_comparison_plot <- renderPlot({
  history <- model_history()
  
  if(nrow(history) == 0) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "Run multiple simulations\nto compare Net ERRs", cex = 1.0)
    return()
  }
  
  par(mar = c(6, 4, 3, 1))
  
  x_labels <- paste0(history$model_method, "\n(Run ", history$run_id, ")")
  colors <- rainbow(nrow(history), alpha = 0.7)
  
  bp <- barplot(history$net_errs,
          names.arg = x_labels,
          col = colors,
          border = "darkgreen",
          main = "Net ERRs Comparison Across Models",
          ylab = "Net ERRs (tCO2/yr)",
          las = 2,
          cex.names = 0.8,
          cex.main = 1.0,
          ylim = c(0, max(history$net_errs) * 1.1))
  
  # Add values on bars
  text(bp, history$net_errs, 
       labels = format(round(history$net_errs, 0), big.mark = ","),
       pos = 3, cex = 0.7, srt = 90, adj = 0)
  
  # Highlight best (maximum) Net ERRs
  best_idx <- which.max(history$net_errs)
  points(bp[best_idx], history$net_errs[best_idx], 
         col = "darkgreen", pch = 8, cex = 3, lwd = 2)
  
  legend("topright",
         legend = c("Highest Net ERRs"),
         pch = 8,
         col = "darkgreen",
         pt.cex = 2,
         bty = "n")
  
}, height = 400)


	
	# Code Tools outputs
	output$download_full_source <- downloadHandler(
		filename = function() { "app.R" },	
		content = function(file) {	
			# This is a safe placeholder, as the execution environment often prevents 
			# reading the app.R file itself for dynamic content delivery.
			writeLines(
				c("# ART-TREES Monte Carlo Tool Source Code",
				  "# Full code contents are not directly readable in this runtime environment.",
				  "# Please copy the code from the 'Source Code' tab for auditing."),
				file)
		}
	)

	observeEvent(input$upload_data_inline, {
		req(input$upload_data_inline)
		showNotification(paste("File uploaded:", input$upload_data_inline$name), type = "message")
	})

	# DYNAMIC THEME STATUS RENDERER
	output$normality_box_output <- renderUI({
		analysis <- current_analysis()
		
		if (is.null(analysis) || nrow(analysis) == 0) {
			theme_color <- "secondary"
			status_text <- "No Analysis"
		} else {
			status_text <- analysis$Normality[1]
			theme_color <- if (status_text == "Normal") {
				"success"
			} else if (status_text == "Non-Normal") {
				"danger"
			} else {
				"secondary"
			}
		}

		bslib::value_box(
			title = "Normality",
			value = status_text,
			showcase = icon("chart-line"),
			theme = theme_color
		)
	})

	# Monte Carlo Simulation value boxes
	output$total_emissions_text <- renderText({
		req(mc_results())
		format(round(mc_results()$simulated_mean, 0), big.mark = ",")
	})

	output$uncertainty_deduction_text <- renderText({
		req(mc_results())
		format(round(mc_results()$uncertainty_deduction, 0), big.mark = ",")
	})

	output$net_err_text <- renderText({
		req(mc_results())
		format(round(mc_results()$net_errs, 0), big.mark = ",")
	})

	output$ci_percent_text <- renderText({
		req(mc_results())
		paste0(round(mc_results()$ci_90_percent_of_mean, 1), "%")
	})


	
	# Distribution Analysis value boxes
	output$shapiro_p_value_text <- renderText({
		analysis <- current_analysis()
		if(is.null(analysis) || nrow(analysis) == 0) {
			"--"
		} else {
			format(as.numeric(analysis$p_value[1]), scientific = FALSE)
		}
	})

	output$sample_size_text <- renderText({
		analysis <- current_analysis()
		if(is.null(analysis) || nrow(analysis) == 0) {
			"--"
		} else {
			as.character(analysis$N[1])
		}
	})

	output$transformation_rec_text <- renderText({
		analysis <- current_analysis()
		if (is.null(analysis) || nrow(analysis) == 0) {
			return("--")
		}	
		
		first_var <- analysis[1, ]
		if("Skewness" %in% names(first_var) && !is.na(first_var$Skewness)) {
			skew_val <- as.numeric(first_var$Skewness)
			if(abs(skew_val) < 0.5) {
				"None Needed"
			} else if(skew_val > 0.5) {
				"Log/Sqrt"
			} else {
				"Square/Cubic"
			}
		} else {
			"None Needed"
		}
	})

	# Hyperparameter Tuning value boxes
	output$tuning_method_text <- renderText({
		req(mc_results())
		results <- mc_results()
		if(results$tuning_applied) results$hyperparameters_used$model_method else "None"
	})

	output$best_performance_text <- renderText({
		req(mc_results())
		results <- mc_results()
		if(results$tuning_applied && !is.null(results$tuning_results)) {
			metric_name <- results$hyperparameters_used$performance_metric
			best_val <- min(results$tuning_results$cv_results[[metric_name]], na.rm = TRUE)
			paste0(round(best_val, 3), " (", metric_name, ")")
		} else {
				"N/A"
		}
	})

	output$bias_correction_text <- renderText({
		req(mc_results())
		as.character(round(mc_results()$bias_correction_factor, 4))	
	})

	output$optimization_status_text <- renderText({
		req(mc_results())
		if(mc_results()$tuning_applied) "✓ Active" else "✗ Disabled"
	})
	
	# ============================================================================
	# DISTRIBUTION ANALYSIS OUTPUTS
	# ============================================================================
	
	output$distribution_plot <- renderPlot({
		analysis <- current_analysis()
		req(analysis)
		
		var_col_map <- list(
			"LIF" = "LIF_tC_km",
			"AGB" = "AGB_tC_ha",
			"BGB" = "BGB_tC_ha",
			"Saplings" = "Saplings_tC_ha",
			"Litter" = "Litter_tC_ha",
			"Standing_Dead" = "Standing_Dead_tC_ha",
			"Lying_Dead" = "Lying_Dead_tC_ha",
			"Soil" = "Soil_tC_ha",
			"LDF" = "LDF_tC_m3"
		)
		
		current_input_var <- input$analysis_variable
		if (is.null(current_input_var) || !(current_input_var %in% names(var_col_map))) {
				current_input_var <- "Litter"
		}
		
		selected_var_name <- var_col_map[[current_input_var]]
		
		data_to_plot <- NULL
		
		if(current_input_var %in% c("AGB", "BGB", "Saplings", "Litter", "Standing_Dead", "Lying_Dead", "Soil")) {
			if(file.exists("AllBiomassData_Prepared.csv")) {
				biomass_data <- read.csv("AllBiomassData_Prepared.csv", stringsAsFactors = FALSE)
				data_to_plot <- as.numeric(biomass_data[, selected_var_name])
			}
		} else if (current_input_var == "LIF") {
			if(file.exists("LIF_Prepared.csv")) {
				lif_data <- read.csv("LIF_Prepared.csv", stringsAsFactors = FALSE)
				data_to_plot <- as.numeric(lif_data$LIF_tC_km)
			}
		} else if (current_input_var == "LDF") {
				if(file.exists("LDF_Prepared.csv")) {
				ldf_data <- read.csv("LDF_Prepared.csv", stringsAsFactors = FALSE)
				data_to_plot <- suppressWarnings(as.numeric(ldf_data$LDF_tC_m3))
			}
		}

		data_to_plot <- na.omit(data_to_plot)

		# Apply zero filter
		if(current_input_var %in% c("Litter", "Standing_Dead", "Lying_Dead", "LDF")) {
				data_to_plot <- data_to_plot[data_to_plot > 0]
		}
		
		if(length(data_to_plot) == 0) {
			plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
			text(1, 1, paste("No non-zero data available for", current_input_var), cex = 0.8)
			return()
		}
		
		var_stats <- analysis
		units <- var_mapping[[current_input_var]]$units
		
		par(mar = c(4, 4, 2, 1))
		
		bins = 30
		h <- hist(data_to_plot, breaks = bins, plot = FALSE)

		p_val_display_plot <- format(as.numeric(var_stats$p_value), scientific = FALSE)
		
		hist(data_to_plot, 
			 main = paste(current_input_var, "Distribution\np =", p_val_display_plot),
			 xlab = paste(current_input_var, "(", units, ")"),
			 ylab = "Frequency",
			 col = ifelse(var_stats$Normality == "Normal", 
						  rgb(0.2, 0.8, 0.2, 0.5), 
						  rgb(0.8, 0.2, 0.2, 0.5)),
			 border = "darkgray",
			 breaks = h$breaks,
			 las = 1,
			 cex.main = 1,
			 cex.lab = 0.9)
		
		if(length(data_to_plot) > 1 && sd(data_to_plot) > 0) {
			x_seq <- seq(min(data_to_plot), max(data_to_plot), length.out = 100)
			y_norm <- dnorm(x_seq, mean = mean(data_to_plot), sd = sd(data_to_plot))
			scale_factor <- max(h$counts) / max(y_norm * diff(h$breaks)[1])	
			lines(x_seq, y_norm * diff(h$breaks)[1] * scale_factor, 
				  col = "blue", lwd = 2)
		}
		
		abline(v = mean(data_to_plot), col = "red", lwd = 2, lty = 2)
		
		legend("topright", 
			   legend = c("Data", "Normal Fit", "Mean"),
			   col = c(NA, "blue", "red"),
			   lty = c(NA, 1, 2),
			   lwd = c(NA, 2, 2),
			   fill = c(ifelse(var_stats$Normality == "Normal", 
							   rgb(0.2, 0.8, 0.2, 0.5), 
							   rgb(0.8, 0.2, 0.2, 0.5)), NA, NA),
			   border = c("darkgray", NA, NA),
			   bty = "n",
			   cex = 0.8)
	}, height = 400)

	output$normality_summary <- renderText({
		analysis <- current_analysis()
		
		if(is.null(analysis) || nrow(analysis) == 0) {
			"No analysis performed yet. Select a variable and click 'Analyze'."
		} else {
			first_var <- analysis[1, ]
			p_val <- as.numeric(first_var$p_value)
			p_val_display <- format(p_val, scientific = FALSE)

			interpretation <- if(p_val > 0.05) {
				"Data appears normally distributed (p > 0.05). **This supports using simpler statistical models (like Linear Models) or simple error propagation for uncertainty estimation.**"
			} else if(p_val > 0.01) {
				"Mild deviation from normality (0.01 < p < 0.05). **Consider using data transformation (e.g., log, square root) or robustness checks to ensure model validity.**"
			} else {
				"Significant deviation from normality (p < 0.01). **Highly recommended to use non-parametric or non-linear models (like Random Forest or GBM) for uncertainty assessment to avoid bias in standard Monte Carlo simulation results.**"
			}
			
			skewness_rec <- if("Skewness" %in% names(first_var) && !is.na(first_var$Skewness)) {
				skew_val <- as.numeric(first_var$Skewness)
				if(p_val > 0.05) {
					"No transformation needed based on Shapiro-Wilk test result."
				} else if(abs(skew_val) > 1.0) {
					if(skew_val > 0) {
						"High positive skew: **Log or Square Root transformation is strongly recommended.**"
					} else {
						"High negative skew: **Square or Cubic transformation is strongly recommended.**"
					}
				} else if(abs(skew_val) > 0.5) {
					"Moderate skew: **Transformation is advisable if non-normal models cannot be used.**"
				} else {
					"Mild skewness: Transformation likely not necessary if a robust model is used."
				}
			} else {
				"Unable to calculate skewness."
			}
					
			paste0(
				"SHAPIRO-WILK NORMALITY TEST SUMMARY\n",
				"====================================\n\n",
				"Variable: ", first_var$Variable, "\n",
				"Sample Size: ", first_var$N, "\n",
				"Mean: ", first_var$Mean, "\n",
				"Std Dev: ", first_var$SD, "\n\n",
				"Test Statistic W: ", first_var$W_statistic, "\n",
				"p-value: ", p_val_display, "\n\n",
				"Interpretation:\n",
				interpretation, "\n\n",
				"Skewness (Fisher): ", ifelse("Skewness" %in% names(first_var), first_var$Skewness, "N/A"), "\n\n",
				"**Monte Carlo Recommendation:**\n",
				skewness_rec
			)
		}
	})
	
	output$descriptive_stats_table <- DT::renderDataTable({
			analyses <- all_analyses()
			
			if(!is.null(analyses) && nrow(analyses) > 0) {
				current_var_name <- current_analysis()$Variable[1]
				first_var <- analyses[analyses$Variable == current_var_name, ]
				
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

	output$all_variables_analysis <- DT::renderDataTable({
		analyses <- all_analyses()
		
		if(is.null(analyses) || nrow(analyses) == 0) {
			DT::datatable(
				data.frame(Message = "Click 'Analyze' to run Shapiro-Wilk tests"),
				options = list(dom = 't'),
				rownames = FALSE
			)
		} else {
			DT::datatable(
				analyses,
				options = list(pageLength = 10, dom = 'frtip', scrollX = TRUE),
				rownames = FALSE
			) %>%
			DT::formatStyle(
				"Normality",
				backgroundColor = DT::styleEqual(
					c("Normal", "Non-Normal"),
					c("#28a745", "#dc3545"),
					default = "#6c757d"
				)
			)
		}
	})

	# ============================================================================
	# MONTE CARLO OUTPUTS
	# ============================================================================
	
output$enhanced_distribution_plot <- renderPlot({
  # Validate inputs
  validate(
    need(mc_results(), "Run Monte Carlo simulation to see results")
  )
  
  results <- mc_results()
  
  validate(
    need(!is.null(results$total_emissions), "No emissions data available"),
    need(length(results$total_emissions) > 0, "Emissions data is empty")
  )
  
  # Extract emissions data
  emissions <- results$total_emissions
  
  # Set up plot margins
  par(mar = c(4, 4, 3, 1))
  
  # Create histogram
  h <- hist(emissions,
       breaks = 50,
       col = rgb(0, 0.48, 1, 0.7),
       border = "darkblue",
       main = paste("Enhanced Monte Carlo Results\nModel:", 
                    if(results$tuning_applied) results$hyperparameters_used$model_method else "Standard"),
       xlab = "Total Emissions (tCO2/yr)",
       ylab = "Frequency",
       las = 1,
       cex.main = 1.0,
       cex.lab = 0.9)
  
  # Add vertical lines for mean and CI
  abline(v = results$simulated_mean, col = "#dc3545", lwd = 2, lty = 2)
  abline(v = results$ci_90_lower, col = "#ffc107", lwd = 2, lty = 3)
  abline(v = results$ci_90_upper, col = "#ffc107", lwd = 2, lty = 3)
  
  # Add legend
  legend("topright",
         legend = c(
           paste("Mean:", format(round(results$simulated_mean, 0), big.mark = ",")),
           paste("90% CI: [", format(round(results$ci_90_lower, 0), big.mark = ","), ",",
                 format(round(results$ci_90_upper, 0), big.mark = ","), "]")
         ),
         col = c("#dc3545", "#ffc107"),
         lty = c(2, 3),
         lwd = 2,
         bty = "n",
         cex = 0.75)
  
}, height = 400, res = 96)
	
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
			"• Equation 11 - UA Factor = 0.524417 \\times (", round(results$ci_90_percent_of_mean, 2), "% / 100) / 1.645006 = ", 
			round(results$ua_factor, 6), "\n",
			"• Equation 10 - Uncertainty Deduction = Gross ERRs \\times UA Factor = ", 
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
			"SENSITIVITY ANALYSIS: ", stratum_name, "\n",
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
	# DATA OVERVIEW TAB OUTPUTS
	# ============================================================================
	
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
					  options = list(scrollX = TRUE, pageLength = 10, dom = 'frtip'),
					  rownames = FALSE) %>%
			DT::formatStyle(
				"Activity Type",
				target = "row",
				backgroundColor = DT::styleEqual(c("Deforestation", "Degradation"), 
												c("#505962", "#454d55"), 
												default = "#222d32")
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
		
		dt <- DT::datatable(display_data, options = list(pageLength = 10, dom = 'tp', scrollX = TRUE), rownames = FALSE) %>%
			DT::formatStyle(
				"Component",
				target = "row",
				backgroundColor = DT::styleEqual(c("Sum C pools w/o litter", "Sum C pools w/ litter", "Sum ALL POOLS"), 
												c("#333b43", "#30373e", "#2b3238"))
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
		
		dt <- DT::datatable(display_data, options = list(pageLength = 10, dom = 'tp', scrollX = TRUE), rownames = FALSE) %>%
			DT::formatStyle(
				"Activity Type",
				target = "row",
				backgroundColor = DT::styleEqual(c("Deforestation", "Degradation"), 
												c("#505962", "#454d55"))
			)
		dt
	})
	
# ============================================================================
# HYPERPARAMETER TUNING OUTPUTS
# ============================================================================

# Hyperparameter results table
output$hyperparameter_results_table <- DT::renderDataTable({
  req(mc_results())
  results <- mc_results()
  
  if(results$tuning_applied && !is.null(results$tuning_results)) {
    cv_results <- results$tuning_results$cv_results
    DT::datatable(cv_results, 
                  options = list(scrollX = TRUE, pageLength = 10, dom = 'frtip'), 
                  rownames = FALSE)
  } else {
    DT::datatable(
      data.frame(Message = "No hyperparameter optimization results available (Model: Standard/Failed)"),
      options = list(dom = 't'),
      rownames = FALSE
    )
  }
})

# Performance comparison plot - FINAL VERSION
output$performance_comparison_plot <- renderPlot({
  
  if(is.null(mc_results())) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "Run Monte Carlo\nsimulation first", cex = 1.0)
    return()
  }
  
  results <- mc_results()
  
  if(!results$tuning_applied) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "No hyperparameter tuning applied.\nUsing standard Monte Carlo.", cex = 1.0)
    return()
  }
  
  if(is.null(results$tuning_results) || is.null(results$tuning_results$cv_results)) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "CV results not available", cex = 1.0)
    return()
  }
  
  cv_results <- results$tuning_results$cv_results
  metric_name <- results$hyperparameters_used$performance_metric
  
  # Linear models only have one configuration (intercept TRUE/FALSE)
  if(nrow(cv_results) <= 1) {
    par(mar = c(4, 4, 3, 1))
    
    # Create a simple bar plot showing the single performance value
    if(metric_name %in% names(cv_results)) {
      perf_value <- cv_results[[metric_name]][1]
      
      barplot(perf_value,
              names.arg = "Model\nPerformance",
              col = rgb(0, 0.48, 1, 0.7),
              border = "darkblue",
              main = paste("Cross-Validation", metric_name),
              ylab = metric_name,
              las = 1,
              cex.main = 1.0,
              ylim = c(0, perf_value * 1.1))
      
      # Add value on top of bar
      text(1, perf_value, 
           labels = format(round(perf_value, 0), big.mark = ","),
           pos = 3, cex = 0.8, col = "darkblue", font = 2)
      
    } else {
      plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
      text(1, 1, "Single configuration:\nNo comparison needed.", cex = 1.0)
    }
    return()
  }
  
  # Multiple configurations - create comparison plot
  if(!metric_name %in% names(cv_results)) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, paste("Metric", metric_name, "\nnot found"), cex = 1.0)
    return()
  }
  
  par(mar = c(4, 4, 3, 1))
  
  x_vals <- 1:nrow(cv_results)
  y_vals <- cv_results[[metric_name]]
  
  plot(x_vals, y_vals,
       type = "b",
       col = "blue",
       lwd = 2,
       pch = 19,
       cex = 1.5,
       main = paste("Hyperparameter Performance by", metric_name),
       xlab = "Configuration",
       ylab = metric_name,
       las = 1,
       cex.main = 1.0,
       cex.lab = 0.9)
  
  grid(nx = NULL, ny = NULL, col = "lightgray", lty = "dotted")
  
  best_idx <- if(metric_name %in% c("RMSE", "MAE")) which.min(y_vals) else which.max(y_vals)
  points(x_vals[best_idx], y_vals[best_idx], col = "red", pch = 19, cex = 2)
  
  legend("topright",
         legend = c("Performance", "Best"),
         col = c("blue", "red"),
         pch = 19,
         cex = 0.8,
         bty = "n")
  
}, height = 300)


# Variable importance plot - FINAL VERSION
output$variable_importance_plot <- renderPlot({
  
  if(is.null(mc_results())) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "Run Monte Carlo simulation first", cex = 1.0)
    return()
  }
  
  results <- mc_results()
  
  if(!results$tuning_applied) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "Variable importance\nnot available:\nNo hyperparameter tuning applied.", cex = 1.0)
    return()
  }
  
  # Linear models don't produce variable importance
  model_method <- results$hyperparameters_used$model_method
  
  if(model_method %in% c("lm", "glm")) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, paste0("Variable importance not available\nfor linear models (", model_method, ").\n\nTry Random Forest or GBM\nfor variable importance."), cex = 0.9)
    return()
  }
  
  # For other models, try to get variable importance
  if(is.null(results$tuning_results$variable_importance)) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "Variable importance\nnot calculated for this model.", cex = 1.0)
    return()
  }
  
  var_imp <- results$tuning_results$variable_importance
  
  tryCatch({
    if(inherits(var_imp, "varImp.train") || "importance" %in% names(var_imp)) {
      imp_df <- if("importance" %in% names(var_imp)) {
        as.data.frame(var_imp$importance)
      } else {
        as.data.frame(var_imp)
      }
      
      if(nrow(imp_df) > 0 && ncol(imp_df) > 0) {
        par(mar = c(4, 8, 3, 1))
        
        imp_vals <- imp_df[, 1]
        var_names <- rownames(imp_df)
        
        ord <- order(imp_vals)
        imp_vals <- imp_vals[ord]
        var_names <- var_names[ord]
        
        barplot(imp_vals,
                names.arg = var_names,
                horiz = TRUE,
                las = 1,
                col = rgb(0, 0.48, 1, 0.7),
                border = "darkblue",
                main = "Variable Importance",
                xlab = "Importance",
                cex.names = 0.8,
                cex.main = 1.0)
      }
    } else {
      plot(var_imp, main = "Variable Importance")
    }
  }, error = function(e) {
    plot(1, 1, type = "n", axes = FALSE, xlab = "", ylab = "")
    text(1, 1, "Error plotting\nvariable importance", cex = 1.0)
  })
  
}, height = 300)

# CV Summary - Already correct, no changes needed
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
      "Training Percentage: 0.75\n",
      "Performance Metric: ", results$hyperparameters_used$performance_metric, "\n",
      "Preprocessing: ", if(results$hyperparameters_used$enable_preprocessing) "Enabled (Center/Scale)" else "Disabled", "\n",
      "Feature Selection: ", if(results$hyperparameters_used$enable_feature_selection) "Enhanced" else "Basic", "\n",
      "Bias Correction Applied: ", round(results$bias_correction_factor, 4), "\n\n",
      "**Impact:** This optimization helps train an uncertainty-reduction model\n",
      "based on synthetic Monte Carlo samples, aiming to reduce the\n",
      "final uncertainty deduction (UA) as allowed by ART TREES V2.0."
    )
  } else {
    paste0(
      "Cross-Validation Summary\n",
      "========================\n\n",
      "Model Method: Standard Monte Carlo\n",
      "Status: Optimization not applied (e.g., small dataset, LM model)\n",
      "Bias Correction: 1.0 (No correction)\n\n",
      "The simulation runs in standard mode, calculating uncertainty\n",
      "directly from input CIs using error propagation rules."
    )
  }
})
	
	# ============================================================================
	# SOURCE CODE TAB OUTPUTS FOR AUDITORS
	# ============================================================================
	

	# NEW: Session Info Log
	output$session_info_log <- renderPrint({
  # Use base R sessionInfo() instead of sessioninfo package
  sessionInfo()
	})


	# NEW: Render the Source Code Content
	output$runtime_source_code <- renderText({
		# Reads and renders the content of the running app.R file.
		tryCatch({
			lines <- readLines("app.R")
			return(paste(lines, collapse = "\n"))
		}, error = function(e) {
			return(paste("Error reading source file:", e$message, 
						 "\n(This usually occurs when the Shiny environment cannot access 'app.R' directly, or the file is not in the working directory.)"))
		})
	})
	
# Report Download Handler (HTML/Formatted Text)
	output$download_report <- downloadHandler(
		filename = function() {
			# Date-stamped HTML filename
			paste0("ART_TREES_Audit_Report_", Sys.Date(), ".html")
		},
		content = function(file) {
			
			# --- 1. Capture Runtime Log ---
			runtime_log <- tryCatch({
				paste(capture.output(sessionInfo()), collapse = "\n")
			}, error = function(e) {
				"ERROR: Failed to capture sessionInfo()."
			})
			
			# --- 2. Capture Source Code ---
			source_code <- tryCatch({
				paste(readLines("app.R"), collapse = "\n")
			}, error = function(e) {
				"ERROR: Could not read 'app.R'. Check file permissions/path."
			})
			
			# --- 3. Compile Full HTML Report ---
			report_html <- paste0(
				"<!DOCTYPE html>",
				"<html lang='en'>",
				"<head>",
				"<meta charset='UTF-8'>",
				"<title>ART TREES Monte Carlo Audit Report</title>",
				"<style>",
				"body { font-family: sans-serif; margin: 20px; background-color: #f4f4f9; color: #333; }",
				"h1, h2 { color: #007bff; border-bottom: 2px solid #007bff; padding-bottom: 5px; }",
				"pre { background-color: #ffffff; color: #000000; border: 1px solid #ccc; padding: 15px; border-radius: 5px; overflow-x: auto; font-family: 'Fira Code', monospace; line-height: 1.2; font-size: 0.9em; }",
				"strong { color: #dc3545; }",
				".metadata { margin-bottom: 20px; border: 1px solid #ddd; padding: 10px; border-radius: 5px; background-color: #fff; }",
				"</style>",
				"</head>",
				"<body>",
				
				"<h1>ART TREES Monte Carlo Audit Report</h1>",
				"<div class='metadata'>",
				"<p><strong>Generated Date:</strong> ", Sys.time(), "</p>",
				"<p><strong>Report Type:</strong> Runtime Environment and Source Code Audit</p>",
				"</div>",
				
				"<h2>1. Runtime Environment (sessionInfo())</h2>",
				"<p>This section documents the specific R version and package libraries used, essential for reproducibility.</p>",
				"<pre>", htmltools::htmlEscape(runtime_log), "</pre>",
				
				"<h2>2. Complete Application Source Code (app.R)</h2>",
				"<p>This is the full R script code executed by the Shiny application for independent verification.</p>",
				"<pre>", htmltools::htmlEscape(source_code), "</pre>",
				
				"</body>",
				"</html>"
			)
			
			# --- 4. Write HTML ---
			writeLines(report_html, file)
		},
		contentType = "text/html" 
	)
	
	
	output$normality_test_code <- renderText({
		"# Shapiro-Wilk Normality Test Implementation

shapiro_wilk_test <- function(data) {
  # Remove NA values
  data_clean <- na.omit(data)
  n <- length(data_clean)
  
  # Shapiro-Wilk requires 3 \\le n \\le 5000
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
AD_Mean_corrected <- AD_Mean \\times bias_correction_factor
EF_Mean_corrected <- EF_Mean \\times bias_correction_factor"
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
   - UA factor \\in [0, 1]
   - Net ERRs \\le Gross ERRs
   - All deductions \\ge 0

4. ART TREES Compliance Check:
   - 90% confidence interval used
   - 10,000 Monte Carlo iterations
   - Equation 10 & 11 correctly applied
   - 5% buffer deduction included

5. Cross-Validation Integrity:
   - Training set \\ne Test set (no data leakage)
   - Performance metrics calculated on held-out data
   - Bias correction within conservative bounds

6. Reproducibility Check:
   - Set seed for random number generation
   - Document R version and package versions
   - Save intermediate results for audit trail"
	})
}

# Run the enhanced application
shinyApp(ui = ui, server = server)
