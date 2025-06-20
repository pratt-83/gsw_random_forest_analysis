#!/usr/bin/env Rscript

########## About
# Author: Nicholas Pratt
# Usage : Train cforest RF models across Gaussian parameters and environmental factors for maize phenotyping.
#         Assumes data files are already subset by year; each run only models a single year.

########## Documentation
# party documentation: https://cran.r-project.org/web/packages/party/party.pdf

########## Preprocessing

# Load libraries for parallel processing, random forest, and data manipulation
suppressPackageStartupMessages({
	library(parallel)		# Enables multicore computations (detectCores, makeCluster, etc.)
	library(party)			# Provides cforest for conditional inference random forest
	library(doParallel)		# Backend for registering parallel computation for foreach
	library(foreach)		# Provides %dopar% for parallel loops
	library(tidyr)			# For tidy reshaping of data (e.g., pivot_longer)
	library(dplyr)			# For general data wrangling
})

# Read command line arguments: input CSV file path and output directory path
args <- commandArgs(trailingOnly = TRUE)
data_csv_file <- args[1]	# First argument: CSV file to process
out_dir <- args[2]			# Second argument: Where to save results

# Define the Gaussian traits to model and the two allele types to analyze
traits <- c("a", "b", "c")
alleles <- c("WT", "mu")

# Load the input CSV data (already subset by year) into a data frame
data_file <- read.csv(data_csv_file)
data_scaled <- data_file	# (Optionally: could scale here, but we're scaling per trait/allele below)

# Extract year from the "Year" column (assume all rows should have same year)
data_year <- data_file$Year[1]

# Warn user if the input file is not truly subset by a single year
if (!all(data_file$Year == data_year)) {
	warning(
		paste(
			"Mismatch in Year column: not all rows are from year:", data_year,
			"\nPlease verify that your input CSV file is subset by year"
		)
	)
}

########## General 5-Fold Cross-Validation Function
# Performs cross-validation using a user-supplied model training function (e.g., cforest_model)
cv_function <- function(data, trait, predictors, k = 5, model_function, ...) {
	# Allocate storage for all predicted values for each fold
	all_predictions <- data.frame(
		Trait = data[[trait]],		# True trait values
		cv_fold = data$cv_fold,		# Assigned fold for each sample
		predicted = NA				# Will be filled in with predicted values
	)
	# Allocate variable importance table: rows = predictors, columns = Fold_1, Fold_2, ..., model
	importance_wide <- data.frame(
		Variable = predictors,
		stringsAsFactors = FALSE
	)
	# Loop over each fold: train on k-1 folds, test on the kth fold
	for (fold in 1:k) {
		train_data <- subset(data, cv_fold != fold)	# Training subset
		test_data  <- subset(data, cv_fold == fold)	# Test subset for this fold

		# Train the model using the provided function (here, cforest_model)
		model <- model_function(train_data, trait, predictors, ...)

		# Predict trait for test data using only predictor columns
		predictions <- predict(model, newdata = test_data[, predictors], type = "response")
		all_predictions$predicted[all_predictions$cv_fold == fold] <- predictions

		# Compute variable importance for this fold
		vi <- varimp(model)	# Named vector of variable importances
		importance_wide[[paste0("Fold_", fold)]] <- as.numeric(vi[predictors])
	}
	# Return all cross-validated predictions and per-fold importance
	return(list(predictions = all_predictions, importance = importance_wide))
}

# Wrapper function to build a cforest model with given parameters
cforest_model <- function(train_data, trait, predictors, ntree, mtry) {
	# Build R formula for the model: Trait ~ predictor1 + predictor2 + ...
	formula <- as.formula(paste(trait, "~", paste(predictors, collapse = "+")))
	# Train the model with unbiased splitting
	model <- cforest(
		formula,
		data = train_data,
		controls = cforest_unbiased(ntree = ntree, mtry = mtry)
	)
	return(model)
}

########## Run Analysis Function 
# Trains and evaluates models for one trait, allele, and year; writes all outputs and calls plotting script.
run_analysis <- function(data, trait, year, allele_type, all_traits, out_directory) {
	# Subset for current year and allele type (row filter)
	data_filtered <- subset(data, Year == year & grepl(allele_type, Allele))
	# If there's no data, print a message and skip this trait/allele
	if (nrow(data_filtered) == 0) {
		message("No data available for allele type ", allele_type, " and trait ", trait)
		return(NULL)
	}

	# Assign folds randomly for 5-fold cross-validation
	set.seed(999)	# For reproducibility
	data_filtered$cv_fold <- sample(rep(1:5, length.out = nrow(data_filtered)))
	y <- data_filtered[[trait]]

	# Drop columns that should NOT be used as predictors; keep trait + fold + predictors
	drop <- c(
		"Hybrid", "Geno", "Allele", "Year", "Date", "PD", "DAP", "GDD",
		"Cumulative_GDD", "days_to_anthesis", "days_to_silking", 
		"Wind_Direction_max", "Wind_Direction_min", "Wind_Direction_var", "Wind_Direction_mean",
		"Leaf_Wetness_max", "Leaf_Wetness_min", "Leaf_Wetness_var", "Leaf_Wetness_mean"
	)
	# Remove drop columns, but keep trait and fold columns
	data_model <- data_filtered[, !(names(data_filtered) %in% drop)]

	# Identify predictor columns: exclude all other traits and the fold column itself
	predictor_cols <- setdiff(names(data_model), c(all_traits, "cv_fold"))

	# SAFETY: Print and check predictors; fail if none remain
	#cat("Predictor columns after dropping and exclusions for trait", trait, ":", predictor_cols, "\n")
	if (length(predictor_cols) == 0) {
		stop("No predictor columns remain after filtering; check your input data and drop list.")
	}

	# Scale all predictor columns (mean 0, sd 1) but leave response and fold untouched
	data_model_scaled <- data_model
	data_model_scaled[, predictor_cols] <- scale(data_model[, predictor_cols])

	# Ensure 'Trait' and 'cv_fold' columns are correct (overwrite with originals)
	data_model_scaled$Trait <- data_model[[trait]]
	data_model_scaled$cv_fold <- data_filtered$cv_fold

	# Define grid of model hyperparameters to search (mtry: number of variables at split, ntree: number of trees)
	mtry_vals <- c(10, 20, 30, 40, 50, 60)
	ntree_vals <- c(500, 1000, 2000)
	param_grid <- expand.grid(mtry = mtry_vals, ntree = ntree_vals)

	# Initialize parallel backend using all available cores minus one (keep system responsive)
	cl <- makeCluster(max(1, detectCores() - 1))
	registerDoParallel(cl)

	# For every parameter combo, train/test in parallel
	results <- foreach(
		i = seq_len(nrow(param_grid)),
		.packages = "party",
		.export   = c("cv_function", "cforest_model")
	) %dopar% {
		this_mtry  <- param_grid$mtry[i]
		this_ntree <- param_grid$ntree[i]
		res <- cv_function(
			data = data_model_scaled,
			trait = "Trait",
			predictors = predictor_cols,
			k = 5,
			model_function = cforest_model,
			ntree = this_ntree,
			mtry  = this_mtry
		)
		# Annotate model parameters for later identification
		model_name <- paste0("mtry_", this_mtry, "_ntree_", this_ntree)
		res$predictions$model  <- model_name
		res$importance$model   <- model_name
		res
	}
	# Stop the parallel backend (clean up cluster)
	stopCluster(cl)

	# Create output directory structure for this trait and allele (safe if already exists)
	output_path <- file.path(out_directory, trait, tolower(allele_type))
	dir.create(output_path, recursive = TRUE, showWarnings = FALSE)

	# 1) Write out combined predictions from all models/folds
	combined_preds <- do.call(rbind, lapply(results, function(x) x$predictions))
	write.table(
		combined_preds,
		file.path(output_path, paste0("combined_predictions_", trait, ".txt")),
		row.names = FALSE, sep = "\t", quote = FALSE
	)

	# 2) Write out combined variable importance (wide, per-model)
	combined_importance <- do.call(rbind, lapply(results, function(x) x$importance))
	write.table(
		combined_importance,
		file.path(output_path, paste0("combined_importance_scores_", trait, ".txt")),
		row.names = FALSE, sep = "\t", quote = FALSE
	)

	# 3) Compute and write out accuracy metrics (RMSE, R2) for each model
	accuracy <- do.call(rbind, lapply(results, function(x) {
		preds <- x$predictions
		rmse <- sqrt(mean((preds$Trait - preds$predicted)^2, na.rm = TRUE))
		r2   <- cor(preds$Trait, preds$predicted, use = "complete.obs")^2
		data.frame(Model = unique(preds$model), RMSE = rmse, R2 = r2)
	}))
	write.table(
		accuracy,
		file.path(output_path, paste0("accuracy_metrics_", trait, ".txt")),
		row.names = FALSE, sep = "\t", quote = FALSE
	)

	# 4) Log the best model (highest R2) to a human-readable text file
	best_model <- accuracy[which.max(accuracy$R2), ]
	log_file <- file.path(output_path, "best_model_log.txt")
	writeLines(
		c(
			paste("Best Model for Trait:", trait),
			paste("Allele Type:", allele_type),
			paste("Year:", year),
			paste("Model Parameters:", best_model$Model),
			paste("RMSE:", best_model$RMSE),
			paste("R2:", best_model$R2)
		),
		con = log_file
	)

	######### Plot Results 
	# Calls plotting script with result file paths as arguments (for this trait/allele/year)
	importance_path <- file.path(output_path, paste0("combined_importance_scores_", trait, ".txt"))
	accuracy_path <- file.path(output_path, paste0("accuracy_metrics_", trait, ".txt"))
	cmd <- paste(
		"Rscript Scripts/02_rf_plot.R",
		shQuote(importance_path),
		shQuote(accuracy_path),
		shQuote(trait),
		shQuote(allele_type),
		shQuote(year)
	)
	cat("Running plot command:\n", cmd, "\n")
	system(cmd)
}

########## Loop through each trait and allele type (main script logic)
# For each trait ('a', 'b', 'c') and each allele ('WT', 'mu'), run the full pipeline
for (trait in traits) {
	for (allele in alleles) {
		run_analysis(data_scaled, trait, data_year, allele, traits, out_dir)
	}
}

