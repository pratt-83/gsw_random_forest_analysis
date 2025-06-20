# Maize Stomatal Conductance Random Forest Pipeline

This repository implements a reproducible pipeline for analyzing Gaussian parameters (`a`, `b`, `c`) of stomatal conductance in maize using Random Forest models.  
It supports parallel cross-validation and variable importance analysis, with clean output and publication-ready plots.

---

## Workflow Overview

The analysis workflow is driven by the Makefile and proceeds in three main steps:

1. **Makefile**  
	- Orchestrates the full pipeline.
	- For each trait (`a`, `b`, `c`) and each allele (`WT`, `mu`), calls `01_rf_train.R` with the appropriate inputs.

2. **01_rf_train.R**  
	- Main training and evaluation script.
	- For each combination of trait and allele:
		- **Data Subsetting:** Selects data for the given year and allele type.
		- **Cross-Validation Setup:** Assigns random 5-fold CV labels to each row.
		- **Scaling:** Scales all predictor (feature) columns.
		- **Model Training (Parallel):**  
			- For each combination of `mtry` and `ntree` hyperparameters:
				- Runs 5-fold cross-validation:
					- For each fold:
						- Trains a cforest (Random Forest) model on k-1 folds, predicts on the held-out fold.
						- Stores predictions and computes variable importance.
				- Aggregates predictions, importance, and accuracy metrics for each model.
		- **Output:**  
			- Writes out all predictions, variable importance scores, and accuracy metrics.
		- **Plotting:**  
			- Calls `02_rf_plot.R` to generate summary plots for each trait/allele/year.

3. **02_rf_plot.R**  
	- Reads variable importance and accuracy files produced by `01_rf_train.R`.
	- Generates the following publication-ready plots (PDF):
		- RMSE and R² accuracy plots for all models.
		- Combined variable importance boxplot for all models.
		- Individual importance plots for each model/hyperparameter combination.
	- All plots are saved alongside the results for traceability.

---

## How to Run

1. **Edit the Makefile** as needed for your dataset locations and desired outputs.
2. Run the full analysis pipeline from the project root:
	```
	make
	```
	- This will process all trait/allele combinations, generate results, and create all plots in the results directory.


