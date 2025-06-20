#!/usr/bin/env Rscript

########## About
# Author: Nicholas Pratt
# Usage : Plot results of RF analysis from importance and accuracy files.

########## Load Libraries
# Load libraries for plotting and data manipulation, suppressing startup messages for clarity.
suppressPackageStartupMessages({
	library(ggplot2)	# For plotting
	library(dplyr)		# For data wrangling
	library(tidyr)		# For data reshaping
})

########## Parse Command Line Arguments
# Expecting: <importance_file> <accuracy_file> <trait> <allele> <year>
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) {
	stop("Usage: Rscript 02_rf_plot.R <importance_file> <accuracy_file> <trait> <allele> <year>")
}
importance_file <- args[1]	# Path to combined_importance_scores_<trait>.txt
accuracy_file   <- args[2]	# Path to accuracy_metrics_<trait>.txt
trait           <- args[3]	# Trait name (a, b, or c)
allele          <- args[4]	# Allele type (WT or mu)
year            <- args[5]	# Year (e.g., 2023)

########## Load Data
# Read variable importance and accuracy metrics data
imp_data <- read.table(importance_file, header = TRUE, sep = "\t", stringsAsFactors = FALSE)
acc_data <- read.table(accuracy_file,   header = TRUE, sep = "\t", stringsAsFactors = FALSE)

########## Prepare Output Directory
# All plots will be saved in the same directory as the input files
out_dir <- dirname(importance_file)

########## Tidy importance into long form
# Convert wide format (Fold_1, Fold_2, ...) to long format for plotting
imp_long <- imp_data %>%
	pivot_longer(
		cols = starts_with("Fold_"),
		names_to  = "Fold",
		values_to = "Importance"
	)

# Order variables by their overall median importance (low to high)
var_order <- imp_long %>%
	group_by(Variable) %>%
	summarize(med = median(Importance, na.rm = TRUE)) %>%
	arrange(med) %>%
	pull(Variable)
imp_long$Variable <- factor(imp_long$Variable, levels = var_order)

########## 1) Accuracy plots (RMSE + R²)
# Plot both RMSE and R² (one page each) for all models
pdf(file.path(out_dir, paste0("accuracy_plots_", trait, ".pdf")), height = 5, width = 10)

# Plot RMSE for each model
print(
	ggplot(acc_data, aes(x = Model, y = RMSE, group = 1)) +
		geom_point() + geom_line() +
		theme_light() +
		labs(x = "Model", y = "RMSE") +
		theme(
			axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
			axis.text.y = element_text(size = 16),
			axis.title.x = element_text(size = 18, face = "bold"),
			axis.title.y = element_text(size = 18, face = "bold")
		)
)

# Plot R² for each model
print(
	ggplot(acc_data, aes(x = Model, y = R2, group = 1)) +
		geom_point() + geom_line() +
		theme_light() +
		labs(x = "Model", y = expression(R^2)) +
		theme(
			axis.text.x = element_text(angle = 90, hjust = 1, size = 16),
			axis.text.y = element_text(size = 16),
			axis.title.x = element_text(size = 18, face = "bold"),
			axis.title.y = element_text(size = 18, face = "bold")
		)
)
dev.off()

########## 2) Combined importance boxplot
# Create boxplot summarizing variable importances across all models and folds
pdf(file.path(out_dir, paste0("importance_boxplot_", trait, ".pdf")), height = 5, width = 10)
print(
	ggplot(imp_long, aes(x = Variable, y = Importance)) +
		geom_boxplot() +
		theme_light() +
		labs(x = "Environmental Parameters", y = "Importance Score") +
		theme(
			axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 16),
			axis.text.y = element_text(size = 16),
			axis.title.x = element_text(size = 18, face = "bold"),
			axis.title.y = element_text(size = 18, face = "bold"),
			text = element_text(size = 16),
			legend.text.align = 0,
			legend.text = element_text(size = rel(0.8))
		)
)
dev.off()

########## 3) One importance plot per model
# For each unique model, create a separate variable importance boxplot PDF
for (mod in unique(imp_data$model)) {
	sub <- filter(imp_long, model == mod)
	pdf(file.path(out_dir, paste0("Importance_", mod, "_cforest.pdf")), height = 5, width = 10, family = "sans")
	print(
		ggplot(sub, aes(x = Variable, y = Importance)) +
			geom_boxplot() +
			theme_light() +
			labs(x = "Environmental Parameters", y = "Importance Score") +
			theme(
				strip.background = element_blank(),
				panel.grid.major = element_blank(),
				panel.grid.minor = element_blank(),
				panel.background = element_blank(),
				strip.text.y = element_text(size = 8, colour = 'black'),
				strip.text.x = element_text(size = 8, colour = 'black'),
				strip.placement = "outside",
				axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0, size = 16),
				axis.text.y = element_text(size = 16),
				axis.title.x = element_text(size = 18, face = "bold"),
				axis.title.y = element_text(size = 18, face = "bold"),
				text = element_text(size = 16),
				legend.text.align = 0,
				legend.text = element_text(size = rel(0.8))
			)
	)
	dev.off()
}

# Print message to user when done
message("Finished plotting RF results for trait=", trait,
	", allele=", allele, ", year=", year)

