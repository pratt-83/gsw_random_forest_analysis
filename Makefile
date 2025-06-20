# Makefile

# path to input data file, expects csv
INPUT = ~/Maize_stoma_proj/Data/head_file_abc_pheno_atmo_2023_ONLY.csv
OUTDIR = Results

# default target, executes when you type 'make' into the terminal
all:
	@mkdir -p $(OUTDIR)
	@for input in $(INPUT); do \
		echo "Begin RF analysis with the following file: $$input ..."; \
		echo "Time: $$(date '+%Y-%m-%d %H:%M:%S')"; \
		Rscript Scripts/01_rf_train.R $$input $(OUTDIR); \
		echo "Finished RF analysis at $$(date '+%Y-%m-%d %H:%M:%S')"; \
	done
