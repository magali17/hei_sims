# script runs R scripts to: 1) estimate annual averages, 2) make UK-pls predictions, 3) evaluate models, 4) generate summary figures and tables
# to run this script in the terminal, enter: bash run_scripts.bash

## run rscripts
rscript 1_annual_avg.R
rscript 2.0_uk_workspace.R
rscript 2_uk_cv.R
rscript 2_uk_test.R
rscript 2_uk_grid.R
rscript 2_uk_combine_predictions.R
rscript 3_model_eval.R

## knit markdown with results 
Rscript -e 'rmarkdown::render("4_results_summary.Rmd", "html_document")'

echo "done running scripts"
