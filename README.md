# Background

Code for paper: "Risk factors associated with plasmid antibiotic resistance gene carriage revealed using large-scale multivariable analysis".<br>

Code runs statistical analysis of associated factors of plasmid antibiotic resistance gene (ARG) carriage: 1) exploratory analysis; 2) unadjusted analysis; 3) adjusted analysis.<br>

Separate logistic GAM models are constructed for 10 major ARG type outcomes.<br>
<br>

# Instructions for running code

```
# exploratory analysis: data wrangling and transformation prior to modelling; investigating associations between variables; descriptive analysis of ARGs

Rscript exploratory_analysis.R
Rscript exploratory_analysis_characteristics.R 
Rscript exploratory_analysis_confounding.R
Rscript arg_descriptive_analysis.R

# unadjusted analysis: calculate and plot unadjusted odds ratios

Rscript unadjusted_oddsratiocalc.R
Rscript unadjusted_plotodds.R
bash run_unadjusted_modelling.sh  # runs unadjusted_modelling.R; specify nthreads using -j argument

# adjusted analysis: construct and visualise multivariate GAM models using the mgcv package; compare output from different models

Rscript adjusted_modelling_exploratory.R
bash run_adjusted_modelling.sh  # runs adjusted_modelling.R; specify nthreads using -j argument; specify models to run in modelnames argument.
Rscript adjusted_modelling_parametric_comparison.R
Rscript adjusted_modelling_smooth_comparison.R

```

# Licence
[MIT License](https://en.wikipedia.org/wiki/MIT_License)