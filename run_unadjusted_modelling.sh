# ---------------------------
# code to run unadjusted_modelling.R in parallel with different models
# adjust -j to set number of threads

modelnames=('log10PlasmidSize' 'InsertionSequenceDensity' 'NumOtherResistanceClasses' 'CollectionDate' 'CollectionDate_NonImputedCollectionDates')
printf '%s\n' "${modelnames[@]}" | parallel -j 5 "Rscript unadjusted_modelling.R {}"
