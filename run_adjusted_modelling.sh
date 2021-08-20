# ---------------------------
# code to run adjusted_modelling.R in parallel with different models
# adjust -j to set number of threads

modelnames=('mainmodel' 'mainmodel_minus_NumOtherResistanceClasses' 'mainmodel_minus_Integron' 'mainmodel_minus_BiocideMetalResistance' 'mainmodel_minus_RepliconCarriage' 'mainmodel_minus_HostTaxonomy' 'mainmodel_minus_log10PlasmidSize' 'mainmodel_minus_InsertionSequenceDensity' 'mainmodel_minus_RepliconCarriage_NumOtherResistanceClasses' 'mainmodel_minus_BiocideMetalResistance_NumOtherResistanceClasses' 'mainmodel_minus_Integron_NumOtherResistanceClasses' 'mainmodel_minus_BiocideMetalResistance_Integron' 'mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity' 'mainmodel_minus_associatedfactorsofConjugativeSystem')  # 14 models in total

printf '%s\n' "${modelnames[@]}" | parallel -j 5 "Rscript adjusted_modelling.R {}"
