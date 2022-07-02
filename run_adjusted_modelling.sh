# ---------------------------
# code to run adjusted_modelling.R in parallel with different models
# adjust -j to set number of threads

modelnames=('mainmodel' 'mainmodel_NonImputedCollectionDates' 'minus_NumOtherResistanceClasses' 'minus_Integron' 'minus_BiocideMetalResistance' 'minus_RepliconCarriage' 'minus_HostTaxonomy' 'minus_log10PlasmidSize' 'minus_ConjugativeSystem' 'minus_RepliconCarriage_NumOtherResistanceClasses' 'minus_BiocideMetalResistance_NumOtherResistanceClasses' 'minus_Integron_NumOtherResistanceClasses' 'minus_BiocideMetalResistance_Integron' 'minus_3associatedfactorsofConjugativeSystem' 'minus_6associatedfactorsofConjugativeSystem')  # 15 models in total

printf '%s\n' "${modelnames[@]}" | parallel -j 5 "Rscript adjusted_modelling.R {}"


#NOTES
# minus_6associatedfactorsofConjugativeSystem: minus log10PlasmidSize, NumOtherResistanceClasses, HostTaxonomy, InsertionSequenceDensity, Integron, RepliconCarriage
# minus_3associatedfactorsofConjugativeSystem: minus log10PlasmidSize, NumOtherResistanceClasses, HostTaxonomy
