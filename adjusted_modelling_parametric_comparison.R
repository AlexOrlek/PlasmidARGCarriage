# Comparison of parametric term log-odds ratios across models

# ---------------------------
# compare adjusted vs unadjusted models
unadjusteddf<-read.table('output_unadjusted/unadjustedodds.tsv',as.is = TRUE,sep='\t',header=TRUE)
adjusteddf<-read.table('output_adjusted/mainmodel/adjustedodds.tsv',as.is = TRUE,sep='\t',header=TRUE)

#is.na(as.numeric(unadjusteddf$logOddsRatio[unadjusteddf$logOddsRatio!='baseline']))
#is.na(unadjusteddf$logOddsRatio) | unadjusteddf$logOddsRatio!='baseline'

# remove 'baseline' rows from unadjusted
unadjusteddf<-unadjusteddf[unadjusteddf$logOddsRatio!='baseline' | is.na(unadjusteddf$logOddsRatio),]

# re-order by factor level
unadjusteddf<-unadjusteddf[order(unadjusteddf$FactorVariable),]
adjusteddf<-adjusteddf[order(adjusteddf$FactorVariable),]

# append comparison to adjusteddf
#adjusteddf$ProbabilityScale_unadjusted<-unadjusteddf$ProbabilityScale
#adjusteddf$ProbabilityScale_difference_vs_unadjusted<-as.numeric(adjusteddf$ProbabilityScale)-as.numeric(adjusteddf$ProbabilityScale_unadjusted)
adjusteddf$logOddsRatio_unadjusted<-unadjusteddf$logOddsRatio
adjusteddf$logOddsRatio_difference_vs_unadjusted<-as.numeric(adjusteddf$logOddsRatio)-as.numeric(adjusteddf$logOddsRatio_unadjusted)

# save
write.table(adjusteddf,file='output_adjusted/mainmodel/adjustedodds_vs_unadjusted.tsv',row.names = FALSE,col.names = TRUE,sep='\t')


# ---------------------------
# compare adjusted alternative models vs 1) unadjusted; 2) adjusted (baseline "mainmodel")

alternativemodelnames<-c('mainmodel_minus_NumOtherResistanceClasses', 'mainmodel_minus_Integron', 'mainmodel_minus_BiocideMetalResistance', 'mainmodel_minus_RepliconCarriage', 'mainmodel_minus_HostTaxonomy', 'mainmodel_minus_log10PlasmidSize', 'mainmodel_minus_InsertionSequenceDensity', 'mainmodel_minus_RepliconCarriage_NumOtherResistanceClasses', 'mainmodel_minus_BiocideMetalResistance_NumOtherResistanceClasses', 'mainmodel_minus_Integron_NumOtherResistanceClasses', 'mainmodel_minus_BiocideMetalResistance_Integron', 'mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity', 'mainmodel_minus_associatedfactorsofConjugativeSystem')
removedparametricvariables<-list(NA,c('Integron'),c('BiocideMetalResistance'),c('RepliconCarriage'),c('HostTaxonomy'),NA,NA,c('RepliconCarriage'),c('BiocideMetalResistance'),c('Integron'),c('BiocideMetalResistance','Integron'),NA,c('HostTaxonomy','RepliconCarriage','Integron'))

for (i in 1:length(alternativemodelnames)) {
  modelname<-alternativemodelnames[i]
  removedvars<-removedparametricvariables[[i]]
  print(modelname)
  #read adjusted output from alternative model
  adjusteddf_alternative<-read.table(gsub('%s',modelname,'output_adjusted/%s/adjustedodds.tsv'),as.is = TRUE,sep='\t',header=TRUE)
  #re-order by factor level
  adjusteddf_alternative<-adjusteddf_alternative[order(adjusteddf_alternative$FactorVariable),]
  #edit unadjusted and adjusted mainmodel output to match rows of adjusted alternative (remove rows of removed factor)
  unadjusteddf_edit<-unadjusteddf
  adjusteddf_edit<-adjusteddf
  if (!is.na(removedparametricvariables[i])) {
    unadjusteddf_edit<-unadjusteddf[!unadjusteddf$FactorVariable %in% removedvars,]
    adjusteddf_edit<-adjusteddf[!adjusteddf$FactorVariable %in% removedvars,]
  }
  #append log-odds comparisons vs unadjusted and adjusted mainmodel
  stopifnot(all((adjusteddf_alternative$ResistanceClass==unadjusteddf_edit$ResistanceClass)==TRUE))
  stopifnot(all((adjusteddf_alternative$FactorVariable==unadjusteddf_edit$FactorVariable)==TRUE))
  stopifnot(all((adjusteddf_alternative$FactorLevel==unadjusteddf_edit$FactorLevel)==TRUE))
  stopifnot(all((adjusteddf_alternative$ResistanceClass==adjusteddf_edit$ResistanceClass)==TRUE))
  stopifnot(all((adjusteddf_alternative$FactorVariable==adjusteddf_edit$FactorVariable)==TRUE))
  stopifnot(all((adjusteddf_alternative$FactorLevel==adjusteddf_edit$FactorLevel)==TRUE))
  #adjusteddf_alternative$ProbabilityScale_unadjusted<-unadjusteddf_edit$ProbabilityScale
  #adjusteddf_alternative$ProbabilityScale_mainmodel<-adjusteddf_edit$ProbabilityScale
  #adjusteddf_alternative$ProbabilityScale_difference_vs_unadjusted<-as.numeric(adjusteddf_alternative$ProbabilityScale)-as.numeric(adjusteddf_alternative$ProbabilityScale_unadjusted)
  #adjusteddf_alternative$ProbabilityScale_difference_vs_mainmodel<-as.numeric(adjusteddf_alternative$ProbabilityScale)-as.numeric(adjusteddf_alternative$ProbabilityScale_mainmodel)
  adjusteddf_alternative$logOddsRatio_unadjusted<-unadjusteddf_edit$logOddsRatio
  adjusteddf_alternative$logOddsRatio_mainmodel<-adjusteddf_edit$logOddsRatio
  adjusteddf_alternative$logOddsRatio_difference_vs_unadjusted<-as.numeric(adjusteddf_alternative$logOddsRatio)-as.numeric(adjusteddf_alternative$logOddsRatio_unadjusted)
  adjusteddf_alternative$logOddsRatio_difference_vs_mainmodel<-as.numeric(adjusteddf_alternative$logOddsRatio)-as.numeric(adjusteddf_alternative$logOddsRatio_mainmodel)
  #save
  write.table(adjusteddf_alternative,file=gsub('%s',modelname,'output_adjusted/%s/adjustedodds_vs_mainmodel.tsv'),row.names = FALSE,col.names = TRUE,sep='\t')
}
