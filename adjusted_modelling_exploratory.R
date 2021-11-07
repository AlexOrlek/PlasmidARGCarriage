pacman::p_load(tidyverse, mgcv, data.table)
source('functions.R')
dir.create('output_adjusted', showWarnings = FALSE)
dir.create('output_adjusted/exploratory', showWarnings = FALSE)
set.seed(42)

finaldftrunc<-read.table('data/plasmiddf_transformed.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")
outcomeclasses<-c('aminoglycoside','phenicol','sulphonamide','tetracycline','macrolide','TEM.1','trimethoprim','ESBL', 'carbapenem','quinolone','colistin')

# convert categorical variables to factors
finaldftrunc$BiocideMetalResistance<-as.factor(finaldftrunc$BiocideMetalResistance)
finaldftrunc$Virulence<-as.factor(finaldftrunc$Virulence)
finaldftrunc$Integron<-as.factor(finaldftrunc$Integron)
finaldftrunc$ConjugativeSystem<-factor(finaldftrunc$ConjugativeSystem, ordered = FALSE,levels = c("non-mobilisable", "mobilisable", "conjugative"))
finaldftrunc$GeographicLocation<-factor(finaldftrunc$GeographicLocation, ordered = FALSE,levels = c("high-income", "middle-income", "China", "United States", "EU", "other"))
finaldftrunc$IsolationSource<-factor(finaldftrunc$IsolationSource, ordered = FALSE,levels = c("human", "livestock","other"))
finaldftrunc$RepliconCarriage<-factor(finaldftrunc$RepliconCarriage, ordered = FALSE,levels = c("untyped", "single-replicon", "multi-replicon"))
finaldftrunc$HostTaxonomy<-factor(finaldftrunc$HostTaxonomy, ordered = FALSE,levels = c("Enterobacteriaceae", "Proteobacteria_other", "Firmicutes","other"))
finaldftrunc$HostTaxonomy<-factor(finaldftrunc$HostTaxonomy,levels=c('Enterobacteriaceae','Proteobacteria (non-Enterobacteriaceae)','Firmicutes','other'))
for (outcomeclass in outcomeclasses) {
  finaldftrunc[,gsub('%s',outcomeclass,'outcome%s')]<-as.factor(finaldftrunc[,gsub('%s',outcomeclass,'outcome%s')])
}


# ---------------------------
# initial modelling (test different possible model structures with smallest outcomeclass [colistin] and largest outcome class [aminoglycoside])

# colistin models
outcomeclass <- 'colistin'
# failed model - inadequate df
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize)+s(InsertionSequenceDensity)+s(NumOtherResistanceClasses%s)+s(CollectionDate)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
model1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
# Error in smooth.construct.tp.smooth.spec(object, dk$data, dk$knots) : 
#   A term has fewer unique covariate combinations than specified maximum degrees of freedom

# initial model that works
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses%s,k=5)+s(CollectionDate,k=5)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
model1.1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.1)  # log10PlasmidSize basis dimension inadequate at p < 0.0001; InsertionSequenceDensity inadequate at p < 0.05

# increasing InsertionSequenceDensity and log10PlasmidSize basis dimensions fails to resolve low p-value (see gam.check output)
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize,k=7)+s(InsertionSequenceDensity,k=7)+s(NumOtherResistanceClasses%s,k=5)+s(CollectionDate,k=5)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
model1.2.1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.2.1)

frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize,k=10)+s(InsertionSequenceDensity,k=10)+s(NumOtherResistanceClasses%s,k=5)+s(CollectionDate,k=5)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
model1.2.2<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.2.2)


# aminoglycoside models
outcomeclass<-'aminoglycoside'

# failed model - inadequate df
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize)+s(InsertionSequenceDensity)+s(NumOtherResistanceClasses%s)+s(CollectionDate)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
model1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
# Error in smooth.construct.tp.smooth.spec(object, dk$data, dk$knots) : 
#   A term has fewer unique covariate combinations than specified maximum degrees of freedom

# initial model that works
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses%s,k=5)+s(CollectionDate,k=5)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
model1.1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.1)  # log10PlasmidSize basis dimension inadequate at p < 0.0001; InsertionSequenceDensity inadequate at p < 0.05

# increasing InsertionSequenceDensity and log10PlasmidSize basis dimensions fails to resolve low p-value (see gam.check output)
frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize,k=7)+s(InsertionSequenceDensity,k=7)+s(NumOtherResistanceClasses%s,k=5)+s(CollectionDate,k=5)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
model1.2.1<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.2.1)

frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize,k=10)+s(InsertionSequenceDensity,k=10)+s(NumOtherResistanceClasses%s,k=5)+s(CollectionDate,k=5)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
model1.2.2<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
gam.check(model1.2.2)


# Because log10PlasmidSize basis dimensionality remains inadequate at k=10, I will choose model1.1
summary(model1.1)
anova.gam(model1.1)  # all terms except log10PlasmidSize and CollectionDate are significant at p<0.05; log10PlasmidSize and CollectionDate show effect for other outcomes (not shown) - I will not simplify model structure
# Now, need to perform model checking across all classes...


# ---------------------------
# all class modelling
dfcolnames<-colnames(finaldftrunc)

# baseline model created for all classes
modellist<-list()
for (outcomeclass in outcomeclasses) {
  print(outcomeclass)
  #need to temporarily rename NumOtherResistanceClasses%s column so that all models have same explanatory variable names
  predictorindx<-which(dfcolnames==gsub('%s',outcomeclass,'NumOtherResistanceClasses%s'))
  colnames(finaldftrunc)[predictorindx]<-'NumOtherResistanceClasses'
  frm <- formula(gsub('%s',outcomeclass,'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=5)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'))
  modellist[[outcomeclass]]<-gam(frm,family='binomial',data=finaldftrunc,method = 'ML')
  colnames(finaldftrunc)[predictorindx]<-gsub('%s',outcomeclass,'NumOtherResistanceClasses%s')
}


# gam.check output
gamcheckfile<-'output_adjusted/exploratory/gamcheck_model1.1.txt'
concurvityfile<-'output_adjusted/exploratory/concurvity_model1.1.txt'
cat('gam.check output\n',file=gamcheckfile,sep = '\n',append = FALSE)
cat('concurvity output\n',file=concurvityfile,sep = '\n',append = FALSE)
for (outcomeclass in outcomeclasses) {
  gamcheckplotfile<-gsub('%s',outcomeclass,'output_adjusted/exploratory/gamcheckplots_model1.1_%s.pdf')
  cat(gsub('%s',outcomeclass,'\n###gam.check for outcome: %s binary'),file=gamcheckfile,sep = '\n',append = TRUE)
  writegamcheck(modellist,outcomeclass,gamcheckfile,gamcheckplotfile)
  cat(gsub('%s',outcomeclass,'\n###concurvity for outcome: %s binary'),file=concurvityfile,sep = '\n',append = TRUE)
  writeconcurvity(modellist,outcomeclass,concurvityfile)
}


# N.B gam.check residual checking p-values are stochastic (but overall, inadequacy of log10PlasmidSize basis dimension is supported)
