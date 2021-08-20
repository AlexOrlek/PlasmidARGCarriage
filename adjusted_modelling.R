pacman::p_load(mgcv, mgcv.helper, RColorBrewer, mgcViz, scales, grid, cowplot, devtools, gsubfn, stringr)
set.seed(42)

finaldftrunc<-read.table('data/plasmiddf_transformed.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")

# ---------------------------
# set parameters
outcomeclasses<-c('aminoglycoside','betalactam_carbapenem','betalactam_ESBL','betalactam_other','macrolide','phenicol','quinolone','sulphonamide','tetracycline','trimethoprim')

args = commandArgs(trailingOnly=TRUE)
modelname<-args[1]
#modelname<-'mainmodel'
#modelname<-'mainmodel_minus_NumOtherResistanceClasses'
#modelname<-'mainmodel_minus_Integron'
#modelname<-'mainmodel_minus_BiocideMetalResistance'
#modelname<-'mainmodel_minus_RepliconCarriage'
#modelname<-'mainmodel_minus_HostTaxonomy'

# create output directories
dir.create(file.path('output_adjusted',modelname), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_adjusted/%s'), 'coefficientplots'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_adjusted/%s'), 'pvalues'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_adjusted/%s/coefficientplots'), 'logoddsscale'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_adjusted/%s/coefficientplots'), 'probscale'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_adjusted/%s/coefficientplots'), 'plotgam'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_adjusted/%s/coefficientplots/plotgam'), 'logoddsscale'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_adjusted/%s/coefficientplots/plotgam'), 'probscale'), showWarnings = FALSE)

if (modelname=='mainmodel') {
  
  #main model
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
  
  #models with single terms removed
} else if (modelname=='mainmodel_minus_NumOtherResistanceClasses') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(CollectionDate,k=3)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_Integron') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_BiocideMetalResistance') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate','Integron','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+Integron+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_RepliconCarriage') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+Integron+BiocideMetalResistance+ConjugativeSystem+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_HostTaxonomy') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_log10PlasmidSize') {
  outputnames_subset<-c('InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_InsertionSequenceDensity') {
  outputnames_subset<-c('log10PlasmidSize','NumOtherResistanceClasses','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
  
  #models with multiple terms removed
} else if (modelname=='mainmodel_minus_RepliconCarriage_NumOtherResistanceClasses') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(CollectionDate,k=3)+Integron+BiocideMetalResistance+ConjugativeSystem+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_BiocideMetalResistance_NumOtherResistanceClasses') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','CollectionDate','Integron','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(CollectionDate,k=3)+Integron+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_Integron_NumOtherResistanceClasses') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','CollectionDate','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(CollectionDate,k=3)+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_BiocideMetalResistance_Integron') {
  outputnames_subset<-c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5)+s(InsertionSequenceDensity,k=5)+s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity') {
  outputnames_subset<-c('NumOtherResistanceClasses','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(NumOtherResistanceClasses,k=5)+s(CollectionDate,k=3)+Integron+BiocideMetalResistance+ConjugativeSystem+RepliconCarriage+HostTaxonomy+Virulence+GeographicLocation+IsolationSource'
} else if (modelname=='mainmodel_minus_associatedfactorsofConjugativeSystem') {
  outputnames_subset<-c('CollectionDate','BiocideMetalResistance','ConjugativeSystem','Virulence','GeographicLocation','IsolationSource')
  frmtext<-'outcome%s~s(CollectionDate,k=3)+BiocideMetalResistance+ConjugativeSystem+Virulence+GeographicLocation+IsolationSource'
} else {
  stop('modelname not recognised')
}

brewerpal1<-brewer.pal(8,'Set1')
brewerpal1[6]<-'#CDCD00'
#plot(1:8,1:8, col=brewerpal1, pch=16,cex=5)
betalactampal<-brewer.pal(9,'Blues')[c(3,6,9)]
#plot(1:3,1:3,col=betalactampal,pch=16,cex=5)
brewerpal2<-c(brewerpal1[1],betalactampal,brewerpal1[3:length(brewerpal1)])
#plot(1:10,1:10, col=brewerpal2, pch=16,cex=5)
numbetalactam<-sum(grepl('^betalactam',outcomeclasses))

if (numbetalactam==3) {
  brewerpal<-brewerpal2
  stopifnot(length(outcomeclasses)==10)
} else if (numbetalactam==1) {
  brewerpal<-brewerpal1
  stopifnot(length(outcomeclasses)==8)
} else {
  stop('Error: unexpected number of beta-lactam resistance subclasses ***')
}


# ---------------------------
# convert categorical variables to factors
finaldftrunc$BiocideMetalResistance<-as.factor(finaldftrunc$BiocideMetalResistance)
finaldftrunc$Virulence<-as.factor(finaldftrunc$Virulence)
finaldftrunc$Integron<-as.factor(finaldftrunc$Integron)
finaldftrunc$ConjugativeSystem<-factor(finaldftrunc$ConjugativeSystem, ordered = FALSE,levels = c("non-mobilisable", "mobilisable", "conjugative"))
finaldftrunc$GeographicLocation<-factor(finaldftrunc$GeographicLocation, ordered = FALSE,levels = c("high-income", "China", "United States", "other", "middle-income", "EU"))
finaldftrunc$IsolationSource<-factor(finaldftrunc$IsolationSource, ordered = FALSE,levels = c("human", "livestock","other"))
finaldftrunc$RepliconCarriage<-factor(finaldftrunc$RepliconCarriage, ordered = FALSE,levels = c("untyped", "single-replicon", "multi-replicon"))
finaldftrunc$HostTaxonomy<-factor(finaldftrunc$HostTaxonomy, ordered = FALSE,levels = c("Enterobacteriaceae", "Proteobacteria_other", "Firmicutes","other"))
for (outcomeclass in outcomeclasses) {
  finaldftrunc[,gsub('%s',outcomeclass,'outcome%s')]<-as.factor(finaldftrunc[,gsub('%s',outcomeclass,'outcome%s')])
}


# ---------------------------
# MODELLING/PLOTTING USING MGCV PLOT.GAM

# create models
dfcolnames<-colnames(finaldftrunc)
modellist<-list()
for (outcomeclass in outcomeclasses) {
  print(outcomeclass)
  #need to temporarily rename NumOtherResistanceClasses%s column so that all models have same predictor names (so mgcviz works)
  predictorindx<-which(dfcolnames==gsub('%s',outcomeclass,'NumOtherResistanceClasses%s'))
  colnames(finaldftrunc)[predictorindx]<-'NumOtherResistanceClasses'
  frm <- formula(gsub('%s',outcomeclass,frmtext))
  modellist[[outcomeclass]]<-gam(frm,family='binomial',data=finaldftrunc,method = 'REML',gamma=1.5)
  colnames(finaldftrunc)[predictorindx]<-gsub('%s',outcomeclass,'NumOtherResistanceClasses%s')
}

saveRDS(modellist,file = gsub('%s',modelname,'output_adjusted/%s/modellist.rds'))

# plot.gam plots
for (outcomeclass in outcomeclasses) {
  pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outcomeclass),'output_adjusted/%1/coefficientplots/plotgam/logoddsscale/%2.pdf'))
  plot(modellist[[outcomeclass]], pages = 1, all.terms = TRUE,residuals = FALSE,shade=TRUE)
  dev.off()
  pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outcomeclass),'output_adjusted/%1/coefficientplots/plotgam/probscale/%2.pdf'))
  plot(modellist[[outcomeclass]], pages = 1, trans = plogis, shift = coef(modellist[[outcomeclass]])[1], all.terms = TRUE,residuals = FALSE,shade=TRUE)
  dev.off()
}

# save model summaries to file
pvaluesfile<-gsub('%s',modelname,'output_adjusted/%s/pvalues/pvalues.txt')
cat('model anova() and summary() p-values',file=pvaluesfile,sep = '\n',append = FALSE)

for (outcomeclass in outcomeclasses) {
  print(outcomeclass)
  cat(gsub('%s',outcomeclass,'\n###p-values for outcome: %s binary'),file=pvaluesfile,sep = '\n',append = TRUE)
  cat('\n#anova() p-values',file=pvaluesfile,sep = '\n',append = TRUE)
  anovaouts<-capture.output(anova(modellist[[outcomeclass]]))
  replace=FALSE
  for (anovaout in anovaouts) {
    if (replace==TRUE) {
      if (startsWith(anovaout,'Approximate')) {
        cat(anovaout,file=pvaluesfile,sep='\n',append=TRUE)
      } else {
        anovaout<-gsub("(?<!\\<)\\s+", "\t", anovaout,perl = TRUE)
        cat(anovaout,file=pvaluesfile,sep='\n',append=TRUE)
      }
    } else {
      cat(anovaout,file=pvaluesfile,sep='\n',append=TRUE)
    }
    if (startsWith(anovaout,'Parametric')) {
      replace=TRUE
    }
  }
  cat('\n#summary() p-values',file=pvaluesfile,sep = '\n',append = TRUE)
  summaryouts<-capture.output(summary(modellist[[outcomeclass]]))
  replace=FALSE
  for (summaryout in summaryouts) {
    if ('Pr(>|z|)' %in% unlist(strsplit(summaryout,'\\s'))) {
      cat('\tEstimate\tStd. Error\tz value\tPr(>|z|)',file=pvaluesfile,sep='\n',append=TRUE)
    } else {
      if (startsWith(summaryout,'(Intercept)')) {
        replace=TRUE
      }
      if (replace==TRUE) {
        summaryout<-gsub("(?<![United|\\<])\\s+", "\t", summaryout,perl = TRUE)
        cat(summaryout,file=pvaluesfile,sep='\n',append=TRUE)
      } else {
        cat(summaryout,file=pvaluesfile,sep='\n',append=TRUE)
      }
      if (startsWith(summaryout,'Approximate')) {
        replace=TRUE
      }
      if (startsWith(summaryout,'---')) {
        replace=FALSE
      }
    }
  }
}

# check model intercept values and transform to probability; write coefficients and CIs to file
interceptslist<-list()
for (outcomeclass in outcomeclasses) {
  interceptlogodds<-as.numeric(coef(modellist[[outcomeclass]])[1])
  interceptprob<-plogis(interceptlogodds)
  interceptslist[[outcomeclass]]<-c(outcomeclass,interceptlogodds,interceptprob)
}
interceptsdf<-as.data.frame(do.call(rbind,interceptslist))
colnames(interceptsdf)<-c('ResistanceClass','InterceptlogOddsScale','InterceptProbabilityScale')
write.table(interceptsdf,file=gsub('%s',modelname,'output_adjusted/%s/pvalues/intercepts.txt'),col.names = TRUE,row.names = FALSE,sep='\t')

# calculate 95% CIs and save coefficients/CIs on logodds/odds/probability scale
flevels<-lapply(modellist[[outcomeclasses[1]]]$xlevels,function(x) tail(x,-1))
numflevels<-lapply(flevels,function(x) length(x))
flevels<-unname(unlist(flevels))
fnames<-vector()
for (i in 1:length(names(modellist[[outcomeclasses[1]]]$xlevels))) {
  fnames<-c(fnames,rep(names(modellist[[outcomeclasses[1]]]$xlevels)[i],numflevels[i]))
}
fcols<-cbind(fnames,flevels)
colnames(fcols)<-c('FactorVariable','FactorLevel')

coefslist<-list()
for (outcomeclass in outcomeclasses) {
  modelout<-modellist[[outcomeclass]]
  modelsummary<-summary(modelout)
  #get pvalues for parametric terms
  `p-value`<-as.vector(modelsummary$p.table[-1,'Pr(>|z|)'])
  #get odds, logodds, intercept-centred probabilities for parametric terms
  out<-as.data.frame(mgcv.helper::confint.gam(modelout))
  intercept<-out[1,2]
  #logodds scale
  out<-out[-1,]
  out$term<-NULL
  out$Statistic<-NULL
  out$`Std. Error`<-NULL
  colnames(out)<-c('logOddsRatio','logLower95CI','logUpper95CI')
  #odds scale
  outexp<-exp(out)
  colnames(outexp)<-c('OddsRatio','Lower95CI','Upper95CI')
  #probability scale
  outplogis<-out
  outplogis<-lapply(outplogis, function(x) plogis(x))
  outplogis<-do.call(cbind,outplogis)
  colnames(outplogis)<-c('ProbabilityScale','ProbabilityScaleLower95CI','ProbabilityScaleUpper95CI')
  #probability scale, intercept centred
  outplogis_int_centred<-out+intercept
  outplogis_int_centred<-lapply(outplogis_int_centred, function(x) plogis(x))
  outplogis_int_centred<-do.call(cbind,outplogis_int_centred)
  colnames(outplogis_int_centred)<-c('ProbabilityScaleInterceptCentred','ProbabilityScaleInterceptCentredLower95CI','ProbabilityScaleInterceptCentredUpper95CI')
  #combine
  ResistanceClass<-rep(outcomeclass,nrow(fcols))
  dfout<-cbind(ResistanceClass,fcols,outexp,out,outplogis,outplogis_int_centred,`p-value`)
  coefslist[[outcomeclass]]<-dfout
}

coefsdf<-do.call(rbind,coefslist)

# re-order
coefsdf<-coefsdf[order(coefsdf$ResistanceClass,coefsdf$FactorVariable),]
if ('GeographicLocation' %in% outputnames_subset) {
  geographiesdf<-coefsdf[coefsdf$FactorVariable=='GeographicLocation',]
  geographiesdf<-do.call(rbind,lapply(split(geographiesdf,geographiesdf$ResistanceClass),function(x) x[match(c("middle-income","EU","China","United States","other"),x$FactorLevel),]))
  coefsdf[coefsdf$FactorVariable=='GeographicLocation',]<-geographiesdf
}

# in cases of complete separation, mask factor level with NAs
coefsdf[coefsdf$Lower95CI==Inf | coefsdf$Upper95CI==Inf,!colnames(coefsdf) %in% c('ResistanceClass','FactorVariable','FactorLevel')]<-NA

# save
write.table(coefsdf,file=gsub('%s',modelname,'output_adjusted/%s/adjustedodds.tsv'),row.names = FALSE,col.names = TRUE,sep='\t')


# ---------------------------
# MODELLING/PLOTTING USING MGCVIZ

# re-label factor levels (alphabetical prefixes to coerce plotting layout) and re-create/re-plot models using mgcViz
levels(finaldftrunc$BiocideMetalResistance)<-c("Aabsence","Bpresence")
levels(finaldftrunc$Virulence)<-c("Aabsence","Bpresence")
levels(finaldftrunc$Integron)<-c("Aabsence","Bpresence")
levels(finaldftrunc$ConjugativeSystem)<-c("Anon-mobilisable", "Bmobilisable", "Cconjugative")
levels(finaldftrunc$GeographicLocation)<-c("Ahigh-income", "BChina", "CUnited States","Dother","Emiddle-income", "FEU")
levels(finaldftrunc$IsolationSource)<-c("Ahuman", "Blivestock","Cother")
levels(finaldftrunc$RepliconCarriage)<-c("Auntyped", "Bsingle-replicon", "Cmulti-replicon")
levels(finaldftrunc$HostTaxonomy)<-c("AEnterobacteriaceae", "BProteobacteria_other", "CFirmicutes","Dother")


# create models
dfcolnames<-colnames(finaldftrunc)
modellist<-list()
for (outcomeclass in outcomeclasses) {
  print(outcomeclass)
  #need to temporarily rename NumOtherResistanceClasses%s column so that all models have same predictor names (so mgcviz works)
  predictorindx<-which(dfcolnames==gsub('%s',outcomeclass,'NumOtherResistanceClasses%s'))
  colnames(finaldftrunc)[predictorindx]<-'NumOtherResistanceClasses'
  frm <- formula(gsub('%s',outcomeclass,frmtext))
  modellist[[outcomeclass]]<-gam(frm,family='binomial',data=finaldftrunc,method = 'REML',gamma=1.5)
  colnames(finaldftrunc)[predictorindx]<-gsub('%s',outcomeclass,'NumOtherResistanceClasses%s')
}
saveRDS(modellist,file = gsub('%s',modelname,'output_adjusted/%s/modellist_prefixedfactorlevels.rds'))

# plot coefficients and smooths using mgcViz; first create plot lists and save; then plot to pdf

labelfunc<-function(x) {
  #relabels facet grid labels (removing alphabetic prefix)
  x<-substr(x,2,nchar(x))
  return(x)
}

outputnames_all<-c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate','Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
ggtitles_all<-c('      log10 Plasmid size (kb)\n      baseline: 10 kb\n','      Insertion sequence density\n      baseline: 0\n','      Number of other resistance gene classes\n      baseline: 0\n','      Collection date\n      baseline: initial year (1994)\n','Integron presence\nbaseline: absence','Biocide/metal resistance gene presence\nbaseline: absence','Conjugative system\nbaseline: non-mobilisable','Replicon carriage\nbaseline: untyped','Host taxonomy\nbaseline: Enterobacteriaceae','Virulence gene presence\nbaseline: absence','Geographic location\nbaseline: high-income','Isolation source\nbaseline: human')
xlabs_all<-c('log10 Plasmid size (centred on 10 kb)','Insertion sequence density (frequency per 10 kb)','Other resistance gene classes','Years since initial collection year','Integron presence','Biocide/metal resistance gene presence','Conjugative system','Replicon carriage','Host taxonomy','Virulence gene presence','Geographic location','Isolation source')
numpanels_all<-c(1,1,1,1,1,1,2,2,3,1,5,2)
paraplot_width_multiplicationfactor_all<-c(NA,NA,NA,NA,1,1,1.2,1.2,1.4,1,1.8,1.2)
#facetedparaplot_width_multiplicationfactor_all<-c(NA,NA,NA,NA,1,1,2,2,3,1,5,2) #use for 4 x 5
numcols=3
if (numcols==5) {
  facetedparaplot_width_multiplicationfactor_all<-c(NA,NA,NA,NA,1.159, 1.159, 2.120, 2.120, 3.087, 1.159, 5.145, 2.120)
  facetedparaplot_height_multiplicationfactor_all<-c(NA,NA,NA,NA,1,1,1,1,1,1,1,1)
} else {
  facetedparaplot_width_multiplicationfactor_all<-c(NA,NA,NA,NA,1.159, 1.159, 2.120, 2.120, 3.087, 1.159, 3.087, 2.120)
  facetedparaplot_height_multiplicationfactor_all<-c(NA,NA,NA,NA,1,1,1,1,1,1,1.61,1)
}


outputnames<-vector()
ggtitles<-vector()
xlabs<-vector()
numpanels<-vector()
smoothplotnames<-vector()
paraplotnames<-vector()
paraplot_width_multiplicationfactor<-vector()
facetedparaplot_width_multiplicationfactor<-vector()
facetedparaplot_height_multiplicationfactor<-vector()
numsmooths<-0
for (i in 1:length(outputnames_all)) {
  outputname<-outputnames_all[i]
  if (outputname %in% outputnames_subset) {
    outputnames<-c(outputnames,outputnames_all[i])
    ggtitles<-c(ggtitles,ggtitles_all[i])
    xlabs<-c(xlabs,xlabs_all[i])
    numpanels<-c(numpanels,numpanels_all[i])
    if (outputname %in% c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate')) {
      numsmooths<-numsmooths+1
      smoothplotnames<-c(smoothplotnames,outputnames_all[i])
    } else {
      paraplotnames<-c(paraplotnames,outputnames_all[i])
      paraplot_width_multiplicationfactor<-c(paraplot_width_multiplicationfactor,paraplot_width_multiplicationfactor_all[i])
      facetedparaplot_width_multiplicationfactor<-c(facetedparaplot_width_multiplicationfactor,facetedparaplot_width_multiplicationfactor_all[i])
      facetedparaplot_height_multiplicationfactor<-c(facetedparaplot_height_multiplicationfactor,facetedparaplot_height_multiplicationfactor_all[i])
    }
  }
}


smoothplotlist1<-list()  # logodds
smoothplotlist2<-list()  # probability
smoothplotlist3<-list()  # probability; y-axis lims=0-1
facetedparaplotlist<-list()  # logodds
# paraplotlist<-list() #probability  #!need to retrieve coefficients from facetedparaplotlist, intercept centre and plogis transform, then re-facet on a per model basis (can't really facet by factor because different models have different probability intercepts, and plotmgcviz doesn't allow this anyway) 

for (i in 1:length(outputnames)) {
  outputname<-outputnames[i]
  myggtitle<-ggtitles[i]
  myxlab<-xlabs[i]
  lowerlim<--5
  upperlim<-5
  if (outputname=='NumOtherResistanceClasses') {
    lowerlim<--3
    upperlim<-7
  }
  if (outputname=='CollectionDate') {
    lowerlim<--8
    upperlim<-2
  }
  if (outputname=='HostTaxonomy') {
    lowerlim<--6
    upperlim<-4
  }
  #create colour vector for faceted plots
  if (numpanels[i]>1) {
    mypal<-vector()
    for (palcol in brewerpal) {
      mypal<-c(mypal,rep(palcol,numpanels[i]))
    }
  } else {
    mypal<-brewerpal
  }
  if (outputname=='log10PlasmidSize' || outputname=='InsertionSequenceDensity' || outputname=='NumOtherResistanceClasses' || outputname=='CollectionDate') {
    #for smooth, create list of smooth plots to be combined in a grid
    smoothplotlistnest1<-list()
    smoothplotlistnest2<-list()
    smoothplotlistnest3<-list()
    for (j in 1:length(outcomeclasses)) {
      o<-getViz(modellist[[j]]) #this produces gamViz object for a given class, from which smooth can be extracted and plotted using sm(o,n)
      gamcoef<-coef(modellist[[j]])[1] #intercept
      o1 <- plot.mgcv.smooth.1D( sm(o, i) ) #for resistance class j, extract smooth for smooth outputname i 
      o1<-o1 + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + l_ciPoly() +l_fitLine(colour = brewerpal[j])
      o1<-o1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) +
        scale_x_continuous(breaks = scales::pretty_breaks())
      o2 <- plot.mgcv.smooth.1D( sm(o, i), trans = function(x) plogis(x+gamcoef))
      o2<-o2 + geom_hline(yintercept = plogis(gamcoef),linetype='solid',colour='light grey',size=0.3) + l_ciPoly() +l_fitLine(colour = brewerpal[j]) +
        scale_x_continuous(breaks = scales::pretty_breaks())
      #remove redundant plot elements
      if (j<=5) {
        o1<-o1 + theme(axis.text.x = element_text(colour='white'),axis.ticks.x = element_blank())
        o2<-o2 + theme(axis.text.x = element_text(colour='white'),axis.ticks.x = element_blank())
      }
      if (!j %in% c(1,6)) {
        o1<-o1+theme(axis.text.y = element_text(colour='white'),axis.ticks.y = element_blank())
        o2<-o2+theme(axis.text.y = element_text(colour='white'),axis.ticks.y = element_blank())
      }
      #extract edf
      edf<-str_extract(o1$ggObj$labels$y,"\\d+\\.*\\d*")
      o1<-o1 + ggtitle(outcomeclasses[j]) + 
        annotate('text',x=-Inf,y=Inf,label=paste('edf =',edf,collapse=''),hjust=-.1,vjust=1.8,colour="#525252",size=4)
      o2<-o2 + ggtitle(outcomeclasses[j]) + 
        annotate('text',x=-Inf,y=Inf,label=paste('edf =',edf,collapse=''),hjust=-.1,vjust=1.8,colour="#525252",size=4)
      #further customise plots
      o1<-o1 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,2,-12,-12),'pt'),axis.text=element_text(size=rel(1.1)))
      o2<-o2 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,2,-12,-12),'pt'),axis.text=element_text(size=rel(1.1)))
      #adjust probability y axis scale
      o3<-o2
      o3<-o3 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,1))
      if (outputname=='NumOtherResistanceClasses') {
        o2<-o2 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,1))
      } else {
        o2<-o2 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,0.2))
      }
      #assign plots to list
      smoothplotlistnest1[[outcomeclasses[j]]]<-o1 #logodds
      smoothplotlistnest2[[outcomeclasses[j]]]<-o2 #prob
      smoothplotlistnest3[[outcomeclasses[j]]]<-o3
    }
    smoothplotlist1[[outputname]]<-smoothplotlistnest1
    smoothplotlist2[[outputname]]<-smoothplotlistnest2
    smoothplotlist3[[outputname]]<-smoothplotlistnest3
  } else {
    #for parametric terms, create log odds scale plots, faceted by factor level, across resistance classes (rather than having 1 plot per resistance class and combining in grid)
    p1<-plot.mgamViz(modellist, select = i,allTerms = T,a.facet = list(labeller=as_labeller(labelfunc),as.table=FALSE,ncol=numcols)) #for paramtric terms, create single plot for all coefficients
    p1<-p1 + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + l_ciBar(col='grey 42',linetype=1,width=0.5,size=0.4) + l_fitPoints(col=mypal,size=2) + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)),axis.title.y = element_text(size=rel(1.3)),axis.title.x = element_blank(),plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
    p1<-p1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim))
    facetedparaplotlist[[outputname]]<-p1
  }
}


# plotting; save plot lists to file; save plots to pdf

# save lists to file
saveRDS(smoothplotlist1,file = gsub('%s',modelname,'output_adjusted/%s/coefficientplots/logoddsscale/smoothplotlist1.rds'))
saveRDS(smoothplotlist2,file = gsub('%s',modelname,'output_adjusted/%s/coefficientplots/probscale/smoothplotlist2.rds'))
saveRDS(smoothplotlist3,file = gsub('%s',modelname,'output_adjusted/%s/coefficientplots/probscale/smoothplotlist3.rds'))
saveRDS(facetedparaplotlist,file = gsub('%s',modelname,'output_adjusted/%s/coefficientplots/logoddsscale/facetedparaplotlist.rds'))
#saveRDS(paraplotlist,file = gsub('%s',modelname,'output_adjusted/%s/coefficientplots/probscale/paraplotlist.rds'))

# smooth plots
#width=20
#height=7.4
#height=4.8
width=12
height=5.6
probwidth=12.7  # need to adjust for longer ylab text
for (i in 1:length(smoothplotnames)) {
  smoothplotname<-smoothplotnames[i]
  print(smoothplotname)
  pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=smoothplotname),'output_adjusted/%1/coefficientplots/logoddsscale/%2.pdf'),width=width,height=height)
  gridPrint(grobs=(smoothplotlist1[[smoothplotname]][order(outcomeclasses)]),nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15)),left=textGrob('Effect on log odds',gp=gpar(fontsize=15),rot=90),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
  dev.off()
  pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=smoothplotname),'output_adjusted/%1/coefficientplots/probscale/%2.pdf'),width=probwidth,height=height)
  gridPrint(grobs=(smoothplotlist2[[smoothplotname]][order(outcomeclasses)]),nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15)),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
  dev.off()
  pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=smoothplotname),'output_adjusted/%1/coefficientplots/probscale/ylim0to1_%2.pdf'),width=probwidth,height=height)
  gridPrint(grobs=(smoothplotlist3[[smoothplotname]][order(outcomeclasses)]),nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15)),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
  dev.off()
}

# parametric plots
#faceted parametric coefficient plots (logodds scale)
#width=4
#height=5
width=3.2
height=4.5
for (i in 1:length(paraplotnames)) {
  paraplotname<-paraplotnames[i]
  print(paraplotname)
  mywidth<-width*facetedparaplot_width_multiplicationfactor[i]
  myheight<-height*facetedparaplot_height_multiplicationfactor[i]
  pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=paraplotname),'output_adjusted/%1/coefficientplots/logoddsscale/%2_faceted.pdf'),width=mywidth,height=myheight)
  print(facetedparaplotlist[[paraplotname]])
  dev.off()
}


# ---------------------------
# OLD CODE

# legend for smooth plots (not needed)
#dummydf<-data.frame(outcomeclasses,brewerpal,1:length(outcomeclasses))
#colnames(dummydf)<-c('classes','cols','data')
#p<-ggplot(dummydf,aes(x=classes,y=data,colour=classes)) + geom_point() + scale_colour_manual(values=as.vector(dummydf$cols)[order(outcomeclasses)]) +labs(colour='Resistance class')  # scale_colour_manual(values=as.vector(dummydf$cols),limits=outcomeclasses)  # this places in outcomeclasses order

#smoothlegend<-get_legend(p)
#pdf(gsub('%s',modelname,'output_adjusted/%s/coefficientplots/smoothplots_legend.pdf'))
#plot(smoothlegend)
#dev.off()
