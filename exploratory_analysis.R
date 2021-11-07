pacman::p_load(tidyverse, gsubfn, knitr, rio, countrycode, cowplot, car, rcompanion, rstatix, eply)
source('functions.R')
dir.create('output_exploratory', showWarnings = FALSE)

# load data
plasmiddf<-convert(in_file = 'data/Appendix_2.xlsx',out_file='data/Appendix_2G.tsv',in_opts=list(sheet=7))  # data is quoted to preserve dates; convert to tsv and re-load data
plasmiddf<-read.table('data/Appendix_2G.tsv',header=TRUE,quote="'",sep='\t',as.is=TRUE)
plasmiddf<-plasmiddf[plasmiddf$InFinalDataset==TRUE,]
#nrow(plasmiddf) #14143
plasmiddf <- plasmiddf %>% rename(carbapenem = betalactam_carbapenem, ESBL = betalactam_ESBL, `TEM-1` = betalactam_TEM.1)
outcomeclasses<-c('aminoglycoside','phenicol','sulphonamide','tetracycline','macrolide','TEM-1','trimethoprim','ESBL', 'carbapenem','quinolone','colistin')

# ---------------------------
# resistance class outcome variables

# get colsums for resistance genes; convert resistance outcomes to binary; write totals/binary totals to file
resgenedf<-plasmiddf[,c('TotalResGenes',outcomeclasses)]
resgenestotal<-as.data.frame(rev(sort(colSums(resgenedf))))
# total resistance genes by class
colnames(resgenestotal)<-'Total resistance genes by class'
write.table(resgenestotal,file='output_exploratory/resgenestotalbyclass.tsv',sep='\t',col.names = TRUE,row.names = TRUE)
# binary totals (i.e. per-plasmid presence/absence) of resistance genes by class
resgenedfbinary<-resgenedf
resgenedfbinary[resgenedfbinary>0]<-1
resgenestotalplasmids<-as.data.frame(rev(sort(colSums(resgenedfbinary))))
colnames(resgenestotalplasmids)<-'Total plasmids with a resistance gene by class'
write.table(resgenestotalplasmids,file='output_exploratory/resgenestotalplasmidsbyclass.tsv',sep='\t',col.names = TRUE,row.names = TRUE)


# plot resistance class interaction heatmap (coloured by proportion of total class)
orderoutcomeclasses<-match(rownames(resgenestotalplasmids),outcomeclasses)
orderoutcomeclasses<-orderoutcomeclasses[!is.na(orderoutcomeclasses)]
#resgeneclasseshm<-resclassheatmap(resgenedfbinary,outcomeclasses[orderoutcomeclasses],twoway=TRUE)
resgeneclasseshm<-resclassheatmap(resgenedfbinary,outcomeclasses,twoway=TRUE)
pdf('output_exploratory/resgeneclasseshm.pdf')
resgeneclasseshm
dev.off()


# ---------------------------
# explanatory variables

PlasmidSize<-as.numeric(plasmiddf$SequenceLength)  # plasmid size (kb)



CollectionDate<-as.numeric(unlist(lapply(strsplit(as.character(plasmiddf$CollectionDate),'-',fixed=TRUE), function(x) x[1])))  # collection date
#sum(is.na(CollectionDate))  # 7768 - ~half of collection dates are missing - impute with createdates
CreateDate<-as.numeric(unlist(lapply(strsplit(as.character(plasmiddf$CreateDate),'-'), function(x) x[1])))  # all accessions have createdate
CollectionDateImputed<-CollectionDate
CollectionDateImputed[is.na(CollectionDateImputed)]<-CreateDate[is.na(CollectionDateImputed)]

# scatter plot of collection date vs createdate
colcreatedatedf<-data.frame(CollectionDate[!is.na(CollectionDate)],CreateDate[!is.na(CollectionDate)])
colnames(colcreatedatedf)<-c('CollectionDate','CreateDate')
colcreatepearson=round(cor(colcreatedatedf$CollectionDate,colcreatedatedf$CreateDate,method='pearson'),3)
sp<-ggplot(colcreatedatedf, aes(x=CollectionDate, y=CreateDate))
sp<-sp + stat_bin2d(bins=50)
spdat<-ggplot_build(sp)
spdat<-spdat[[1]][[1]][,'count']
sp<-sp + scale_fill_gradient(low="lightblue", high="red", limits=c(0, max(spdat))) + 
  ggtitle(gsubfn('%1|%2',list('%1'=colcreatepearson,'%2'=nrow(colcreatedatedf)),"Pearson's correlation (r)=%1; n=%2")) + theme_bw() +
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab('Collection Date') + ylab('Create Date') + scale_x_continuous(limits=c(1910,2020),breaks=c(1930,1950,1970,1990,2010))

pdf('output_exploratory/CollectionDate_CreateDate_scatter.pdf',5,5)
print(sp)
dev.off()

CollectionDate<-CollectionDateImputed  # re-assign imputed collection date as CollectionDate


BiocideMetalResistance<-plasmiddf$NumBacMetGenes  # biocide/metal resistance genes detected using BacMet database
BiocideMetalResistance[BiocideMetalResistance>0]<-1  # convert to binary
BiocideMetalResistance<-as.factor(BiocideMetalResistance)


Virulence<-plasmiddf$NumVFDBgenes  # virulence genes detected using VFDB database
Virulence[Virulence>0]<-1  # convert to binary
Virulence<-as.factor(Virulence)


Integron<-plasmiddf$NumIntegron  # complete integrons
Integron[Integron>0]<-1  # convert to binary
Integron<-as.factor(Integron)


NumInsertionSequences<-plasmiddf$NumIS
InsertionSequenceDensity<-plasmiddf$NumIS/(plasmiddf$SequenceLength/10000)


conj<-plasmiddf$NumMobType
conj[conj>0]<-"mobilisable"
conj[conj=="mobilisable" & plasmiddf$NumConjType>0]<-"conjugative"
conj[conj==0]<-"non-mobilisable"
ConjugativeSystem<-factor(conj, ordered = FALSE,levels = c("non-mobilisable", "mobilisable", "conjugative"))


# country coding: include frequent countries (US, mainland China) and EU as categories; and include "Regions as defined in the World Bank Development Indicators":
# WB_HI (labelled high-income) and WB_LMI/WB_UMI (labelled middle-income); all others categorised as 'other'
countries<-plasmiddf$Country
countries[(grepl(glob2rx("R*union") , countries))]<-'Reunion'  # handle Reunion corrupted accent

regions<-countrycode(sourcevar = countries,origin = 'country.name', destination='region')  # Some values were not matched unambiguously...coded as NA
for (i in 1:length(regions)) {
  if (is.na(regions[i])) {
    regions[i]<-'-'
  }
}

countryregiondf<-as.data.frame(cbind(countries,regions))

eucountries<-countrycode(sourcevar = countries,origin = 'country.name', destination='eu28')
for (i in 1:length(eucountries)) {
  if (is.na(eucountries[i])) {
    eucountries[i]<-'-'
  }
}

# GlobalHealthObservatoryMetadata downloaded from https://apps.who.int/gho/data/node.metadata.COUNTRY?lang=en; minor edits to conform to country names in plasmiddf e.g. Bahamas -> The Bahamas; Gambia -> The Gambia; 
# Additinal discrepancies handled in code below:
# Raas Cabaad, Somalia ; change to Somalia
# Antarctica; change to '-'
# Svalbard and Jan Mayen; change to Norway
# United States Minor Outlying Islands; change to United States
# Reunion: remove accent

incomedf<-read.table('data/GlobalHealthObservatoryMetadata_edit.tsv',header=TRUE,sep='\t',as.is = TRUE,quote = "",comment.char = "")
incomedf$DisplayString<-sapply(incomedf$DisplayString, function(x) unquote(x))

missingincomes<-vector()
incomes<-vector()
for (country in countries) {
  if (country!='-' && country!='Antarctica') {
    if (country=='Raas Cabaad, Somalia') {
      country='Somalia'
    }
    if (country=='Svalbard and Jan Mayen') {
      country='Norway'
    }
    if (country=='United States Minor Outlying Islands') {
      country='United States'
    }
    
    incomeindx<-which(incomedf$DisplayString==country)
    stopifnot(length(incomeindx)==1)
    income<-incomedf$WORLD_BANK_INCOME_GROUP_CODE[incomeindx]
    if (nchar(income)==0) {
      missingincomes<-c(missingincomes,country)
      income<-'-'
    }
    incomes<-c(incomes,income)
  } else {
    incomes<-c(incomes,'-')
  }
}

countryregiondf<-as.data.frame(cbind(countries,eucountries,incomes))
#rev(sort(table(as.character(countryregiondf[countryregiondf$incomes=='-',1])))) #countries associated with no income category...
#-     Taiwan Antarctica    Reunion    Mayotte Martinique 
#7284         93         58         10          3          2

#create final geography vector
WB_UMI_count<-0
WB_LMI_count<-0
geographies<-character(nrow(countryregiondf))
for (i in 1:nrow(countryregiondf)) {
  #print(i)
  country<-as.character(countryregiondf[i,"countries"])
  eucountry<-as.character(countryregiondf[i,"eucountries"])
  income<-as.character(countryregiondf[i,"incomes"])
  if (country=='-') {
    geographies[i]<-'other'
  } 
  else if (country=='United States' || country=='China') {
    geographies[i]<-country
  }
  else if (eucountry=='EU') {
    geographies[i]<-'EU'
  }
  else {
    if (income=='WB_LI') {
      income='other'
    }
    else if (income=='WB_UMI' || income=='WB_LMI') {
      if (income=='WB_UMI') {
        WB_UMI_count<-WB_UMI_count+1
      } else {
        WB_LMI_count<-WB_LMI_count+1
      }
      income='middle-income'
    }
    else if (income=='WB_HI') {
      income='high-income'
    }
    else {
      income='other'
    }
    geographies[i]<-income
  }
}
#> WB_UMI_count
#[1] 593
#> WB_LMI_count
#[1] 362

GeographicLocation<-factor(geographies, ordered = FALSE,levels = c("high-income", "middle-income", "China", "United States", "EU", "other"))


# isolation source; refactored original variable (IsolationSource_unmergedfactorlevels), merging livestock categories and assigning non-livestock agriculture as "other", to avoid low numbers
isolationsources<-vector()
livestockpattern<-paste(c('aquaculture','cow','chicken','pig','turkey','sheep','poultry','goat','goose','duck','cow;pig'),collapse='|')
for (i in 1:nrow(plasmiddf)) {
  source<-plasmiddf$IsolationSource_unmergedfactorlevels[i]
  if (source=='human' || source=='sewage') {
    isolationsources<-c(isolationsources,'human')
  }
  else if (source=='agriculture') {
    isolationsources<-c(isolationsources,'other')
  }
  else if (length(grep(livestockpattern,source,value=T))==1) {
    isolationsources<-c(isolationsources,'livestock')
  }
  else {
    isolationsources<-c(isolationsources,'other')
  }
}
IsolationSource<-factor(isolationsources, ordered = FALSE,levels = c("human", "livestock","other"))


# replicon carriage: untyped/single/multi-replicon type carriage
reptypes<-vector()
for (i in 1:length(plasmiddf$RepliconType)) {
  reptype<-plasmiddf$RepliconType[i]
  #check for untyped plasmids
  if (reptype=='-') {
    reptypes<-c(reptypes,'untyped')
    next
  }
  #assign single/multi rep type - based on whether there are multiple unique replicon types (not family level)
  reptypeunique<-unique(unlist(strsplit(reptype,',')))
  if (length(reptypeunique)==1) {
    reptypes<-c(reptypes,'single-replicon')
  }
  else {
    reptypes<-c(reptypes,'multi-replicon')
  }
}
RepliconCarriage<-factor(reptypes, ordered = FALSE,levels = c("untyped", "single-replicon", "multi-replicon"))


# host taxonomy
taxa<-vector()
for (i in 1:nrow(plasmiddf)) {
  family<-plasmiddf$Family[i]
  phylum<-plasmiddf$Phylum[i]
  if (family=='Enterobacteriaceae') {
    taxa<-c(taxa,family)
  }
  else if (phylum=='Proteobacteria' || phylum=='Firmicutes') {
    taxa<-c(taxa,phylum)
  }
  else {
    taxa<-c(taxa,'other')
  }
}

taxa[taxa=='Proteobacteria']<-'Proteobacteria_other'  # (non-Enterobacteriaceae Proteobacteria)
HostTaxonomy<-factor(taxa, ordered = FALSE,levels = c("Enterobacteriaceae", "Proteobacteria_other", "Firmicutes","other"))


# ---------------------------
# save explanatory and outcome variables to file; create per resistance class NumOtherResistanceClasses explanatory variables

outcomeclasseslist<-vector('list',length(outcomeclasses))
NumOtherResistanceClasseslist<-vector('list',length(outcomeclasses))

for (i in 1:length(outcomeclasses)) {
  #outcome
  outcomeclass<-outcomeclasses[i]
  outcomeclassbinary<-as.factor(resgenedfbinary[,outcomeclass])
  outcomeclasseslist[[i]]<-outcomeclassbinary
  #number of other resistance gene classes (for a given outcome class)
  otherresgenedf<-resgenedf[,-which(colnames(resgenedf) %in% c('TotalResGenes',outcomeclass))]
  otherresgenedf[otherresgenedf>1]<-1
  NumOtherResistanceClasses<-rowSums(otherresgenedf)
  NumOtherResistanceClasseslist[[i]]<-NumOtherResistanceClasses
}

outcomeclasseslistdf<-do.call(cbind.data.frame, outcomeclasseslist)
NumOtherResistanceClasseslistdf<-do.call(cbind.data.frame, NumOtherResistanceClasseslist)
colnames(outcomeclasseslistdf)<-paste('outcome',outcomeclasses,sep='')
colnames(NumOtherResistanceClasseslistdf)<-paste('NumOtherResistanceClasses',outcomeclasses,sep='')

# create final dataframe of explanatory and outcome variables
Accession<-plasmiddf$Accession
log10PlasmidSize<-log10(PlasmidSize)  # this will be over-written after data transformation
finaldf<-data.frame(Accession,PlasmidSize,log10PlasmidSize,CollectionDate,NumInsertionSequences,InsertionSequenceDensity,BiocideMetalResistance,Virulence,Integron,ConjugativeSystem,GeographicLocation,IsolationSource,RepliconCarriage,HostTaxonomy)
finaldf<-cbind(finaldf,outcomeclasseslistdf,NumOtherResistanceClasseslistdf)

# check if any values are NA
for (i in 1:ncol(finaldf)) {
  if (any(is.na(finaldf[,i]))==TRUE)
    print(colnames(finaldf)[i])
}


# cross-tabulate categorical predictors and outcome variables (to ensure factor levels balance bias-variance trade-off and for later calcualtion of unadjusted odds ratios)
factorvars<-c('BiocideMetalResistance','ConjugativeSystem','GeographicLocation','Integron','IsolationSource','RepliconCarriage','HostTaxonomy','Virulence')
counter<-0
crosstablelist<-list()
for (outcomeclass in outcomeclasses) {
  for (factorvar in factorvars) {
    counter<-counter+1
    crosstablelist[[counter]]<-crosstabulate(outcomeclass,factorvar)
  }
}

crosstabledf<-as.data.frame(do.call(rbind,crosstablelist))
write.table(crosstabledf,'output_exploratory/crosstabulation.tsv',col.names = TRUE,row.names = FALSE,sep='\t',quote=FALSE)


# ---------------------------
# pre-modelling transformations; truncate back (=winsorise) outliers for continuous variables: PlasmidSize (truncate both tails' log10; centre on 100 kb), NumInsertionSequences (truncate right tail), CollectionDate (truncate left tail); NumOtherResistanceClasses (truncate right tail)
finaldftrunc<-finaldf
PlasmidSize95pct<-round(quantile(finaldftrunc$PlasmidSize, probs = c(0.05, 0.95))[2])
finaldftrunc$PlasmidSize[finaldftrunc$PlasmidSize>PlasmidSize95pct]<-PlasmidSize95pct
PlasmidSize5pct<-round(quantile(finaldftrunc$PlasmidSize, probs = c(0.05, 0.95))[1])
finaldftrunc$PlasmidSize[finaldftrunc$PlasmidSize<PlasmidSize5pct]<-PlasmidSize5pct
finaldftrunc$log10PlasmidSize<-log10(finaldftrunc$PlasmidSize)  # log10PlasmidSize derived from truncated PlasmidSize
NumInsertionSequences95pct<-round(quantile(finaldftrunc$NumInsertionSequences, probs = c(0.05, 0.95))[2])
finaldftrunc$NumInsertionSequences[finaldftrunc$NumInsertionSequences>NumInsertionSequences95pct]<-NumInsertionSequences95pct

InsertionSequenceDensity95pct<-round(quantile(finaldftrunc$InsertionSequenceDensity, probs = c(0.05, 0.95)),2)[2]
finaldftrunc$InsertionSequenceDensity[finaldftrunc$InsertionSequenceDensity>InsertionSequenceDensity95pct]<-InsertionSequenceDensity95pct

CollectionDate5pct<-round(quantile(finaldftrunc$CollectionDate, probs = c(0.05, 0.95))[1])
finaldftrunc$CollectionDate[finaldftrunc$CollectionDate<CollectionDate5pct]<-CollectionDate5pct

NumOtherResistanceClassestrunccutoffs<-list()
for (outcomeclass in outcomeclasses) {
  NumOtherResistanceClasses95pct<-round(quantile(finaldftrunc[,gsub('%s',outcomeclass,'NumOtherResistanceClasses%s')], probs = c(0.05, 0.95))[2])
  finaldftrunc[,gsub('%s',outcomeclass,'NumOtherResistanceClasses%s')][finaldftrunc[,gsub('%s',outcomeclass,'NumOtherResistanceClasses%s')]>NumOtherResistanceClasses95pct]<-NumOtherResistanceClasses95pct
  NumOtherResistanceClassestrunccutoffs[[outcomeclass]]<-NumOtherResistanceClasses95pct
}


# express collection date as date since 5% percentile (1994)
finaldftrunc$CollectionDate<-finaldftrunc$CollectionDate-round(as.numeric(CollectionDate5pct))

# express log10PlasmidSize as log10PlasmidSize - 4 (i.e. centred on 10kb since 0kb is not meaningful)
finaldftrunc$log10PlasmidSize<-finaldftrunc$log10PlasmidSize - 4

# write to file
write.table(finaldftrunc,'data/plasmiddf_transformed.tsv',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)  # this file is appendix 2H

# write truncation thresholds to file
trunccutofffilename<-'data/truncationcutoffs.tsv'
cat('',file=trunccutofffilename,append=FALSE)
cat(gsub('%s',PlasmidSize95pct,'PlasmidSize95pct: %s\n'),file=trunccutofffilename,append=TRUE)
cat(gsub('%s',PlasmidSize5pct,'PlasmidSize5pct: %s\n'),file=trunccutofffilename,append=TRUE)
cat(gsub('%s',NumInsertionSequences95pct,'NumInsertionSequences95pct: %s\n'), file=trunccutofffilename,append=TRUE)  
cat(gsub('%s',InsertionSequenceDensity95pct,'InsertionSequenceDensity95pct: %s\n'), file=trunccutofffilename,append=TRUE) 
cat(gsub('%s',CollectionDate5pct,'CollectionDate5pct: %s\n'),file=trunccutofffilename,append=TRUE)
for (outcomeclass in outcomeclasses) {
  NumOtherResistanceClasses95pct<-NumOtherResistanceClassestrunccutoffs[[outcomeclass]]
  cat(gsubfn('%1|%2',list('%1'=outcomeclass,'%2'=NumOtherResistanceClasses95pct),'%1 NumOtherResistanceClasses95pct: %2\n'),file=trunccutofffilename,append=TRUE)
}


# ---------------------------
# With truncated dataset, flag potential issues with collinearity using association statistics
print('exploring collinearity...')
finaldftrunccopy<-finaldftrunc
for (binaryvar in c('Integron','BiocideMetalResistance','Virulence')) {
  finaldftrunccopy[,binaryvar]<-as.numeric(as.character(finaldftrunccopy[,binaryvar]))
}

# Matrix outputs: spearman's correlation; Kramer's V for association between categorical variables (0=no association; 1=perfect association); Kruskal-Wallis H test eta-squared statistic for association between continuoius and categorical

# Numeric: log10PlasmidSize, NumInsertionSequences, InsertionSequenceDensity, NumOtherResistanceClasses, CollectionDate
# Categorical:
#    Binary: Integron, BiocideMetalResistance, Virulence
#    Nominal: ConjugativeSystem, RepliconCarriage, HostTaxonomy, GeographicLocation, IsolationSource

# Numeric-numeric (Spearman's correlation)
# Generate vector of NumOtherResistanceClasses variable across outcomes
NumOtherResistanceClassesvec<-vector()
for (outcomeclass in outcomeclasses) {
  NumOtherResistanceClasses<-gsub('%s',outcomeclass,'NumOtherResistanceClasses%s')
  NumOtherResistanceClassesvec<-c(NumOtherResistanceClassesvec,NumOtherResistanceClasses)
}

numericvars<-c('log10PlasmidSize','NumInsertionSequences','InsertionSequenceDensity',NumOtherResistanceClassesvec,'CollectionDate')
spearmans_matrix<-matrix(NA,nrow=length(numericvars),ncol=length(numericvars))
colnames(spearmans_matrix)<-numericvars
rownames(spearmans_matrix)<-numericvars

for (i in 1:length(numericvars)) {
  numericvar1<-numericvars[i]
  for (j in 1:length(numericvars)) {
    if (i==j || j<i) {
      next
    }
    numericvar2<-numericvars[j]
    out<-as.numeric(cor(finaldftrunc[,numericvar1],finaldftrunc[,numericvar2],method = 'spearman'))
    spearmans_matrix[i,j]<-out
  }
}
spearmans_df<-as.data.frame(spearmans_matrix)
spearmans_df_abs<-abs(spearmans_df)
write.table(spearmans_df,file='output_exploratory/cor_spearmans_numeric_numeric_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)
write.table(spearmans_df_abs,file='output_exploratory/cor_spearmans_numeric_numeric_abs_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)

# Categorical-categorical (Cramers V)
categoricalvars<-c('Integron','BiocideMetalResistance','Virulence','ConjugativeSystem','RepliconCarriage','HostTaxonomy','GeographicLocation','IsolationSource')

cramerV_matrix<-matrix(NA,nrow=length(categoricalvars),ncol=length(categoricalvars))
colnames(cramerV_matrix)<-categoricalvars
rownames(cramerV_matrix)<-categoricalvars

for (i in 1:length(categoricalvars)) {
  categoricalvar1<-categoricalvars[i]
  for (j in 1:length(categoricalvars)) {
    if (i==j || j<i) {
      next
    }
    categoricalvar2<-categoricalvars[j]
    crosstable<-table(finaldftrunc[,categoricalvar1],finaldftrunc[,categoricalvar2])
    out<-as.numeric(cramerV(crosstable,bias.correct = TRUE))
    cramerV_matrix[i,j]<-out
  }
}
cramerV_df<-as.data.frame(cramerV_matrix)
write.table(cramerV_df,file='output_exploratory/cor_cramerV_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)

# Numeric-categorical (Kruskal-Wallis H)
kruskalH_matrix<-matrix(NA,nrow=length(numericvars),ncol=length(categoricalvars))
colnames(kruskalH_matrix)<-categoricalvars
rownames(kruskalH_matrix)<-numericvars
kruskalpvalue_matrix<-matrix(NA,nrow=length(numericvars),ncol=length(categoricalvars))
colnames(kruskalpvalue_matrix)<-categoricalvars
rownames(kruskalpvalue_matrix)<-numericvars
kruskaletasqrd_matrix<-matrix(NA,nrow=length(numericvars),ncol=length(categoricalvars))
colnames(kruskaletasqrd_matrix)<-categoricalvars
rownames(kruskaletasqrd_matrix)<-numericvars

for (i in 1:length(numericvars)) {
  numericvar<-numericvars[i]
  for (j in 1:length(categoricalvars)) {
    #print(c(i,j))
    categoricalvar<-categoricalvars[j]
    out<-kruskal.test(finaldftrunc[,numericvar],finaldftrunc[,categoricalvar])
    Hstatistic<-as.numeric(out$statistic)
    pvalue<-as.numeric(out$p.value)
    kruskalH_matrix[i,j]<-Hstatistic
    kruskalpvalue_matrix[i,j]<-pvalue
    kruskaletasqrd_matrix[i,j]<-as.numeric(epsilonSquared(finaldftrunc[,numericvar],finaldftrunc[,categoricalvar]))
  }
}


kruskalH_matrix<-as.data.frame(kruskalH_matrix)
write.table(kruskalH_matrix,file='output_exploratory/cor_kruskalH_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)
kruskalpvalue_matrix<-as.data.frame(kruskalpvalue_matrix)
write.table(kruskalpvalue_matrix,file='output_exploratory/cor_kruskalpvalue_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)
kruskaletasqrd_matrix<-as.data.frame(kruskaletasqrd_matrix)
write.table(kruskaletasqrd_matrix,file='output_exploratory/cor_kruskaletasqrd_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)


# Numeric-binary
#https://www.statology.org/point-biserial-correlation-r/
#Can't use point-biserial because normality assumption is violated; use Spearman's instead

binaryvars<-c('Integron','BiocideMetalResistance','Virulence')

biserialcor_matrix<-matrix(NA,nrow=length(numericvars),ncol=length(binaryvars))
colnames(biserialcor_matrix)<-binaryvars
rownames(biserialcor_matrix)<-numericvars

for (i in 1:length(numericvars)) {
  numericvar<-numericvars[i]
  for (j in 1:length(binaryvars)) {
    binaryvar<-binaryvars[j]
    out<-cor.test(as.numeric(as.character(finaldftrunc[,numericvar])),as.numeric(finaldftrunc[,binaryvar]),method = 'spearman')
    out<-as.numeric(out$estimate)
    biserialcor_matrix[i,j]<-out
  }
}

biserialcor_matrix<-as.data.frame(biserialcor_matrix)
biserialcor_matrix_abs<-abs(biserialcor_matrix)
write.table(biserialcor_matrix,file='output_exploratory/cor_spearmans_numeric_binary_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)
write.table(biserialcor_matrix_abs,file='output_exploratory/cor_spearmans_numeric_binary_abs_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)

