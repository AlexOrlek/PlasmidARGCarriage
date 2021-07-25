library(ggplot2)
library(gsubfn)
library(knitr)
library(rio)
library(countrycode)
library(cowplot)
library(car)
library(rcompanion)
library(rstatix)
library(eply)
setwd("C:/Users/alexo/Dropbox/Oxford_2015/thesis/CORRECTIONS/Chapter_5/resistance_prediction")
source('exploratory_analysis_functions.R')


#outcomeclasses<-c('aminoglycoside','betalactam','macrolide','phenicol','quinolone','sulphonamide','tetracycline','trimethoprim')
outcomeclasses<-c('aminoglycoside','betalactam_carbapenem','betalactam_ESBL','betalactam_other','macrolide','phenicol','quinolone','sulphonamide','tetracycline','trimethoprim')

plotexploratory=TRUE  #should exploratory plots be generated?
plotpairs=TRUE  #should pairs plots be produced (slow)

###load data
#?convert
#plasmiddf<-convert(in_file = 'data/Appendix_D_Table_6_splitbetalactam.xlsx',out_file='data/Appendix_D_Table_6_splitbetalactam.tsv',in_opts=list(sheet=1))
#rm(plasmiddf)
plasmiddf<-read.table('data/Appendix_D_Table_6_splitbetalactam.tsv',header=TRUE,quote="'",sep='\t',as.is=TRUE)
plasmiddf<-plasmiddf[plasmiddf$InDeduplicatedDataset==TRUE,]
#nrow(plasmiddf) #14143


###plot resistance gene class stats (class counts, class intersections)

##get colsums for resistance genes; convert to binary; write totals/binary totals to file
resgenedf<-plasmiddf[,c('TotalResGenes',outcomeclasses)]
resgenestotal<-as.data.frame(rev(sort(colSums(resgenedf))))
#total resistance genes by class
colnames(resgenestotal)<-'Total resistance genes by class'
kable(resgenestotal)
write.table(resgenestotal,file='output_exploratory/resgenestotalbyclass.tsv',sep='\t',col.names = TRUE,row.names = TRUE)
#binary totals (i.e. per-plasmid presence/absence) of resistance genes by class
resgenedfbinary<-resgenedf
resgenedfbinary[resgenedfbinary>0]<-1
resgenestotalplasmids<-as.data.frame(rev(sort(colSums(resgenedfbinary))))
colnames(resgenestotalplasmids)<-'Total plasmids with a resistance gene by class'
kable(resgenestotalplasmids)
write.table(resgenestotalplasmids,file='output_exploratory/resgenestotalplasmidsbyclass.tsv',sep='\t',col.names = TRUE,row.names = TRUE)


#plot resistance class interaction heatmap (coloured by proportion of total class)
orderoutcomeclasses<-match(rownames(resgenestotalplasmids),outcomeclasses)
orderoutcomeclasses<-orderoutcomeclasses[!is.na(orderoutcomeclasses)]
resgeneclasseshm<-resclassheatmap(resgenedfbinary,outcomeclasses[orderoutcomeclasses],twoway=TRUE)
print(resgeneclasseshm)
pdf('output_exploratory/resgeneclasseshm.pdf')
resgeneclasseshm
dev.off()



###assign predictors; plot coldate/createdate scatter plot

#size
plasmidsize<-as.numeric(plasmiddf$SequenceLength)
colnames(plasmiddf)
#all(is.na(plasmidsize)==F) #TRUE

#collection date
coldate<-as.numeric(unlist(lapply(strsplit(plasmiddf$CollectionDate,'-',fixed=TRUE), function(x) x[1])))

#sum(is.na(coldate)) #7768 - half of collection dates are missing - impute with createdates
createdate<-as.numeric(unlist(lapply(strsplit(plasmiddf$CreateDate,'-'), function(x) x[1]))) #all accessions have createdate
coldateimputed<-coldate
coldateimputed[is.na(coldateimputed)]<-createdate[is.na(coldateimputed)]

#scatter plot of collection date vs createdate
colcreatedatedf<-data.frame(coldate[!is.na(coldate)],createdate[!is.na(coldate)])
colnames(colcreatedatedf)<-c('CollectionDate','CreateDate')
#min(colcreatedatedf$CollectionDate) #1911
#min(colcreatedatedf$CreateDate) #1998
colcreatepearson=round(cor(colcreatedatedf$CollectionDate,colcreatedatedf$CreateDate,method='pearson'),3)
colcreatenrow<-nrow(colcreatedatedf) #6375
sp<-ggplot(colcreatedatedf, aes(x=CollectionDate, y=CreateDate))
sp<-sp + stat_bin2d(bins=50)
spdat<-ggplot_build(sp)
spdat<-spdat[[1]][[1]][,'count']

sp<-sp + scale_fill_gradient(low="lightblue", high="red", limits=c(0, max(spdat))) + 
  ggtitle(gsubfn('%1|%2',list('%1'=colcreatepearson,'%2'=colcreatenrow),"Pearson's correlation (r)=%1; n=%2")) + 
  theme(plot.title = element_text(size=10),panel.grid.major = element_blank(),panel.grid.minor = element_blank()) + xlab('Collection Date') + ylab('Create Date') + scale_x_continuous(limits=c(1910,2020),breaks=c(1930,1950,1970,1990,2010))
print(sp)
pdf('output_exploratory/CollectionDate_CreateDate_scatter.pdf',5,5)
print(sp)
dev.off()

#bacmet
numbacmet<-plasmiddf$NumBacMetGenes
bacmetbinary<-plasmiddf$NumBacMetGenes
bacmetbinary[bacmetbinary>0]<-1
bacmetbinary<-as.factor(bacmetbinary)

#vfdb
#rev(sort(table(plasmiddf$NumVFDBgenes)))  #biologically, having more virulence genes is not necessarily going to increase virulence relative to having one potent virulence gene - so makes sense to model as binary (and there are too few >0 anyway)
numvfdb<-plasmiddf$NumVFDBgenes
vfdbbinary<-plasmiddf$NumVFDBgenes
vfdbbinary[vfdbbinary>0]<-1
vfdbbinary<-as.factor(vfdbbinary)

#integron
#rev(sort(table(plasmiddf$NumIntegron)))  #few above 1
numintegron<-plasmiddf$NumIntegron
integronbinary<-plasmiddf$NumIntegron
integronbinary[integronbinary>0]<-1
integronbinary<-as.factor(integronbinary)

#IS
#rev(sort(table(plasmiddf$NumIS)))
numis<-plasmiddf$NumIS
isscore<-plasmiddf$NumIS/(plasmiddf$SequenceLength/10000)
isbinary<-plasmiddf$NumIS
isbinary[isbinary>0]<-1
isbinary<-as.factor(isbinary)

#conjscan
conj<-plasmiddf$NumMobType
conj[conj>0]<-1
conj[conj==1 & plasmiddf$NumConjType>0]<-2  #conj==1 i.e. mob>0 is redundant because all conj systems have mob
conj[conj==0]<-"non-mob"
conj[conj==1]<-"mob"
conj[conj==2]<-"conj"
conj<-factor(conj, ordered = FALSE,levels = c("non-mob", "mob", "conj"))

###COUNTRY CODING - FINAL APPROACH: include common countries (US, mainland China) and EU as categories; and include "Regions as defined in the World Bank Development Indicators": WB_HI (labelled High-income) and WB_LMI/WB_UMI (labelled Middle-income); all others categorised as 'other'

countries<-plasmiddf$Country
countries[countries=='RÃ©union']<-'Réunion'
#rev(sort(table(countries)))  #common countries: US, China [only countries making up >=10%] [South Korea, Japan - east asia?]
#length(table(countries)) #125
#any(is.na(countries))


regions<-countrycode(sourcevar = countries,origin = 'country.name', destination='region') #Some values were not matched unambiguously: -, Antarctica, United States Minor Outlying Islands; coded as NA
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

#GlobalHealthObservatoryMetadata downloaded from https://apps.who.int/gho/data/node.metadata.COUNTRY?lang=en; minor edit to conform to country names in plasmiddf e.g. Bahamas -> The Bahamas; Gambia -> The Gambia; saved edited doc as GlobalHealthObservatoryMetadata_edit.tsv
#Additinal discrepancies are handled in code below, and described as follows:
#Raas Cabaad, Somalia ; change to Somalia
#Antarctica; change to '-'
#Svalbard and Jan Mayen; change to Norway
#United States Minor Outlying Islands; change to United States

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
#-     Taiwan Antarctica    Réunion    Mayotte Martinique 
#7284         93         58         10          3          2
unique(countryregiondf$incomes)
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
      income='Middle-income'
    }
    else if (income=='WB_HI') {
      income='High-income'
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

geographies<-factor(geographies, ordered = FALSE,levels = c("High-income", "Middle-income", "China", "United States", "EU", "other"))

###

#isolation source; refactored original factor (IsolationSource_unmergedfactorlevels), merging livestock and assigning (non-livestock) "agriculture" as other, to avoid low numbers
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
isolationsources<-factor(isolationsources, ordered = FALSE,levels = c("human", "livestock","other"))
table(isolationsources)

#replicon type: group by untyped/single/multi replicon
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
    reptypes<-c(reptypes,'single')
  }
  else {
    reptypes<-c(reptypes,'multi')
  }
}
reptypes<-factor(reptypes, ordered = FALSE,levels = c("untyped", "single", "multi"))

#taxa
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

taxa<-factor(taxa, ordered = FALSE,levels = c("Enterobacteriaceae", "Proteobacteria", "Firmicutes","other"))


###save predictor and outcome varaibles to file (need to create per resistance class outcome and otherresgenes predictor varaibles)

outcomeclasseslist<-vector('list',length(outcomeclasses))
numotherresgeneslist<-vector('list',length(outcomeclasses))
otherresgenesbinarylist<-vector('list',length(outcomeclasses))

for (i in 1:length(outcomeclasses)) {
  #outcome
  outcomeclass<-outcomeclasses[i]
  outcomeclassbinary<-as.factor(resgenedfbinary[,outcomeclass])
  outcomeclasseslist[[i]]<-outcomeclassbinary
  #other resistance genes - as binary and as count (by class)
  otherresgenedf<-resgenedf[,-which(colnames(resgenedf) %in% c('TotalResGenes',outcomeclass))]
  otherresgenedf[otherresgenedf>1]<-1
  numotherresgenes<-rowSums(otherresgenedf)
  otherresgenesbinary<-numotherresgenes
  otherresgenesbinary[numotherresgenes>1]<-1
  numotherresgeneslist[[i]]<-numotherresgenes
  otherresgenesbinarylist[[i]]<-as.factor(otherresgenesbinary)
}

outcomeclasseslistdf<-do.call(cbind.data.frame, outcomeclasseslist)
otherresgenesbinarylistdf<-do.call(cbind.data.frame, otherresgenesbinarylist)
numotherresgeneslistdf<-do.call(cbind.data.frame, numotherresgeneslist)
colnames(outcomeclasseslistdf)<-paste('outcome',outcomeclasses,sep='')
colnames(numotherresgeneslistdf)<-paste('numotherresgenes',outcomeclasses,sep='')
colnames(otherresgenesbinarylistdf)<-paste('otherresgenesbinary',outcomeclasses,sep='') #not used in modelling (used numotherresgeneslistdf)

#create final dataframe of predictors and response
accession<-plasmiddf$Accession
logplasmidsize<-log10(plasmidsize)
finaldf<-data.frame(accession,plasmidsize,logplasmidsize,coldateimputed,numis,isscore,numbacmet,numvfdb,numintegron,bacmetbinary,vfdbbinary,integronbinary,conj,geographies,isolationsources,reptypes,taxa)
finaldf<-cbind(finaldf,outcomeclasseslistdf,numotherresgeneslistdf,otherresgenesbinarylistdf)

#check if any values are NA
for (i in 1:ncol(finaldf)) {
  if (any(is.na(finaldf[,i]))==TRUE)
    print(colnames(finaldf)[i])
}


###exploratory analysis - histograms to determine truncation (=winsorising); cross-tabulations to ensure factor levels balance bias-variance trade-off

#cross-tabulate categorical predictors and outcome variables (for later calcualtion of unadjusted odds ratios)
factorvars<-c('bacmetbinary','conj','geographies','integronbinary','isolationsources','reptypes','taxa','vfdbbinary')
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


#additional exploratory analysis
if (plotexploratory==TRUE) {
  #histograms of all variables that were originally on continuous scale
  writefilepath='output_exploratory/NumericPredictors_histograms.pdf'
  pdf(writefilepath,10,15)
  par(mfrow=c(3,3))
  hist(finaldf$plasmidsize,main=NULL,xlab='Size (kb)')
  hist(finaldf$logplasmidsize,main=NULL,xlab='Size (log10 kb)')
  hist(finaldf$coldateimputed,main=NULL,xlab='Collection dates (with imputation using NCBI CreateDate)')
  hist(finaldf$numis,main=NULL,xlab='Number of ISs')
  hist(finaldf$isscore,main=NULL,xlab='ISs per 10 kb')
  hist(finaldf$numbacmet,main=NULL,xlab='Number of biocide/metal resistance genes')
  hist(finaldf$numvfdb,main=NULL,xlab='Number of virulence genes')
  hist(finaldf$numintegron,main=NULL,xlab='Number of complete integrons')
  mtext(gsub('%s',outcomeclass,'(Numeric predictors; some predictors will be modelled as binary'), side = 3, line = -2, outer = TRUE)
  dev.off()
  
  #histograms of other resgenes
  writefilepath='output_exploratory/NumericPredictors_otherresgenes_histograms.pdf'
  pdf(writefilepath,10,15)
  par(mfrow=c(3,3))
  for (i in 1:length(outcomeclasses)) {
    outcomeclass<-outcomeclasses[i]
    hist(finaldf[,gsub('%s',outcomeclass,'numotherresgenes%s')],xlab='Other resistance gene classes',main=outcomeclass)
  }
  dev.off()
  
  #other exploratory analysis, saved to resistance class folders
  for (i in 1:length(outcomeclasses)) {
    outcomeclass<-outcomeclasses[i]
    dir.create(gsub('%s',outcomeclass,'output_exploratory/%s'),showWarnings = FALSE)
    
    #relationship between continuous predictors and outcome
    outcomeclassbinary<-finaldf[,gsub('%s',outcomeclass,'outcome%s')]
    p1<-contpredictorplot(finaldf$plasmidsize,outcomeclassbinary,outcomeclass,'Plasmid size (kb)')
    p2<-contpredictorplot(finaldf$logplasmidsize,outcomeclassbinary,outcomeclass,'log10 plasmid size (kb)')
    p3<-contpredictorplot(finaldf$coldateimputed,outcomeclassbinary,outcomeclass,'Collection date (with imputation using NCBI CreateDate)')
    p4<-contpredictorplot(finaldf$numis,outcomeclassbinary,outcomeclass,'Number of ISs')
    p5<-contpredictorplot(finaldf$isscore,outcomeclassbinary,outcomeclass,'ISs per 10 kb')
    p6<-contpredictorplot(finaldf[,gsub('%s',outcomeclass,'numotherresgenes%s')],outcomeclassbinary,outcomeclass,'Number of other resistance gene classes')
    outputname<-'NumericPredictors_v_outcomeclass.pdf'
    writefilepath=gsub('%s',outcomeclass,'output_exploratory/%s/NumericPredictors_v_outcomeclass.pdf')
    pdf(writefilepath,15,10)
    print(plot_grid(p1,p2,p3,p4,p5,p6))
    dev.off()
    
    #relationship between categorical predictors and outcome (cross tabulation with proportions and odds ratios)
    proptablefilename<-gsub('%s',outcomeclass,'output_exploratory/%s/proptable.tsv')
    cat(gsub('%s',outcomeclass,'relationship between categorical predictors and outcome: %s binary'),file=proptablefilename,sep = '\n')
    writeproptablelist(myproptable(outcomeclassbinary,finaldf$bacmetbinary,'bacmetbinary'),proptablefilename)
    writeproptablelist(myproptable(outcomeclassbinary,finaldf$conj,'conj'),proptablefilename)
    writeproptablelist(myproptable(outcomeclassbinary,finaldf$geographies,'geographies'),proptablefilename)
    writeproptablelist(myproptable(outcomeclassbinary,finaldf$integronbinary,'integronbinary'),proptablefilename)
    writeproptablelist(myproptable(outcomeclassbinary,finaldf$isolationsources,'isolationsources'),proptablefilename)
    writeproptablelist(myproptable(outcomeclassbinary,finaldf[,gsub('%s',outcomeclass,'otherresgenesbinary%s')],'otherresgenesbinary'),proptablefilename)
    writeproptablelist(myproptable(outcomeclassbinary,finaldf$reptypes,'reptypes'),proptablefilename)
    writeproptablelist(myproptable(outcomeclassbinary,finaldf$taxa,'taxa'),proptablefilename)
    writeproptablelist(myproptable(outcomeclassbinary,finaldf$vfdbbinary,'vfdbbinary'),proptablefilename)
  }
}


###pre-modelling truncation; refactoring collection date to express as date since 5th percentile; refactoring log10plasmidsize to centre on 100kb

#truncate outliers for continuous variables plasmidsize (truncate both tails), numis (truncate right tail), coldateimputed (truncate left tail); numotherresgenes (truncate right tail)
finaldftrunc<-finaldf
plasmidsize95pct<-round(quantile(finaldftrunc$plasmidsize, probs = c(0.05, 0.95))[2])
finaldftrunc$plasmidsize[finaldftrunc$plasmidsize>plasmidsize95pct]<-plasmidsize95pct
plasmidsize5pct<-round(quantile(finaldftrunc$plasmidsize, probs = c(0.05, 0.95))[1])
finaldftrunc$plasmidsize[finaldftrunc$plasmidsize<plasmidsize5pct]<-plasmidsize5pct
finaldftrunc$logplasmidsize<-log10(finaldftrunc$plasmidsize) #log10plasmidsize derived from truncated plasmidsize
numis95pct<-round(quantile(finaldftrunc$numis, probs = c(0.05, 0.95))[2])
finaldftrunc$numis[finaldftrunc$numis>numis95pct]<-numis95pct

isscore95pct<-round(quantile(finaldftrunc$isscore, probs = c(0.05, 0.95)),2)[2]
finaldftrunc$isscore[finaldftrunc$isscore>isscore95pct]<-isscore95pct
#hist(finaldftrunc$isscore)

coldateimputed5pct<-round(quantile(finaldftrunc$coldateimputed, probs = c(0.05, 0.95))[1])
finaldftrunc$coldateimputed[finaldftrunc$coldateimputed<coldateimputed5pct]<-coldateimputed5pct

numotherresgenestrunccutoffs<-list()
for (outcomeclass in outcomeclasses) {
  numotherresgenes95pct<-round(quantile(finaldftrunc[,gsub('%s',outcomeclass,'numotherresgenes%s')], probs = c(0.05, 0.95))[2])
  finaldftrunc[,gsub('%s',outcomeclass,'numotherresgenes%s')][finaldftrunc[,gsub('%s',outcomeclass,'numotherresgenes%s')]>numotherresgenes95pct]<-numotherresgenes95pct
  numotherresgenestrunccutoffs[[outcomeclass]]<-numotherresgenes95pct
}


#re-express collection date as date since 5% percentile (1994)
finaldftrunc$coldateimputed<-finaldftrunc$coldateimputed-round(as.numeric(coldateimputed5pct))

#re-express log10plasmidsize as log10plasmidsize - 4 (i.e. centred on 10kb since 0kb is not meaningful)
finaldftrunc$logplasmidsize<-finaldftrunc$logplasmidsize - 4


#write to file
write.table(finaldftrunc,'data/Appendix_D_Table_6_splitbetalactam_truncated.tsv',col.names=TRUE,row.names=FALSE,sep='\t',quote=FALSE)

#write truncation thresholds to file
trunccutofffilename<-'data/truncationcutoffs.tsv'
cat('',file=trunccutofffilename,append=FALSE)
cat(gsub('%s',plasmidsize95pct,'plasmidsize95pct: %s\n'),file=trunccutofffilename,append=TRUE)
cat(gsub('%s',plasmidsize5pct,'plasmidsize5pct: %s\n'),file=trunccutofffilename,append=TRUE)
cat(gsub('%s',numis95pct,'numis95pct: %s\n'), file=trunccutofffilename,append=TRUE)  
cat(gsub('%s',isscore95pct,'isscore95pct: %s\n'), file=trunccutofffilename,append=TRUE) 
cat(gsub('%s',coldateimputed5pct,'coldateimputed5pct: %s\n'),file=trunccutofffilename,append=TRUE)
for (outcomeclass in outcomeclasses) {
  numotherresgenes95pct<-numotherresgenestrunccutoffs[[outcomeclass]]
  cat(gsubfn('%1|%2',list('%1'=outcomeclass,'%2'=numotherresgenes95pct),'%1 numotherresgenes95pct: %2\n'),file=trunccutofffilename,append=TRUE)
}


#With truncated dataset, flag potential issues with collinearity using association metrics and VIFs of glm models
print('exploring collinearity...')
finaldftrunccopy<-finaldftrunc
for (binaryvar in c('integronbinary','bacmetbinary','vfdbbinary')) {
  finaldftrunccopy[,binaryvar]<-as.numeric(as.character(finaldftrunccopy[,binaryvar]))
}

#pairs plot for spearman's rank correlation between numeric variables
numotherresgenesvec<-vector()
for (outcomeclass in outcomeclasses) {
  numotherresgenes<-gsub('%s',outcomeclass,'numotherresgenes%s')
  numotherresgenesvec<-c(numotherresgenesvec,numotherresgenes)
}

if (plotpairs==TRUE) {
  pairsfilepath=gsub('%s',outcomeclass,'output_exploratory/cor_predictors_pairsplot_scatter.png')
  png(pairsfilepath,4000,4000)  #png rather than vector graphic is more efficient for plots with many points
  #pairs(finaldftrunccopy[,c('logplasmidsize','numis', 'isscore', numotherresgenesvec,'coldateimputed','integronbinary','bacmetbinary','vfdbbinary')],upper.panel=panel.smooth,diag.panel = panel.hist,lower.panel = panel.cor,cex.axis=4)
  pairs(finaldftrunccopy[,c('logplasmidsize','numis','isscore',numotherresgenesvec,'coldateimputed')],upper.panel=panel.smooth,diag.panel = panel.hist,lower.panel = panel.cor,cex.axis=4)
  dev.off()
  pairsfilepath=gsub('%s',outcomeclass,'output_exploratory/cor_predictors_pairsplot.png')
  png(pairsfilepath,4000,4000)  #png rather than vector graphic is more efficient for plots with many points
  #pairs(finaldftrunccopy[,c('logplasmidsize','numis', 'isscore', numotherresgenesvec,'coldateimputed','integronbinary','bacmetbinary','vfdbbinary')],upper.panel=NULL,diag.panel = panel.hist,lower.panel = panel.cor,cex.axis=4)
  pairs(finaldftrunccopy[,c('logplasmidsize','numis', 'isscore', numotherresgenesvec,'coldateimputed')],upper.panel=NULL,diag.panel = panel.hist,lower.panel = panel.cor,cex.axis=4)
  dev.off()
}


#Matrix outputs: spearman's correlation; Kramer's V for association between categorical variables (0=no association; 1=perfect association); Kruskal-Wallis H test statistic for association between continuoius and categorical

#Numeric: logplasmidsize, numis, isscore, numotherresgenes, coldateimputed
#Categorical:
##Binary: integron, bacmet, vfdb
##Nominal: conj, reptypes, taxa, geographies, isolationsources

#numeric-numeric (Spearman's correlation)
#(finaldftrunc)
numericvars<-c('logplasmidsize','numis','isscore',numotherresgenesvec,'coldateimputed')
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

#categorical-categorical (Cramers V)
categoricalvars<-c('integronbinary','bacmetbinary','vfdbbinary','conj','reptypes','taxa','geographies','isolationsources')

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

#numeric-categorical (Kruskal-Wallis H)

#out<-kruskal.test(numotherresgenesaminoglycoside ~ integronbinary, data=finaldftrunc)

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

finaldftrunc %>% kruskal_effsize(numotherresgenesaminoglycoside ~ integronbinary)


kruskalH_matrix<-as.data.frame(kruskalH_matrix)
write.table(kruskalH_matrix,file='output_exploratory/cor_kruskalH_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)
kruskalpvalue_matrix<-as.data.frame(kruskalpvalue_matrix)
write.table(kruskalpvalue_matrix,file='output_exploratory/cor_kruskalpvalue_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)
kruskaletasqrd_matrix<-as.data.frame(kruskaletasqrd_matrix)
write.table(kruskaletasqrd_matrix,file='output_exploratory/cor_kruskaletasqrd_matrix.tsv',sep='\t',row.names = TRUE,col.names = NA)


#numeric-binary (kruskalH isn't really an association metric; point-biserial is, but only applies to binary categorical)
#https://www.statology.org/point-biserial-correlation-r/
#Can't use point-biserial because normality assumption is violated; use Spearman's instead

binaryvars<-c('integronbinary','bacmetbinary','vfdbbinary')

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


##build glm and calculate VIFs
# viffilename<-'output_exploratory/glm_VIFs.tsv'
# viffilename2<-'output_exploratory/glm_VIFs_logplasmidsizeremoved.tsv'
# cat('VIF output\n',file=viffilename,sep = '\n',append = FALSE)
# cat('VIF output\n',file=viffilename2,sep = '\n',append = FALSE)
# for (outcomeclass in outcomeclasses) {
#   ##initial model including logplasmidsize
#   frm <- formula(gsub('%s',outcomeclass,'outcome%s~logplasmidsize+numis+numotherresgenes%s+coldateimputed+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
#   myglm<-glm(frm,family='binomial',data=finaldftrunc)
#   cat(gsub('%s',outcomeclass,'\n###VIFs for outcome: %s binary'),file=viffilename,sep = '\n',append = TRUE)
#   vifouts<-capture.output(car::vif(myglm))
#   for (vifout in vifouts) {
#     vifout<-gsub("*\\s+", "\t", vifout)
#     cat(vifout,sep='\n',append = TRUE,file = viffilename)
#   }
#   ##after identifying potential collinearity, removed logplamsidsize; re-run vif calculations
#   frm <- formula(gsub('%s',outcomeclass,'outcome%s~numis+numotherresgenes%s+coldateimputed+integronbinary+bacmetbinary+conj+reptypes+taxa+vfdbbinary+geographies+isolationsources'))
#   myglm<-glm(frm,family='binomial',data=finaldftrunc)
#   cat(gsub('%s',outcomeclass,'\n###VIFs for outcome: %s binary'),file=viffilename2,sep = '\n',append = TRUE)
#   vifouts<-capture.output(car::vif(myglm))
#   for (vifout in vifouts) {
#     vifout<-gsub("*\\s+", "\t", vifout)
#     cat(vifout,sep='\n',append = TRUE,file = viffilename2)
#   }
# }


# 
# #regression tree analysis after truncation
# library(rpart)
# library(rpart.plot)
# #head(finaldftrunc[,c(4,5,6,8,12,14,15,16,17,19,20,23)])
# for (outcomeclass in outcomeclasses) {
#   writefilepath=gsub('%s',outcomeclass,'output_exploratory/%s/rparttree.pdf')
#   #pdf(writefilepath,10,10)
#   fit<-rpart(finaldf[,gsub('%s',outcomeclass,'outcome%s')]~.,data=finaldftrunc[,c('logplasmidsize','coldateimputed','numis','bacmetbinary','vfdbbinary','integronbinary','conj','geographies','isolationsources','reptypes','taxa',gsub('%s',outcomeclass,'numotherresgenes%s'))],method='class',model=T)
#   rpart.plot(fit)
#   #dev.off()
# }
# 
