library(epitools)
setwd("C:/Users/alexo/Dropbox/Oxford_2015/thesis/CORRECTIONS/Chapter_5/resistance_prediction")

oddsdf<-read.table('output_exploratory/crosstabulation.tsv',as.is = TRUE,sep='\t',header=TRUE)

resclasses<-unique(oddsdf$ResistanceClass)
factorvars<-unique(oddsdf$FactorVariable)

ors<-vector()
lowercis<-vector()
uppercis<-vector()
logors<-vector()
loglowercis<-vector()
loguppercis<-vector()
probors<-vector()
problowercis<-vector()
probuppercis<-vector()
pvalueschi<-vector()
pvaluesfisher<-vector()
for (i in 1:length(resclasses)) {
  for (j in 1:length(factorvars)) {
    resvardata<-oddsdf[oddsdf$ResistanceClass==resclasses[i] & oddsdf$FactorVariable==factorvars[j],]
    for (r in 1:nrow(resvardata)) {
      if (r==1) {
        baselinenonresistance<-resvardata[r,'Count_NonResistancePlasmids']
        baselineresistance<-resvardata[r,'Count_ResistancePlasmids']
        ors<-c(ors,'baseline')
        lowercis<-c(lowercis,'baseline')
        uppercis<-c(uppercis,'baseline')
        pvalueschi<-c(pvalueschi,'baseline')
        pvaluesfisher<-c(pvaluesfisher,'baseline')
        logors<-c(logors,'baseline')
        loglowercis<-c(loglowercis,'baseline')
        loguppercis<-c(loguppercis,'baseline')
        probors<-c(probors,'baseline')
        problowercis<-c(problowercis,'baseline')
        probuppercis<-c(probuppercis,'baseline')
      } else {
        exposednonresistance<-resvardata[r,'Count_NonResistancePlasmids']
        exposedresistance<-resvardata[r,'Count_ResistancePlasmids']
        ORtable<-matrix(c(baselinenonresistance,exposednonresistance,baselineresistance,exposedresistance),nrow = 2, ncol = 2)
        ORoutput<-oddsratio.wald(ORtable)
        #check for complete separation
        if (any(ORtable==0)) { #complete separation - fill with NA
          ors<-c(ors,NA)
          lowercis<-c(lowercis,NA)
          uppercis<-c(uppercis,NA)
          pvalueschi<-c(pvalueschi,NA)
          pvaluesfisher<-c(pvaluesfisher,NA)
          logors<-c(logors,NA)
          loglowercis<-c(loglowercis,NA)
          loguppercis<-c(loguppercis,NA)
          probors<-c(probors,NA)
          probuppercis<-c(probuppercis,NA)
          problowercis<-c(problowercis,NA)
        } else {
          or<-ORoutput$measure['Exposed2','estimate']
          lowerci<-ORoutput$measure['Exposed2','lower']
          upperci<-ORoutput$measure['Exposed2','upper']
          pvaluechi<-ORoutput$p.value['Exposed2','chi.square']
          pvaluefisher<-ORoutput$p.value['Exposed2','fisher.exact']
          ors<-c(ors,or)
          lowercis<-c(lowercis,lowerci)
          uppercis<-c(uppercis,upperci)
          pvalueschi<-c(pvalueschi,pvaluechi)
          pvaluesfisher<-c(pvaluesfisher,pvaluefisher)
          logors<-c(logors,log(or))
          loglowercis<-c(loglowercis,log(lowerci))
          loguppercis<-c(loguppercis,log(upperci))
          probors<-c(probors,plogis(log(or)))
          problowercis<-c(problowercis,plogis(log(lowerci)))
          probuppercis<-c(probuppercis,plogis(log(upperci)))
        }
      }
    }
  }
}

oddsdf$OddsRatio<-ors
oddsdf$Lower95CI<-lowercis
oddsdf$Upper95CI<-uppercis
oddsdf$logOddsRatio<-logors
oddsdf$logLower95CI<-loglowercis
oddsdf$logUpper95CI<-loguppercis
oddsdf$ProbabilityScale<-probors
oddsdf$ProbabilityScaleLower95CI<-problowercis
oddsdf$ProbabilityScaleUpper95CI<-probuppercis
oddsdf$`p-valueChi`<-pvalueschi
oddsdf$`p-valueFisher`<-pvaluesfisher

#re-order factor levels
geographiesdf<-oddsdf[oddsdf$FactorVariable=='geographies',]
geographiesdf<-do.call(rbind,lapply(split(geographiesdf,geographiesdf$ResistanceClass),function(x) x[match(c("High-income","Middle-income","EU","China","United States","other"),x$FactorLevel),]))
oddsdf[oddsdf$FactorVariable=='geographies',]<-geographiesdf

#save
write.table(oddsdf,file='output_unadjusted/unadjustedodds.tsv',quote=FALSE,sep='\t',row.names = FALSE)

