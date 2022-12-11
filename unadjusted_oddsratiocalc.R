pacman::p_load(epitools, tidyverse, readxl)
source('functions.R')
dir.create('output_unadjusted', showWarnings = FALSE)
dir.create(file.path('output_unadjusted','firmicutes_analysis'), showWarnings = FALSE)


# ---------------------------
# load data: for main analysis load previously calculated crosstabulation.tsv; for sub-analyses generate new cross-tabulations

# main analysis
main_crosstabulation<-read.table('output_exploratory/crosstabulation.tsv',as.is = TRUE,sep='\t',header=TRUE)

# analysis of macrolide resistance across taxa and most frequent Firmicutes species
finaldftrunc <- read.table('data/plasmiddf_transformed.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")
plasmiddf <- read_excel("data/Data_S1.xlsx", sheet = 'G')
plasmiddf <- plasmiddf %>% mutate(across(everything(), ~ str_replace_all(.x, "^\'|\'$", "")))
colnames(plasmiddf) <- colnames(plasmiddf) %>% str_replace_all("^\'|\'$", "")
plasmiddf <- plasmiddf %>% mutate(InFinalDataset = as.logical(InFinalDataset)) %>% filter(InFinalDataset == TRUE) %>% select(Accession, Species)
finaldftrunc <- finaldftrunc %>% left_join(plasmiddf, by = 'Accession')

crosstablelist <- list()
crosstabledf <- as.data.frame(crosstabulate(finaldftrunc,'macrolide','HostTaxonomy'))
crosstablelist[['all_taxa']] <- crosstabledf

top_firmicutes_species <- c('Staphylococcus aureus', 'Bacillus thuringiensis', 'Lactobacillus plantarum', 'Enterococcus faecium', 'Lactococcus lactis', 'Bacillus cereus')
for (firmicutes_species in top_firmicutes_species) {
  plasmids_subset <- finaldftrunc %>% filter(HostTaxonomy == 'Firmicutes' & Species == firmicutes_species) %>% mutate(HostTaxonomy = ifelse(HostTaxonomy == 'Firmicutes', str_c('Firmicutes: ', firmicutes_species), HostTaxonomy))
  crosstabledf <- as.data.frame(crosstabulate(plasmids_subset,'macrolide','HostTaxonomy'))
  crosstablelist[[firmicutes_species]] <- crosstabledf
}

firmicutes_crosstabulation <- as.data.frame(do.call(rbind,crosstablelist))
firmicutes_crosstabulation <- firmicutes_crosstabulation %>% mutate(ResistanceClass = as.character(ResistanceClass), FactorVariable = as.character(FactorVariable), FactorLevel = as.character(FactorLevel), Count_NonResistancePlasmids = as.integer(as.character(Count_NonResistancePlasmids)), Count_ResistancePlasmids = as.integer(as.character(Count_ResistancePlasmids)))


# ---------------------------
# calculate odds from crosstabulations
calculate_odds <- function(crosstabledf) {
  
  resclasses<-unique(crosstabledf$ResistanceClass)
  factorvars<-unique(crosstabledf$FactorVariable)
  
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
      resvardata<-crosstabledf[crosstabledf$ResistanceClass==resclasses[i] & crosstabledf$FactorVariable==factorvars[j],]
      for (r in 1:nrow(resvardata)) {
        if (r==1) {
          baselinenonresistance<-resvardata[r,'Count_NonResistancePlasmids']
          baselineresistance<-resvardata[r,'Count_ResistancePlasmids']
          ors<-c(ors,'reference')
          lowercis<-c(lowercis,'reference')
          uppercis<-c(uppercis,'reference')
          pvalueschi<-c(pvalueschi,'reference')
          pvaluesfisher<-c(pvaluesfisher,'reference')
          logors<-c(logors,'reference')
          loglowercis<-c(loglowercis,'reference')
          loguppercis<-c(loguppercis,'reference')
          probors<-c(probors,'reference')
          problowercis<-c(problowercis,'reference')
          probuppercis<-c(probuppercis,'reference')
        } else {
          exposednonresistance<-resvardata[r,'Count_NonResistancePlasmids']
          exposedresistance<-resvardata[r,'Count_ResistancePlasmids']
          ORtable<-matrix(c(baselinenonresistance,exposednonresistance,baselineresistance,exposedresistance),nrow = 2, ncol = 2)
          ORoutput<-oddsratio.wald(ORtable)
          #check for complete separation
          if (any(ORtable==0)) {  # complete separation - fill with NA
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
  
  crosstabledf$OddsRatio<-ors
  crosstabledf$Lower95CI<-lowercis
  crosstabledf$Upper95CI<-uppercis
  crosstabledf$logOddsRatio<-logors
  crosstabledf$logLower95CI<-loglowercis
  crosstabledf$logUpper95CI<-loguppercis
  crosstabledf$ProbabilityScale<-probors
  crosstabledf$ProbabilityScaleLower95CI<-problowercis
  crosstabledf$ProbabilityScaleUpper95CI<-probuppercis
  crosstabledf$`p-valueChi`<-pvalueschi
  crosstabledf$`p-valueFisher`<-pvaluesfisher
  
  # sort by outcome class
  crosstabledf <- crosstabledf[order(crosstabledf$ResistanceClass,crosstabledf$FactorVariable),]
  return(crosstabledf)
}


main_oddsdf <- calculate_odds(main_crosstabulation)
firmicutes_oddsdf <- calculate_odds(firmicutes_crosstabulation)


# save
write.table(main_oddsdf, file='output_unadjusted/unadjustedodds.tsv',quote=FALSE,sep='\t',row.names = FALSE)
write.table(firmicutes_oddsdf, file='output_unadjusted/firmicutes_analysis/unadjustedodds.tsv',quote=FALSE,sep='\t',row.names = FALSE)


