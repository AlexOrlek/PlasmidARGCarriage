pacman::p_load(tidyverse,scales,RColorBrewer,rwantshue,ggtext)
theme_set(theme_bw())

oddsdf<-read.table('output_unadjusted/unadjustedodds.tsv',as.is = TRUE,sep='\t',header=TRUE)
baselineindices<-oddsdf$OddsRatio=='reference'
baselineindices[is.na(baselineindices)]<-FALSE
oddsdf<-oddsdf[!baselineindices,]
oddsdf$logOddsRatio<-as.numeric(oddsdf$logOddsRatio)
oddsdf$logLower95CI<-as.numeric(oddsdf$logLower95CI)
oddsdf$logUpper95CI<-as.numeric(oddsdf$logUpper95CI)
resclasses<-unique(oddsdf$ResistanceClass)

brewerpal1<-brewer.pal(8,'Set1')
brewerpal<-c(brewerpal1,'#e0bb6d','#90EE90','#add8e6')
brewerpal[6]<-'#ffd700'
#plot(1:11,1:11,col=brewerpal,pch=16,cex=6)

if (length(resclasses)>11) {
  scheme <- iwanthue(seed = 42, force_init = TRUE)
  brewerpal <- scheme$hex(length(resclasses))
}


# split by factor variables, facet by factor level, resistance class on x axis
oddsdfsplit<-split(oddsdf,oddsdf$FactorVariable)

outputnames<-c('Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
ggtitles<-c('Integron presence\nreference: absence','Biocide/metal resistance gene presence\nreference: absence','Conjugative system\nreference: non-mobilisable','Replicon carriage\nreference: untyped','Host taxonomy\nreference: Enterobacteriaceae','Virulence gene presence\nreference: absence','Geographic location\nreference: high-income','Isolation source\nreference: human')
factorlevels<-list(c('presence'),c('presence'),c('mobilisable','conjugative'),c('single-replicon','multi-replicon'),c("Proteobacteria_non-Enterobacteriaceae", "Firmicutes","other"),c('presence'),c("China", "United States", "other","middle-income", "EU & UK"),c("livestock","other"))
factorlabels<-list(c('presence'),c('presence'),c('mobilisable','conjugative'),c('single-replicon','multi-replicon'),c("Proteobacteria (non-Enterobacteriaceae)", "Firmicutes","other"),c('presence'),c("China", "United States", "other","middle-income", "EU & UK"),c("livestock","other"))
outcomeclasses<-c('aminoglycoside','sulphonamide','tetracycline','phenicol','macrolide','trimethoprim','ESBL', 'carbapenem','quinolone','colistin')

numpanels<-c(1,1,2,2,3,1,5,2)
numcols=3
if (numcols==5) {
  width_multiplicationfactor<-c(1.159, 1.159, 2.120, 2.120, 3.087, 1.159, 5.145, 2.120)
  height_multiplicationfactor<-c(1,1,1,1,1,1,1,1)
} else {  # ncol==3
  width_multiplicationfactor<-c(1.159, 1.159, 2.120, 2.120, 3.087, 1.159, 3.087, 2.120)
  height_multiplicationfactor<-c(1,1,1,1,1,1,1.61,1)
}


width=3.2
height=4.5

for (i in 1:length(outputnames)) {
  outputname<-outputnames[i]
  print(outputname)
  #get data and fix factor levels and labelling
  factorvardf<-oddsdfsplit[[outputnames[i]]]
  factorvardf$FactorLevel<-factor(factorvardf$FactorLevel, ordered = FALSE,levels = factorlevels[[i]])
  labels<-factorlabels[[i]]
  names(labels)<-factorlevels[[i]]
  #construct plot
  myggtitle<-ggtitles[i]
  lowerlim<--5
  upperlim<-5
  if (outputname=='Integron') {
    lowerlim<--2
    upperlim<-8
  }
  if (outputname=='HostTaxonomy') {
    lowerlim<--6
    upperlim<-4
  }
  p<-ggplot(factorvardf,aes(x=ResistanceClass,y=logOddsRatio,color=ResistanceClass)) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + geom_errorbar(aes(ymin=logLower95CI,ymax=logUpper95CI),colour='grey 42',linetype=1,width=0.5,size=0.4) + geom_point(size=2)
  p<-p + facet_wrap(~ FactorLevel, labeller=labeller(FactorLevel=labels), as.table=FALSE,ncol=numcols)
  p<-p + theme(legend.position="none",axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
  p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) + scale_x_discrete(limits=outcomeclasses) + scale_colour_manual(values=brewerpal[match(sort(outcomeclasses),outcomeclasses)])
  #plot
  mywidth<-width*width_multiplicationfactor[i]
  myheight<-height*height_multiplicationfactor[i]
  pdf(gsub('%s',outputnames[i],'output_unadjusted/%s_logoddsscale_faceted_unadjusted.pdf'),mywidth,myheight)
  print(p)
  dev.off()
}



# plot firmicutes / macrolide odds

oddsdf_firmicutes<-read.table('output_unadjusted/firmicutes_analysis/unadjustedodds.tsv',as.is = TRUE,sep='\t',header=TRUE)
baselineindices<-oddsdf_firmicutes$OddsRatio=='reference'
baselineindices[is.na(baselineindices)]<-FALSE
oddsdf_firmicutes<-oddsdf_firmicutes[!baselineindices,]
oddsdf_firmicutes$logOddsRatio<-as.numeric(oddsdf_firmicutes$logOddsRatio)
oddsdf_firmicutes$logLower95CI<-as.numeric(oddsdf_firmicutes$logLower95CI)
oddsdf_firmicutes$logUpper95CI<-as.numeric(oddsdf_firmicutes$logUpper95CI)

oddsdf_firmicutes <- oddsdf_firmicutes %>% filter(!FactorLevel %in% c('other', 'Proteobacteria_non-Enterobacteriaceae')) %>% mutate(FactorLevel = recode(FactorLevel, 'Firmicutes' = 'all Firmicutes')) %>% mutate(FactorLevel = str_replace_all(FactorLevel, "Firmicutes: ", "")) %>% mutate(FactorLevel = ifelse(FactorLevel == 'all Firmicutes', FactorLevel, str_c('*',FactorLevel,'*')))
oddsdf_firmicutes$FactorLevel <- factor(oddsdf_firmicutes$FactorLevel, levels = c('all Firmicutes', head(sort(oddsdf_firmicutes$FactorLevel),-1)))

p<-ggplot(oddsdf_firmicutes,aes(x=FactorLevel,y=logOddsRatio,color="#FF7F00")) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + geom_errorbar(aes(ymin=logLower95CI,ymax=logUpper95CI),colour='grey 42',linetype=1,width=0.5,size=0.4) + geom_point(size=2)

p<-p + theme(legend.position="none",axis.text.x = element_markdown(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle("Host taxonomy: Firmicutes; macrolide resistance carriage\nreference: Enterobacteriaceae")
p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(-5,2))
pdf('output_unadjusted/firmicutes_analysis/firmicutes_macrolide_unadjusted.pdf', width = 7, height = 7/1.159)
print(p)
dev.off()


