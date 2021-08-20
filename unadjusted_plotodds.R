pacman::p_load(ggplot2,scales,RColorBrewer)
theme_set(theme_bw())

oddsdf=read.table('output_unadjusted/unadjustedodds.tsv',as.is = TRUE,sep='\t',header=TRUE)
baselineindices<-oddsdf$OddsRatio=='baseline'
baselineindices[is.na(baselineindices)]<-FALSE
oddsdf<-oddsdf[!baselineindices,]
oddsdf$logOddsRatio<-as.numeric(oddsdf$logOddsRatio)
oddsdf$logLower95CI<-as.numeric(oddsdf$logLower95CI)
oddsdf$logUpper95CI<-as.numeric(oddsdf$logUpper95CI)
resclasses<-unique(oddsdf$ResistanceClass)

brewerpal1<-brewer.pal(8,'Set1')
brewerpal1[6]<-'#CDCD00'
#plot(1:8,1:8, col=brewerpal1, pch=16,cex=5)
betalactampal<-brewer.pal(9,'Blues')[c(3,6,9)]
#plot(1:3,1:3,col=betalactampal,pch=16,cex=5)
brewerpal2<-c(brewerpal1[1],betalactampal,brewerpal1[3:length(brewerpal1)])
#plot(1:10,1:10, col=brewerpal2, pch=16,cex=5)
numbetalactam<-sum(grepl('^betalactam',resclasses))

if (numbetalactam==3) {
  brewerpal<-brewerpal2
  stopifnot(length(resclasses)==10)
} else if (numbetalactam==1) {
  brewerpal<-brewerpal1
  stopifnot(length(resclasses)==8)
} else {
  stop('Error: unexpected number of beta-lactam resistance subclasses ***')
}


# split by factor variables, facet by factor level, resistance class on x axis
oddsdfsplit<-split(oddsdf,oddsdf$FactorVariable)

outputnames<-c('Integron','BiocideMetalResistance','ConjugativeSystem','RepliconCarriage','HostTaxonomy','Virulence','GeographicLocation','IsolationSource')
ggtitles<-c('Integron presence\nbaseline: absence','Biocide/metal resistance gene presence\nbaseline: absence','Conjugative system\nbaseline: non-mobilisable','Replicon type carriage\nbaseline: untyped','Host taxonomy\nbaseline: Enterobacteriaceae','Virulence gene presence\nbaseline: absence','Geographic location\nbaseline: high-income','Isolation source\nbaseline: human')
factorlevels<-list(c('1'),c('1'),c('mobilisable','conjugative'),c('single-replicon','multi-replicon'),c("Proteobacteria_other", "Firmicutes","other"),c('1'),c("China", "United States", "other","middle-income", "EU"),c("livestock","other"))
factorlabels<-list(c('presence'),c('presence'),c('mobilisable','conjugative'),c('single-replicon','multi-replicon'),c("Proteobacteria_other", "Firmicutes","other"),c('presence'),c("China", "United States", "other","middle-income", "EU"),c("livestock","other"))

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
  #get colours
  if (numpanels[i]>1) {
    mypal<-vector()
    for (palcol in brewerpal) {
      mypal<-c(mypal,rep(palcol,numpanels[i]))
    }
  } else {
    mypal<-brewerpal
  }
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
    lowerlim<--8
    upperlim<-2
  }
  p<-ggplot(factorvardf,aes(x=ResistanceClass,y=logOddsRatio)) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + geom_errorbar(aes(ymin=logLower95CI,ymax=logUpper95CI),colour='grey 42',linetype=1,width=0.5,size=0.4) + geom_point(colour=mypal,size=2)
  p<-p + facet_wrap(~ FactorLevel, labeller=labeller(FactorLevel=labels), as.table=FALSE,ncol=numcols)
  p<-p + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
  p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim))
  #plot
  mywidth<-width*width_multiplicationfactor[i]
  myheight<-height*height_multiplicationfactor[i]
  pdf(gsub('%s',outputnames[i],'output_unadjusted/%s_logoddsscale_faceted_unadjusted.pdf'),mywidth,myheight)
  print(p)
  dev.off()
}
