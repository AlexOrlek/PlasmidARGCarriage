pacman::p_load(tidyverse, RColorBrewer, reshape2)
theme_set(theme_bw())

# Comparison of parametric term log-odds ratios across models

# ---------------------------
# compare adjusted vs unadjusted models
unadjusteddf<-read.table('output_unadjusted/unadjustedodds.tsv',as.is = TRUE,sep='\t',header=TRUE)
adjusteddf<-read.table('output_adjusted/mainmodel/adjustedodds.tsv',as.is = TRUE,sep='\t',header=TRUE)

# remove 'reference' rows from unadjusted
unadjusteddf<-unadjusteddf[unadjusteddf$logOddsRatio!='reference' | is.na(unadjusteddf$logOddsRatio),]

# re-order by factor level
unadjusteddf<-unadjusteddf[order(unadjusteddf$FactorVariable),]
adjusteddf<-adjusteddf[order(adjusteddf$FactorVariable),]

# append unadjusted and unadjusted comparison to adjusteddf
adjusteddf <- adjusteddf %>% inner_join(unadjusteddf, by = c('ResistanceClass','FactorVariable','FactorLevel'), suffix = c('','_unadjusted')) %>% select(-c('Count_NonResistancePlasmids', 'Count_ResistancePlasmids')) %>% rename(p.valueChi_unadjusted = p.valueChi, p.valueFisher_unadjusted = p.valueFisher) %>% mutate(OddsRatio_difference_vs_unadjusted = as.numeric(OddsRatio) - as.numeric(OddsRatio_unadjusted), logOddsRatio_difference_vs_unadjusted = as.numeric(logOddsRatio) - as.numeric(logOddsRatio_unadjusted))

# save
write.table(adjusteddf,file='output_adjusted/mainmodel/adjustedodds_vs_unadjusted.tsv',row.names = FALSE,col.names = TRUE,sep='\t')


# ---------------------------
# compare adjusted alternative models vs 1) unadjusted; 2) adjusted (reference "mainmodel")

alternativemodelnames<-c('minus_NumOtherResistanceClasses', 'minus_Integron', 'minus_BiocideMetalResistance', 'minus_RepliconCarriage', 'minus_HostTaxonomy', 'minus_log10PlasmidSize', 'minus_RepliconCarriage_NumOtherResistanceClasses', 'minus_BiocideMetalResistance_NumOtherResistanceClasses', 'minus_Integron_NumOtherResistanceClasses', 'minus_BiocideMetalResistance_Integron', 'minus_3associatedfactorsofConjugativeSystem','minus_6associatedfactorsofConjugativeSystem')
removedparametricvariables<-list(NA,c('Integron'),c('BiocideMetalResistance'),c('RepliconCarriage'),c('HostTaxonomy'),NA,c('RepliconCarriage'),c('BiocideMetalResistance'),c('Integron'),c('BiocideMetalResistance','Integron'),c('HostTaxonomy'),c('HostTaxonomy','RepliconCarriage','Integron'))

for (i in 1:length(alternativemodelnames)) {
  modelname<-alternativemodelnames[i]
  removedvars<-removedparametricvariables[[i]]
  print(modelname)
  #create output directory
  dir.create(file.path('output_adjusted',modelname,'parametric_comparison'), showWarnings = FALSE)
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


# ---------------------------
# Create combined plots of log-odds from alternative adjusted model, main adjusted model, unadjusted analysis

# ---------------------------
# functions and parameters
outcomeclasses<-c('aminoglycoside','sulphonamide','tetracycline','phenicol','macrolide','trimethoprim','ESBL', 'carbapenem','quinolone','colistin')

ncols<-3  # for faceted figures

recolourdf<-function(df,pal=brewerpal) {
  # plot main model with filled coloured points, alternative model with empty coloured points, unadjusted with grey empty points
  df[['colour']][df[['model']]=='logOddsRatio_mainmodel']<-pal  # colour main model
  df[['shape']][df[['model']]=='logOddsRatio_unadjusted']<-4
  df[['shape']][df[['model']]=='logOddsRatio']<-2
  return(df)
}

setparameters<-function(outputname) {
  lowerlim<--5
  upperlim<-5
  if (outputname=='Integron') {
    lowerlim<--2
    upperlim<-8
  }
  if (outputname=='HostTaxonomy') {
    lowerlim<--8
    upperlim<-4
  }
  width=3.2
  height=4.5
  width_multiplicationfactor<-1.159
  height_multiplicationfactor<-1
  if (outputname %in% c('ConjugativeSystem','RepliconCarriage','IsolationSource')) {
    width_multiplicationfactor<-2.120
  }
  if (outputname %in% c('HostTaxonomy','GeographicLocation')) {
    width_multiplicationfactor<-3.087
  }
  if (outputname=='GeographicLocation') {
    height_multiplicationfactor<-1.61
  }
  mywidth<-width*width_multiplicationfactor
  myheight<-height*height_multiplicationfactor
  return(c(lowerlim,upperlim,mywidth,myheight))
}

brewerpal1<-brewer.pal(8,'Set1')
brewerpal<-c(brewerpal1,'#e0bb6d','#90EE90','#add8e6')
brewerpal[6]<-'#ffd700'
brewerpal<-brewerpal[1:length(outcomeclasses)]

labelfunc<-function(x) {
  #relabels facet grid labels
  x<-str_replace(x,'Proteobacteria_non-Enterobacteriaceae', 'Proteobacteria (non-Enterobacteriaceae)')
  return(x)
}


# ---------------------------
# plot log-odds across models

# HostTaxonomy: alternative model with RepliconCarriage and NumOtherResistanceClasses removed
myggtitle<-'Host taxonomy\nreference: Enterobacteriaceae'
adjustedalternativedf<-read.table("output_adjusted/minus_RepliconCarriage_NumOtherResistanceClasses/adjustedodds_vs_mainmodel.tsv",sep='\t',header=TRUE,as.is=TRUE)
factorvardf<-adjustedalternativedf[adjustedalternativedf$FactorVariable=='HostTaxonomy',]
factorvardf<-factorvardf[order(match(factorvardf$ResistanceClass,outcomeclasses)),]

# create long data
factorvardflong<-melt(factorvardf, id.vars=c('ResistanceClass','FactorVariable','FactorLevel'), measure.vars = c('logOddsRatio','logOddsRatio_unadjusted','logOddsRatio_mainmodel'), variable.name = 'model',value.name='logOddsRatio')

# append colour, shape
factorvardflong$colour<-"gray30"
factorvardflong$shape<-16
factorvardflongsplit<-split(factorvardflong,factorvardflong$FactorLevel)
factorvardflong<-do.call(rbind,lapply(factorvardflongsplit,recolourdf))


# set factor levels
factorvardflong$FactorLevel<-factor(factorvardflong$FactorLevel,levels=c('Enterobacteriaceae', 'Proteobacteria_non-Enterobacteriaceae', 'Firmicutes','other'))
factorvardflong$model<-factor(factorvardflong$model,levels=c('logOddsRatio', 'logOddsRatio_mainmodel', 'logOddsRatio_unadjusted'))
factorvardflong<-factorvardflong[order(factorvardflong$model), ]  # re-order dataframe based on model

# set parameters
params<-setparameters('HostTaxonomy')
lowerlim<-params[1];upperlim<-params[2];mywidth<-params[3];myheight<-params[4]

# plot
p<-ggplot(factorvardflong,aes(x=ResistanceClass,y=logOddsRatio)) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + facet_wrap(~ FactorLevel, as.table=FALSE, ncol = ncols, labeller = as_labeller(labelfunc)) + geom_point(shape=factorvardflong$shape,colour=factorvardflong$colour,size=2)
p<-p + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) + scale_x_discrete(limits=outcomeclasses,labels=gsub('.','-',outcomeclasses,fixed=TRUE))

pdf('output_adjusted/minus_RepliconCarriage_NumOtherResistanceClasses/parametric_comparison/HostTaxonomy_faceted_vs_mainmodel_vs_unadjusted.pdf',mywidth,myheight)
p
dev.off()


# BiocideMetalResistance: alternative model with Integron and NumOtherResistanceClasses removed
myggtitle<-'Biocide/metal resistance gene presence\nreference: absence'
adjustedalternativedf<-read.table("output_adjusted/minus_Integron_NumOtherResistanceClasses/adjustedodds_vs_mainmodel.tsv",sep='\t',header=TRUE,as.is=TRUE)
factorvardf<-adjustedalternativedf[adjustedalternativedf$FactorVariable=='BiocideMetalResistance',]
factorvardf$FactorLevel<-'presence'
factorvardf<-factorvardf[order(match(factorvardf$ResistanceClass,outcomeclasses)),]

# create long data
factorvardflong<-melt(factorvardf, id.vars=c('ResistanceClass','FactorVariable','FactorLevel'), measure.vars = c('logOddsRatio','logOddsRatio_unadjusted','logOddsRatio_mainmodel'), variable.name = 'model',value.name='logOddsRatio')

# append colour, shape
factorvardflong$colour<-"gray30"
factorvardflong$shape<-16
factorvardflongsplit<-split(factorvardflong,factorvardflong$FactorLevel)
factorvardflong<-do.call(rbind,lapply(factorvardflongsplit,recolourdf))

# set factor levels
factorvardflong$FactorLevel<-factor(factorvardflong$FactorLevel,levels=c('absence','presence'))
factorvardflong$model<-factor(factorvardflong$model,levels=c('logOddsRatio', 'logOddsRatio_mainmodel', 'logOddsRatio_unadjusted'))
factorvardflong<-factorvardflong[order(factorvardflong$model), ]  # re-order dataframe based on model

# set parameters
params<-setparameters('BiocideMetalResistance')
lowerlim<-params[1];upperlim<-params[2];mywidth<-params[3];myheight<-params[4]

# plot
p<-ggplot(factorvardflong,aes(x=ResistanceClass,y=logOddsRatio)) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + facet_wrap(~ FactorLevel, as.table=FALSE,ncol = ncols) + geom_point(shape=factorvardflong$shape,colour=factorvardflong$colour,size=2)
p<-p + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) + scale_x_discrete(limits=outcomeclasses,labels=gsub('.','-',outcomeclasses,fixed=TRUE))

pdf('output_adjusted/minus_Integron_NumOtherResistanceClasses/parametric_comparison/BiocideMetalResistance_faceted_vs_mainmodel_vs_unadjusted.pdf',mywidth,myheight)
p
dev.off()


# Integron: alternative model with BiocideMetalResistance and NumOtherResistanceClasses removed
myggtitle<-'Integron presence\nreference: absence'
adjustedalternativedf<-read.table("output_adjusted/minus_BiocideMetalResistance_NumOtherResistanceClasses/adjustedodds_vs_mainmodel.tsv",sep='\t',header=TRUE,as.is=TRUE)
factorvardf<-adjustedalternativedf[adjustedalternativedf$FactorVariable=='Integron',]
factorvardf$FactorLevel<-'presence'
factorvardf<-factorvardf[order(match(factorvardf$ResistanceClass,outcomeclasses)),]

# create long data
factorvardflong<-melt(factorvardf, id.vars=c('ResistanceClass','FactorVariable','FactorLevel'), measure.vars = c('logOddsRatio','logOddsRatio_unadjusted','logOddsRatio_mainmodel'), variable.name = 'model',value.name='logOddsRatio')

# append colour, shape
factorvardflong$colour<-"gray30"
factorvardflong$shape<-16
factorvardflongsplit<-split(factorvardflong,factorvardflong$FactorLevel)
factorvardflong<-do.call(rbind,lapply(factorvardflongsplit,recolourdf))

# set factor levels
factorvardflong$FactorLevel<-factor(factorvardflong$FactorLevel,levels=c('absence','presence'))
factorvardflong$model<-factor(factorvardflong$model,levels=c('logOddsRatio', 'logOddsRatio_mainmodel', 'logOddsRatio_unadjusted'))
factorvardflong<-factorvardflong[order(factorvardflong$model), ]  # re-order dataframe based on model

# set parameters
params<-setparameters('Integron')
lowerlim<-params[1];upperlim<-params[2];mywidth<-params[3];myheight<-params[4]

# plot
p<-ggplot(factorvardflong,aes(x=ResistanceClass,y=logOddsRatio)) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + facet_wrap(~ FactorLevel, as.table=FALSE,ncol = ncols) + geom_point(shape=factorvardflong$shape,colour=factorvardflong$colour,size=2)
p<-p + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) + scale_x_discrete(limits=outcomeclasses,labels=gsub('.','-',outcomeclasses,fixed=TRUE))

pdf('output_adjusted/minus_BiocideMetalResistance_NumOtherResistanceClasses/parametric_comparison/Integron_faceted_vs_mainmodel_vs_unadjusted.pdf',mywidth,myheight)
p
dev.off()


# ConjugativeSystem: alternative model with log10PlasmidSize removed
myggtitle<-'Conjugative system\nreference: non-mobilisable'
adjustedalternativedf<-read.table("output_adjusted/minus_log10PlasmidSize/adjustedodds_vs_mainmodel.tsv",sep='\t',header=TRUE,as.is=TRUE)
factorvardf<-adjustedalternativedf[adjustedalternativedf$FactorVariable=='ConjugativeSystem',]
factorvardf<-factorvardf[order(match(factorvardf$ResistanceClass,outcomeclasses)),]

# create long data
factorvardflong<-melt(factorvardf, id.vars=c('ResistanceClass','FactorVariable','FactorLevel'), measure.vars = c('logOddsRatio','logOddsRatio_unadjusted','logOddsRatio_mainmodel'), variable.name = 'model',value.name='logOddsRatio')

# append colour, shape
factorvardflong$colour<-"gray30"
factorvardflong$shape<-16
factorvardflongsplit<-split(factorvardflong,factorvardflong$FactorLevel)
factorvardflong<-do.call(rbind,lapply(factorvardflongsplit,recolourdf))

# set factor levels
factorvardflong$FactorLevel<-factor(factorvardflong$FactorLevel,levels=c('mobilisable','conjugative'))
factorvardflong$model<-factor(factorvardflong$model,levels=c('logOddsRatio', 'logOddsRatio_mainmodel', 'logOddsRatio_unadjusted'))
factorvardflong<-factorvardflong[order(factorvardflong$model), ]  # re-order dataframe based on model

# set parameters
params<-setparameters('ConjugativeSystem')
lowerlim<-params[1];upperlim<-params[2];mywidth<-params[3];myheight<-params[4]

# plot
p<-ggplot(factorvardflong,aes(x=ResistanceClass,y=logOddsRatio)) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + facet_wrap(~ FactorLevel, as.table=FALSE,ncol = ncols) + geom_point(shape=factorvardflong$shape,colour=factorvardflong$colour,size=2)
p<-p + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) + scale_x_discrete(limits=outcomeclasses,labels=gsub('.','-',outcomeclasses,fixed=TRUE))

pdf('output_adjusted/minus_log10PlasmidSize/parametric_comparison/ConjugativeSystem_faceted_vs_mainmodel_vs_unadjusted.pdf',mywidth,myheight)
p
dev.off()


# ConjugativeSystem: alternative model with 3 associated factors removed
myggtitle<-'Conjugative system\nreference: non-mobilisable'
adjustedalternativedf<-read.table("output_adjusted/minus_3associatedfactorsofConjugativeSystem/adjustedodds_vs_mainmodel.tsv",sep='\t',header=TRUE,as.is=TRUE)
factorvardf<-adjustedalternativedf[adjustedalternativedf$FactorVariable=='ConjugativeSystem',]
factorvardf<-factorvardf[order(match(factorvardf$ResistanceClass,outcomeclasses)),]

# create long data
factorvardflong<-melt(factorvardf, id.vars=c('ResistanceClass','FactorVariable','FactorLevel'), measure.vars = c('logOddsRatio','logOddsRatio_unadjusted','logOddsRatio_mainmodel'), variable.name = 'model',value.name='logOddsRatio')

# append colour, shape
factorvardflong$colour<-"gray30"
factorvardflong$shape<-16
factorvardflongsplit<-split(factorvardflong,factorvardflong$FactorLevel)
factorvardflong<-do.call(rbind,lapply(factorvardflongsplit,recolourdf))

# set factor levels
factorvardflong$FactorLevel<-factor(factorvardflong$FactorLevel,levels=c('mobilisable','conjugative'))
factorvardflong$model<-factor(factorvardflong$model,levels=c('logOddsRatio', 'logOddsRatio_mainmodel', 'logOddsRatio_unadjusted'))
factorvardflong<-factorvardflong[order(factorvardflong$model), ]  # re-order dataframe based on model

# set parameters
params<-setparameters('ConjugativeSystem')
lowerlim<-params[1];upperlim<-params[2];mywidth<-params[3];myheight<-params[4]

# plot
p<-ggplot(factorvardflong,aes(x=ResistanceClass,y=logOddsRatio)) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + facet_wrap(~ FactorLevel, as.table=FALSE,ncol = ncols) + geom_point(shape=factorvardflong$shape,colour=factorvardflong$colour,size=2)
p<-p + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) + scale_x_discrete(limits=outcomeclasses,labels=gsub('.','-',outcomeclasses,fixed=TRUE))

pdf('output_adjusted/minus_3associatedfactorsofConjugativeSystem/parametric_comparison/ConjugativeSystem_faceted_vs_mainmodel_vs_unadjusted.pdf',mywidth,myheight)
p
dev.off()


# ConjugativeSystem: alternative model with 6 associated factors removed
myggtitle<-'Conjugative system\nreference: non-mobilisable'
adjustedalternativedf<-read.table("output_adjusted/minus_6associatedfactorsofConjugativeSystem/adjustedodds_vs_mainmodel.tsv",sep='\t',header=TRUE,as.is=TRUE)
factorvardf<-adjustedalternativedf[adjustedalternativedf$FactorVariable=='ConjugativeSystem',]
factorvardf<-factorvardf[order(match(factorvardf$ResistanceClass,outcomeclasses)),]

# create long data
factorvardflong<-melt(factorvardf, id.vars=c('ResistanceClass','FactorVariable','FactorLevel'), measure.vars = c('logOddsRatio','logOddsRatio_unadjusted','logOddsRatio_mainmodel'), variable.name = 'model',value.name='logOddsRatio')

# append colour, shape
factorvardflong$colour<-"gray30"
factorvardflong$shape<-16
factorvardflongsplit<-split(factorvardflong,factorvardflong$FactorLevel)
factorvardflong<-do.call(rbind,lapply(factorvardflongsplit,recolourdf))

# set factor levels
factorvardflong$FactorLevel<-factor(factorvardflong$FactorLevel,levels=c('mobilisable','conjugative'))
factorvardflong$model<-factor(factorvardflong$model,levels=c('logOddsRatio', 'logOddsRatio_mainmodel', 'logOddsRatio_unadjusted'))
factorvardflong<-factorvardflong[order(factorvardflong$model), ]  # re-order dataframe based on model

# set parameters
params<-setparameters('ConjugativeSystem')
lowerlim<-params[1];upperlim<-params[2];mywidth<-params[3];myheight<-params[4]

# plot
p<-ggplot(factorvardflong,aes(x=ResistanceClass,y=logOddsRatio)) + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + facet_wrap(~ FactorLevel, as.table=FALSE,ncol = ncols) + geom_point(shape=factorvardflong$shape,colour=factorvardflong$colour,size=2)
p<-p + theme(axis.text.x = element_text(angle=45,vjust=1,hjust=1,size=rel(1.3)),axis.text.y = element_text(size=rel(1.3)), strip.text.x = element_text(size=rel(1.3)), panel.grid.major=element_blank(), panel.grid.minor=element_blank(),axis.title.y = element_text(size=rel(1.3)), axis.title.x = element_blank(), plot.title = element_text(size=12)) + ylab('\nlog odds ratio (95% CI)') + ggtitle(myggtitle)
p<-p + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) + scale_x_discrete(limits=outcomeclasses,labels=gsub('.','-',outcomeclasses,fixed=TRUE))

pdf('output_adjusted/minus_6associatedfactorsofConjugativeSystem/parametric_comparison/ConjugativeSystem_faceted_vs_mainmodel_vs_unadjusted.pdf',mywidth,myheight)
p
dev.off()


