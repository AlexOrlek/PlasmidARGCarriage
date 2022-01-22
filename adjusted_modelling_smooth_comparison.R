pacman::p_load(tidyverse, mgcv, mgcv.helper, RColorBrewer, mgcViz, grid, gsubfn)

alternativemodelnames<-c('mainmodel_minus_NumOtherResistanceClasses','mainmodel_minus_Integron','mainmodel_minus_BiocideMetalResistance','mainmodel_minus_RepliconCarriage','mainmodel_minus_HostTaxonomy','mainmodel_minus_log10PlasmidSize','mainmodel_minus_RepliconCarriage_NumOtherResistanceClasses','mainmodel_minus_BiocideMetalResistance_NumOtherResistanceClasses','mainmodel_minus_Integron_NumOtherResistanceClasses','mainmodel_minus_BiocideMetalResistance_Integron')

outcomeclasses<-c('aminoglycoside','phenicol','sulphonamide','tetracycline','macrolide','TEM.1','trimethoprim','ESBL', 'carbapenem','quinolone','colistin')
outputnames<-c('log10PlasmidSize','InsertionSequenceDensity','NumOtherResistanceClasses','CollectionDate')
xlabs<-c('log10 Plasmid size (centred on 10 kb)','Insertion sequence density (frequency per 10 kb)','Other resistance gene classes','Years since reference collection year')
ggtitles<-c('        log10 Plasmid size (kb)\n        reference: 10 kb\n','        Insertion sequence density\n        reference: 0\n','        Number of other resistance gene classes\n        reference: 0\n','        Collection date\n        reference: collection year 1994\n')


brewerpal1<-brewer.pal(8,'Set1')
brewerpal<-c(brewerpal1,'#e0bb6d','#90EE90','#add8e6')
brewerpal[6]<-'#ffd700'

# ---------------------------
# N.B loading probability scale model data because this stored both logodds scale (y) and probability scale (ty) coordinate data

# first get smooth coordinate data from mainmodel
mainmodelsmoothplotlist_prob<-readRDS('output_adjusted/mainmodel/coefficientplots/probscale/smoothplotlist2.rds')
mainmodeldatalist_prob<-list()
for (i in 1:length(outputnames)) {
  outputname<-outputnames[i]
  mainmodeldatalist_prob[[outputname]]<-list()
  for (j in 1:length(outcomeclasses)) {
    outcomeclass<-outcomeclasses[j]
    mainmodeldatalist_prob[[outputname]][[outcomeclass]]<-mainmodelsmoothplotlist_prob[[outputname]][[outcomeclass]]$data$fit
    mainmodeldatalist_prob[[outputname]][[outcomeclass]]$modeltype<-'adjusted'
    mainmodeldatalist_prob[[outputname]][[outcomeclass]]$colour<-brewerpal[j]
  }
}

# get smooth coordinate data from unadjusted models
unadjustedmodeldatalist_prob<-list()
for (i in 1:length(outputnames)) {
  outputname<-outputnames[i]
  unadjustedmodelsmoothplotlist_prob<-readRDS(gsub('%s',outputname,'output_unadjusted/%s/coefficientplots/probscale/smoothplotlist2.rds'))
  unadjustedmodeldatalist_prob[[outputname]]<-list()
  for (j in 1:length(outcomeclasses)) {
    outcomeclass<-outcomeclasses[j]
    unadjustedmodeldatalist_prob[[outputname]][[outcomeclass]]<-unadjustedmodelsmoothplotlist_prob[[outputname]][[outcomeclass]]$data$fit
    unadjustedmodeldatalist_prob[[outputname]][[outcomeclass]]$modeltype<-'unadjusted'
    unadjustedmodeldatalist_prob[[outputname]][[outcomeclass]]$colour<-'black'
  }
}


# get smooth coordinate data from alternative models and append data from main model
# store data in nested list: alternativemodelname -> outputname -> outcomeclass
datalist_prob<-list()
for (modelname in alternativemodelnames) {
  print(modelname)
  smoothplotlist_prob<-readRDS(gsub('%s',modelname,'output_adjusted/%s/coefficientplots/probscale/smoothplotlist2.rds'))
  datalist_prob[[modelname]]<-list()
  for (i in 1:length(outputnames)) {
    outputname<-outputnames[i]
    #skip where outputname not applicable to model
    if (modelname %in% c('mainmodel_minus_NumOtherResistanceClasses','mainmodel_minus_RepliconCarriage_NumOtherResistanceClasses','mainmodel_minus_BiocideMetalResistance_NumOtherResistanceClasses','mainmodel_minus_Integron_NumOtherResistanceClasses') && outputname=='NumOtherResistanceClasses') {
      next
    }
    if (modelname %in% c('mainmodel_minus_log10PlasmidSize','mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity') && outputname=='log10PlasmidSize') {
      next
    }
    if (modelname %in% c('mainmodel_minus_InsertionSequenceDensity','mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity') && outputname=='InsertionSequenceDensity') {
      next
    }
    datalist_prob[[modelname]][[outputname]]<-list()
    for (j in 1:length(outcomeclasses)) {
      outcomeclass<-outcomeclasses[j]
      datalist_prob[[modelname]][[outputname]][[outcomeclass]]<-smoothplotlist_prob[[outputname]][[outcomeclass]]$data$fit
      datalist_prob[[modelname]][[outputname]][[outcomeclass]]$modeltype<-'alternative_adjusted'
      datalist_prob[[modelname]][[outputname]][[outcomeclass]]$colour<-'black'
      datalist_prob[[modelname]][[outputname]][[outcomeclass]]<-rbind(datalist_prob[[modelname]][[outputname]][[outcomeclass]],mainmodeldatalist_prob[[outputname]][[outcomeclass]],unadjustedmodeldatalist_prob[[outputname]][[outcomeclass]])
    }
  }
}


# ---------------------------
# create plots
smoothplotcomparisonlist_logodds<-list()
smoothplotcomparisonlist_prob<-list()
for (modelname in alternativemodelnames) {
  smoothplotcomparisonlist_logodds[[modelname]]<-list()
  smoothplotcomparisonlist_prob[[modelname]]<-list()
  for (i in 1:length(outputnames)) {
    outputname<-outputnames[i]
    #skip where outputname not applicable to model
    if (modelname %in% c('mainmodel_minus_NumOtherResistanceClasses','mainmodel_minus_RepliconCarriage_NumOtherResistanceClasses','mainmodel_minus_BiocideMetalResistance_NumOtherResistanceClasses','mainmodel_minus_Integron_NumOtherResistanceClasses') && outputname=='NumOtherResistanceClasses') {
      next
    }
    if (modelname %in% c('mainmodel_minus_log10PlasmidSize','mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity') && outputname=='log10PlasmidSize') {
      next
    }
    if (modelname %in% c('mainmodel_minus_InsertionSequenceDensity','mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity') && outputname=='InsertionSequenceDensity') {
      next
    }
    #logodds limits
    lowerlim_logodds<--5
    upperlim_logodds<-5
    if (outputname=='InsertionSequenceDensity') {
      lowerlim_logodds<--4
      upperlim_logodds<-6
    }
    if (outputname=='NumOtherResistanceClasses') {
      lowerlim_logodds<--4
      upperlim_logodds<-8
    }
    if (outputname=='CollectionDate') {
      lowerlim_logodds<--2
      upperlim_logodds<-10
    }
    #prob limits
    lowerlim_prob<-0
    upperlim_prob<-0.2
    if (outputname=='NumOtherResistanceClasses') {
      upperlim_prob<-1
    }
    smoothplotcomparisonlist_logodds[[modelname]][[outputname]]<-list()
    smoothplotcomparisonlist_prob[[modelname]][[outputname]]<-list()
    for (j in 1:length(outcomeclasses)) {
      outcomeclass<-outcomeclasses[j]
      #log-odds plot
      p_logodds<-ggplot(datalist_prob[[modelname]][[outputname]][[outcomeclass]],aes(x=x,y=y,colour=modeltype,linetype=modeltype)) +
        geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + geom_line() + theme_bw() + xlab('') + ylab('') + ggtitle(gsub('.','-',outcomeclass,fixed=TRUE)) +
        guides(colour=FALSE,linetype=FALSE) + 
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim_logodds,upperlim_logodds)) +
        scale_x_continuous(breaks=scales::pretty_breaks()) +
        scale_colour_manual(values=datalist_prob[[modelname]][[outputname]][[outcomeclass]] %>% group_by(modeltype) %>% distinct(colour) %>% arrange(match(modeltype, c('adjusted','alternative_adjusted','unadjusted'))) %>% ungroup() %>% select(colour) %>% unlist() %>% unname()) +
        scale_linetype_manual(values=c('solid','dashed', 'dotted'))
      #remove redundant plot elements
      if (j<=5) {
        p_logodds<-p_logodds+theme(axis.text.x = element_text(colour='white'),axis.ticks.x = element_blank())
      }
      if (!j %in% c(1,7)) {
        p_logodds<-p_logodds+theme(axis.text.y = element_text(colour='white'),axis.ticks.y = element_blank())
      }
      #further customise plots
      if (j == 1) {
        p_logodds <- p_logodds + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-12,-4),'pt'),axis.text=element_text(size=rel(1.1)))
      } else if (j == 7) {
        p_logodds <- p_logodds + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-4,-4),'pt'),axis.text=element_text(size=rel(1.1)))
      } else if (j == 6) {
        p_logodds <- p_logodds + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,12,-12,-14),'pt'),axis.text=element_text(size=rel(1.1)))
      } else if (j > 7) {
        p_logodds <- p_logodds + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-4,-14),'pt'),axis.text=element_text(size=rel(1.1)))
      } else {
        p_logodds <- p_logodds + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-12,-14),'pt'),axis.text=element_text(size=rel(1.1)))
      }
      #prob plot (N.B no geom_hline since different models linked to different baseline probabilities)
      p_prob<-ggplot(datalist_prob[[modelname]][[outputname]][[outcomeclass]],aes(x=x,y=ty,colour=modeltype,linetype=modeltype)) +
        geom_line() + theme_bw() + xlab('') + ylab('') + ggtitle(gsub('.','-',outcomeclass,fixed=TRUE)) +
        guides(colour=FALSE,linetype=FALSE) + 
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim_prob,upperlim_prob)) +
        scale_x_continuous(breaks=scales::pretty_breaks()) +
        scale_colour_manual(values=datalist_prob[[modelname]][[outputname]][[outcomeclass]] %>% group_by(modeltype) %>% distinct(colour) %>% arrange(match(modeltype, c('adjusted','alternative_adjusted','unadjusted'))) %>% ungroup() %>% select(colour) %>% unlist() %>% unname()) +
        scale_linetype_manual(values=c('solid','dashed', 'dotted'))
      #remove redundant plot elements
      if (j<=5) {
        p_prob<-p_prob+theme(axis.text.x = element_text(colour='white'),axis.ticks.x = element_blank())
      }
      if (!j %in% c(1,6)) {
        p_prob<-p_prob+theme(axis.text.y = element_text(colour='white'),axis.ticks.y = element_blank())
      }
      #further customise plot
      if (j == 1) {
        p_prob <- p_prob + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-12,-4),'pt'),axis.text=element_text(size=rel(1.1)))
      } else if (j == 7) {
        p_prob <- p_prob + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-4,-4),'pt'),axis.text=element_text(size=rel(1.1)))
      } else if (j == 6) {
        p_prob <- p_prob + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,12,-12,-14),'pt'),axis.text=element_text(size=rel(1.1)))
      } else if (j > 7) {
        p_prob <- p_prob + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-4,-14),'pt'),axis.text=element_text(size=rel(1.1)))
      } else {
        p_prob <- p_prob + theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-12,-14),'pt'),axis.text=element_text(size=rel(1.1)))
      }
      #assign plots to list
      smoothplotcomparisonlist_logodds[[modelname]][[outputname]][[outcomeclass]]<-p_logodds
      smoothplotcomparisonlist_prob[[modelname]][[outputname]][[outcomeclass]]<-p_prob
    } 
  }
}






# ---------------------------
# save plots as gridplots
width=12
#height=4.8  #5.2
height=5.6
#width=10
#height=3.75
for (modelname in alternativemodelnames) {
  for (i in 1:length(outputnames)) {
    outputname<-outputnames[i]
    #skip where outputname not applicable to model
    if (modelname %in% c('mainmodel_minus_NumOtherResistanceClasses','mainmodel_minus_RepliconCarriage_NumOtherResistanceClasses','mainmodel_minus_BiocideMetalResistance_NumOtherResistanceClasses','mainmodel_minus_Integron_NumOtherResistanceClasses') && outputname=='NumOtherResistanceClasses') {
      next
    }
    if (modelname %in% c('mainmodel_minus_log10PlasmidSize','mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity') && outputname=='log10PlasmidSize') {
      next
    }
    if (modelname %in% c('mainmodel_minus_InsertionSequenceDensity','mainmodel_minus_log10PlasmidSize_InsertionSequenceDensity') && outputname=='InsertionSequenceDensity') {
      next
    }
    #log-odds
    pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outputname),'output_adjusted/%1/coefficientplots/logoddsscale/%2_vs_mainmodel.pdf'),width=width,height=height)
    gridPrint(grobs=smoothplotcomparisonlist_logodds[[modelname]][[outputname]][outcomeclasses],nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15), vjust=-0.1),left=textGrob('Effect on log odds',gp=gpar(fontsize=15),rot=90, vjust=1),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
    dev.off()
    #prob
    pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outputname),'output_adjusted/%1/coefficientplots/probscale/%2_vs_mainmodel.pdf'),width=width,height=height)
    gridPrint(grobs=smoothplotcomparisonlist_prob[[modelname]][[outputname]][outcomeclasses],nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15), vjust=-0.1),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90, vjust=1),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
    dev.off()
    #prob y0to1
    pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outputname),'output_adjusted/%1/coefficientplots/probscale/ylim0to1_%2_vs_mainmodel.pdf'),width=width,height=height)
    rescaledplots<-lapply(smoothplotcomparisonlist_prob[[modelname]][[outputname]], function(x) x + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,1)))
    gridPrint(grobs=rescaledplots[outcomeclasses],nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15), vjust=-0.1),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90, vjust=1),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
    dev.off()
  }
}

