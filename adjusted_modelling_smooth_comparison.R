setwd("C:/Users/alexo/Dropbox/Oxford_2015/thesis/CORRECTIONS/Chapter_5/resistance_prediction")
library(mgcv)
library(RColorBrewer)
library(mgcViz)
library(mgcv.helper)
library(gsubfn)
library(grid)


alternativemodelnames<-c('mainmodel_minusnumotherresgenes','mainmodel_minusintegron','mainmodel_minusbacmet','mainmodel_minusreptypes','mainmodel_minustaxa',
                         'mainmodel_minusreptypes_numotherresgenes','mainmodel_minusbacmet_numotherresgenes','mainmodel_minusintegron_numotherresgenes','mainmodel_minusbacmet_integron')
outcomeclasses<-c('aminoglycoside','betalactam_carbapenem','betalactam_ESBL','betalactam_other','macrolide','phenicol','quinolone','sulphonamide','tetracycline','trimethoprim')
outputnames<-c('logplasmidsize','isscore','numotherresgenes','coldateimputed')
xlabs<-c('log10 Plasmid size (centred on 10 kb)','Insertion sequence density (frequency per 10 kb)','Other resistance gene classes','Years since initial collection year')
ggtitles<-c('      log10 Plasmid size (kb)\n      baseline: 10 kb\n','      Insertion sequence density\n      baseline: 0\n','      Number of other resistance gene classes\n      baseline: 0\n','      Collection date\n      baseline: initial year (1994)\n')


brewerpal1<-brewer.pal(8,'Set1')
brewerpal1[6]<-'#CDCD00'
betalactampal<-brewer.pal(9,'Blues')[c(3,6,9)]
brewerpal2<-c(brewerpal1[1],betalactampal,brewerpal1[3:length(brewerpal1)])
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


###first get smooth coordinate data from mainmodel
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


###get smooth coordinate data from alternative models and append data from main model
for (modelname in alternativemodelnames) {
  print(modelname)
  smoothplotlist_prob<-readRDS(gsub('%s',modelname,'output_adjusted/%s/coefficientplots/probscale/smoothplotlist2.rds'))
  datalist_prob<-list()
  for (i in 1:length(outputnames)) {
    outputname<-outputnames[i]
    if (modelname %in% c('mainmodel_minusnumotherresgenes','mainmodel_minusreptypes_numotherresgenes','mainmodel_minusbacmet_numotherresgenes','mainmodel_minusintegron_numotherresgenes') && outputname=='numotherresgenes') {
      next
    }
    datalist_prob[[outputname]]<-list()
    for (j in 1:length(outcomeclasses)) {
      outcomeclass<-outcomeclasses[j]
      datalist_prob[[outputname]][[outcomeclass]]<-smoothplotlist_prob[[outputname]][[outcomeclass]]$data$fit
      datalist_prob[[outputname]][[outcomeclass]]$modeltype<-'alternative_adjusted'
      datalist_prob[[outputname]][[outcomeclass]]$colour<-'black'
      datalist_prob[[outputname]][[outcomeclass]]<-rbind(datalist_prob[[outputname]][[outcomeclass]],mainmodeldatalist_prob[[outputname]][[outcomeclass]])
    }
  }
}

###plot
for (modelname in alternativemodelnames) {
  smoothplotcomparisonlist_logodds<-list()
  smoothplotcomparisonlist_prob<-list()
  for (i in 1:length(outputnames)) {
    outputname<-outputnames[i]
    print(outputname)
    if (modelname %in% c('mainmodel_minusnumotherresgenes','mainmodel_minusreptypes_numotherresgenes','mainmodel_minusbacmet_numotherresgenes','mainmodel_minusintegron_numotherresgenes') && outputname=='numotherresgenes') {
      next
    }
    smoothplotcomparisonlistnest_logodds<-list()
    smoothplotcomparisonlistnest_prob<-list()
    outputname<-outputnames[i]
    myggtitle<-ggtitles[i]
    myxlab<-xlabs[i]
    #logodds limits
    lowerlim_logodds<--5
    upperlim_logodds<-5
    if (outputname=='numotherresgenes') {
      lowerlim_logodds<--3
      upperlim_logodds<-7
    }
    if (outputname=='coldateimputed') {
      lowerlim_logodds<--8
      upperlim_logodds<-2
    }
    #prob limits
    lowerlim_prob<-0
    upperlim_prob<-0.2
    if (outputname=='numotherresgenes') {
      upperlim_prob<-1
    }
    for (j in 1:length(outcomeclasses)) {
      outcomeclass<-outcomeclasses[j]
      ###log-odds plot
      p_logodds<-ggplot(datalist_prob[[outputname]][[outcomeclass]],aes(x=x,y=y,colour=modeltype,linetype=modeltype)) +
        geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + geom_line() + theme_bw() +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,2,-12,-12),'pt'),axis.text=element_text(size=rel(1.1))) +
        xlab('') + ylab('') + ggtitle(outcomeclass) +
        guides(colour=FALSE,linetype=FALSE) + 
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim_logodds,upperlim_logodds)) +
        scale_x_continuous(breaks=scales::pretty_breaks()) +
        scale_colour_manual(values=rev(unique(datalist_prob[[outputname]][[outcomeclass]]$colour))) +
        scale_linetype_manual(values=c('solid','dotted'))
      #remove redundant plot elements
      if (j<=5) {
        p_logodds<-p_logodds+theme(axis.text.x = element_text(colour='white'),axis.ticks.x = element_blank())
      }
      if (!j %in% c(1,6)) {
        p_logodds<-p_logodds+theme(axis.text.y = element_text(colour='white'),axis.ticks.y = element_blank())
      }
      ###prob plot (N.B no geom_hline since different models linked to different baseline probabilities)
      p_prob<-ggplot(datalist_prob[[outputname]][[outcomeclass]],aes(x=x,y=ty,colour=modeltype,linetype=modeltype)) +
        geom_line() + theme_bw() +
        theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,2,-12,-12),'pt'),axis.text=element_text(size=rel(1.1))) +
        xlab('') + ylab('') + ggtitle(outcomeclass) +
        guides(colour=FALSE,linetype=FALSE) + 
        scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim_prob,upperlim_prob)) +
        scale_x_continuous(breaks=scales::pretty_breaks()) +
        scale_colour_manual(values=rev(unique(datalist_prob[[outputname]][[outcomeclass]]$colour))) +
        scale_linetype_manual(values=c('solid','dotted'))
      #remove redundant plot elements
      if (j<=5) {
        p_prob<-p_prob+theme(axis.text.x = element_text(colour='white'),axis.ticks.x = element_blank())
      }
      if (!j %in% c(1,6)) {
        p_prob<-p_prob+theme(axis.text.y = element_text(colour='white'),axis.ticks.y = element_blank())
      }
      
      #assign plots to list
      smoothplotcomparisonlistnest_logodds[[outcomeclass]]<-p_logodds
      smoothplotcomparisonlistnest_prob[[outcomeclass]]<-p_prob
    }
    smoothplotcomparisonlist_logodds[[outputname]]<-smoothplotcomparisonlistnest_logodds
    smoothplotcomparisonlist_prob[[outputname]]<-smoothplotcomparisonlistnest_prob
  }
  
  width=12
  #height=4.8  #5.2
  height=5.6
  #width=10
  #height=3.75
  for (i in 1:length(outputnames)) {
    outputname<-outputnames[i]
    if (modelname %in% c('mainmodel_minusnumotherresgenes','mainmodel_minusreptypes_numotherresgenes','mainmodel_minusbacmet_numotherresgenes','mainmodel_minusintegron_numotherresgenes') && outputname=='numotherresgenes') {
      next
    }
    ###log-odds
    pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outputname),'output_adjusted/%1/coefficientplots/logoddsscale/%2_vs_mainmodel.pdf'),width=width,height=height)
    gridPrint(grobs=smoothplotcomparisonlist_logodds[[outputname]][order(outcomeclasses)],nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15)),left=textGrob('Effect on log odds',gp=gpar(fontsize=15),rot=90),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
    dev.off()
    ###prob
    pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outputname),'output_adjusted/%1/coefficientplots/probscale/%2_vs_mainmodel.pdf'),width=width,height=height)
    gridPrint(grobs=smoothplotcomparisonlist_prob[[outputname]][order(outcomeclasses)],nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15)),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
    dev.off()
    ###prob y0to1
    pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outputname),'output_adjusted/%1/coefficientplots/probscale/ylim0to1_%2_vs_mainmodel.pdf'),width=width,height=height)
    rescaledplots<-lapply(smoothplotcomparisonlist_prob[[outputname]], function(x) x + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,1)))
    gridPrint(grobs=rescaledplots[order(outcomeclasses)],nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15)),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90),top = grid::textGrob(ggtitles[i], x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=13,lineheight=1)))
    dev.off()
  }
}



###TESTING - check bug with minusreptypes_numotherresgenes


###TESTING

#test<-lapply(smoothplotcomparisonlist_prob[[outputname]], function(x) x + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,1)))
#gridPrint(grobs=test[order(outcomeclasses)],nrow=2,bottom=textGrob(xlabs[i],gp=gpar(fontsize=15)),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90))


###TESTING

# #problem: data stored in probability scale plots is the same...transformation is applied to the data when I plot. Use logodds data and re-apply transformation
# #N.B !!!ty column is transformed y - use this for probability plots.
# 
# modelname<-"mainmodel_minusintegron"
# smoothplotlist_prob<-readRDS(gsub('%s',modelname,'output_adjusted/%s/coefficientplots/probscale/smoothplotlist2.rds'))
# 
# smoothplotlist_prob$isscore$aminoglycoside$data$fit  #!data is relative to intercept
# 
# head(smoothplotlist_prob$isscore$betalactam_carbapenem$data$fit)
# 
# 
# plot(smoothplotlist_prob$isscore$betalactam_carbapenem$data$fit$x,smoothplotlist_prob$isscore$betalactam_carbapenem$data$fit$ty)
# 
# smoothplotlist_logodds<-readRDS(gsub('%s',modelname,'output_adjusted/%s/coefficientplots/logoddsscale/smoothplotlist1.rds'))
# 
# head(smoothplotlist_logodds$isscore$betalactam_carbapenem$data$fit)
# 
# plot(smoothplotlist_logodds$isscore$betalactam_carbapenem$data$fit$x,smoothplotlist_logodds$isscore$betalactam_carbapenem$data$fit$y+smoothplotlist_logodds$isscore$betalactam_carbapenem$data$fit$ty)
# plot(smoothplotlist_logodds$isscore$betalactam_carbapenem$data$fit$x,smoothplotlist_logodds$isscore$betalactam_carbapenem$data$fit$y)
# datalist_prob<-list()
# for (i in 1:length(outputnames)) {
#   outputname<-outputnames[i]
#   if (modelname=='mainmodel_minusnumotherresgenes' && outputname=='numotherresgenes') {
#     next
#   }
#   datalist_prob[[outputname]]<-list()
#   for (j in 1:length(outcomeclasses)) {
#     outcomeclass<-outcomeclasses[j]
#     datalist_prob[[outputname]][[outcomeclass]]<-smoothplotlist_prob[[outputname]][[outcomeclass]]$data$fit
#     datalist_prob[[outputname]][[outcomeclass]]$colour<-brewerpal[j]
#     datalist_prob[[outputname]][[outcomeclass]]<-rbind(datalist_prob[[outputname]][[outcomeclass]],mainmodeldatalist_prob[[outputname]][[outcomeclass]])
#   }
# }
# 
# datalist_prob$isscore$aminoglycoside
# 
# ##plot
# smoothplotcomparisonlist_prob<-list()
# i<-1
# outputname<-outputnames[i]
# if (modelname=='mainmodel_minusnumotherresgenes' && outputname=='numotherresgenes') {
#   next
# }
# smoothplotcomparisonlistnest_prob<-list()
# outputname<-outputnames[i]
# myggtitle<-ggtitles[i]
# myxlab<-xlabs[i]
# lowerlim<--5
# upperlim<-5
# if (outputname=='numotherresgenes') {
#   lowerlim<--3
#   upperlim<-7
# }
# if (outputname=='coldateimputed') {
#   lowerlim<--8
#   upperlim<-2
# }
# j<-1
# outcomeclass<-outcomeclasses[j]
# 
# ###prob plot
# p<-ggplot(datalist_prob[[outputname]][[outcomeclass]],aes(x=x,y=y,colour=colour,linetype=colour)) +
#   geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + geom_line() + theme_bw() +
#   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,2,-12,-12),'pt'),axis.text=element_text(size=rel(1.1))) + 
#   xlab('') + ylab('') + ggtitle(outcomeclass) +
#   guides(colour=FALSE,linetype=FALSE) + 
#   scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) +
#   scale_x_continuous(breaks=scales::pretty_breaks()) +
#   scale_colour_manual(values=unique(datalist_prob[[outputname]][[outcomeclass]]$colour)) +
#   scale_linetype_manual(values=c('solid','dotted'))
# p
#   
# datalist_prob[[outputname]][[outcomeclass]]



###TESTING

# datalist$isscore
# 
# smoothplotlist$isscore$aminoglycoside$data$fit
# 
# out<-ggplot_build(smoothplotlist$isscore$aminoglycoside$ggObj)
# out$data[[1]]
# out$data[[2]]
# out$data[[1]]$y>out$data[[2]]$y
# 
# 
# ggplot(out$data[[2]],aes(x=x,y=y,colour=colour)) + geom_ribbon(aes(ymin=out$data, ymax=),alpha=0.4) + geom_line(alpha=0.4) 
# 
# 
# test$isscore$aminoglycoside$data$fit
# test$isscore$aminoglycoside$ggObj$data
# 
# out<-test$isscore$aminoglycoside$data$fit
# 
# 
# 
# i<-1
# lowerlim=-5
# upperlim=5
# xlabs<-c('isscore')
# outcomeclass<-'aminoglycoside'
# ggplot(out,aes(x=x,y=y,colour='red')) + geom_line(alpha=0.4) + theme_bw() + 
#   theme(panel.grid.major = element_blank(),panel.grid.minor = element_blank(),plot.title=element_text(vjust=-7,hjust=0,colour="#525252",size=12)) + guides(colour=FALSE) + geom_hline(yintercept = 0,linetype='dashed',colour='black',size=0.3) + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) + xlab(xlabs[i]) + ylab(gsub('%s',xlabs[i],'s(%s)')) + ggtitle(paste(' ',outcomeclass))
# 
# #replace black with grey shade
# 
# test<-readRDS(gsub('%s',modelname,'output_adjusted/%s/modellist.rds'))
# 
# test$aminoglycoside
# confint.gam(test$aminoglycoside)  #only parametric terms
# 
# o<-getViz(test$aminoglycoside) #this produces gamViz object for a given class, from which smooth can be extracted and plotted using sm(o,n)
# gamcoef<-coef(test$aminoglycoside)[1]
# o1 <- plot.mgcv.smooth.1D( sm(o, i) ) #for resistance class j, extract smooth for smooth outputname i 
# out<-ggplot_build(o1$ggObj)
# out$plot$data  #can't find confidence interval data
# 
# 
# gam.vcomp(test$aminoglycoside)


##OLD
#colour="#525252"  #outcomeclass label colour
#margin=theme(plot.margin = unit(c(-10,2,-10,2), "pt"))
#gridPrint(grobs=lapply(smoothplotcomparisonlist[[outputname]][order(outcomeclasses)],'+',margin),nrow=2)


