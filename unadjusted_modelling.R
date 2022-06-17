pacman::p_load(mgcv, mgcv.helper, RColorBrewer, mgcViz, scales, grid, cowplot, devtools, gsubfn, stringr)
set.seed(42)
source('functions.R')

finaldftrunc<-read.table('data/plasmiddf_transformed.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")

# ---------------------------
# set parameters
outcomeclasses<-c('aminoglycoside','sulphonamide','tetracycline','phenicol','macrolide','trimethoprim','ESBL', 'carbapenem','quinolone','colistin')

args = commandArgs(trailingOnly=TRUE)
modelname<-args[1]
#modelname<-'log10PlasmidSize'
#modelname<-'InsertionSequenceDensity'
#modelname<-'NumOtherResistanceClasses'
#modelname<-'CollectionDate'

# create output directories
dir.create(file.path('output_unadjusted',modelname), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_unadjusted/%s'), 'coefficientplots'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_unadjusted/%s'), 'pvalues'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_unadjusted/%s/coefficientplots'), 'logoddsscale'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_unadjusted/%s/coefficientplots'), 'probscale'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_unadjusted/%s/coefficientplots'), 'plotgam'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_unadjusted/%s/coefficientplots/plotgam'), 'logoddsscale'), showWarnings = FALSE)
dir.create(file.path(gsub('%s',modelname,'output_unadjusted/%s/coefficientplots/plotgam'), 'probscale'), showWarnings = FALSE)


if (modelname=='log10PlasmidSize') {
  frmtext<-'outcome%s~s(log10PlasmidSize,k=5,pc=0)'
  myggtitle<-'       reference: 10 kb\n'
  myxlab<-'Plasmid size (kb)'
} else if (modelname=='InsertionSequenceDensity') {
  frmtext<-'outcome%s~s(InsertionSequenceDensity,k=5,pc=0)'
  myggtitle<-'       reference: 0 insertion sequences\n'
  myxlab<-'Insertion sequence density (frequency per 10 kb)'
} else if (modelname=='NumOtherResistanceClasses') {
  frmtext<-'outcome%s~s(NumOtherResistanceClasses,k=5,pc=0)'
  myggtitle<-'       reference: 0 other ARG types\n'
  myxlab<-'Number of other ARG types'
} else if (modelname=='CollectionDate') {
  frmtext<-'outcome%s~s(CollectionDate,k=5,pc=0)'
  myggtitle<-'       reference: collection year 1994\n'
  myxlab<-'Collection date'
} else {
  stop('modelname not recognised')
}


brewerpal1<-brewer.pal(8,'Set1')
brewerpal<-c(brewerpal1,'#e0bb6d','#90EE90','#add8e6')
brewerpal[6]<-'#ffd700'
#plot(1:11,1:11,col=brewerpal,pch=16,cex=6)

if (length(outcomeclasses)>11) {
  pacman::p_load(rwantshue)
  scheme <- iwanthue(seed = 42, force_init = TRUE)
  brewerpal <- scheme$hex(length(resclasses))
}
brewerpal<-brewerpal[1:length(outcomeclasses)]

# ---------------------------
# convert categorical variables to factors
finaldftrunc$BiocideMetalResistance<-as.factor(finaldftrunc$BiocideMetalResistance)
finaldftrunc$Virulence<-as.factor(finaldftrunc$Virulence)
finaldftrunc$Integron<-as.factor(finaldftrunc$Integron)
finaldftrunc$ConjugativeSystem<-factor(finaldftrunc$ConjugativeSystem, ordered = FALSE,levels = c("non-mobilisable", "mobilisable", "conjugative"))
finaldftrunc$GeographicLocation<-factor(finaldftrunc$GeographicLocation, ordered = FALSE,levels = c("high-income", "China", "United States", "other", "middle-income", "EU & UK"))
finaldftrunc$IsolationSource<-factor(finaldftrunc$IsolationSource, ordered = FALSE,levels = c("human", "livestock","other"))
finaldftrunc$RepliconCarriage<-factor(finaldftrunc$RepliconCarriage, ordered = FALSE,levels = c("untyped", "single-replicon", "multi-replicon"))
finaldftrunc$HostTaxonomy<-factor(finaldftrunc$HostTaxonomy, ordered = FALSE,levels = c("Enterobacteriaceae", "Proteobacteria_non-Enterobacteriaceae", "Firmicutes","other"))
for (outcomeclass in outcomeclasses) {
  print(outcomeclass)
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

saveRDS(modellist,file = gsub('%s',modelname,'output_unadjusted/%s/modellist.rds'))

# plot.gam plots
for (outcomeclass in outcomeclasses) {
  pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outcomeclass),'output_unadjusted/%1/coefficientplots/plotgam/logoddsscale/%2.pdf'))
  plot(modellist[[outcomeclass]], pages = 1, all.terms = TRUE,residuals = FALSE,shade=TRUE)
  dev.off()
  pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=outcomeclass),'output_unadjusted/%1/coefficientplots/plotgam/probscale/%2.pdf'))
  plot(modellist[[outcomeclass]], pages = 1, trans = plogis, shift = coef(modellist[[outcomeclass]])[1], all.terms = TRUE,residuals = FALSE,shade=TRUE)
  dev.off()
}

# save model summaries to file
pvaluesfile<-gsub('%s',modelname,'output_unadjusted/%s/pvalues/pvalues.txt')
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
write.table(interceptsdf,file=gsub('%s',modelname,'output_unadjusted/%s/pvalues/intercepts.txt'),col.names = TRUE,row.names = FALSE,sep='\t')


# ---------------------------
# MODELLING/PLOTTING USING MGCVIZ

smoothplotlist1<-list()  # logodds
smoothplotlist2<-list()  # probability
smoothplotlist3<-list()  # probability; y-axis lims=0-1

if (modelname=='log10PlasmidSize') {
  lowerlim<--5
  upperlim<-5
  upperlim_prob<-0.4
}
if (modelname=='InsertionSequenceDensity') {
  lowerlim<--4
  upperlim<-6
  upperlim_prob<-0.4
}
if (modelname=='NumOtherResistanceClasses') {
  lowerlim<--3
  upperlim<-9
  upperlim_prob<-1
}
if (modelname=='CollectionDate') {
  lowerlim<--2
  upperlim<-10
  upperlim_prob<-0.4
}

#create list of smooth plots to be combined in a grid
smoothplotlistnest1<-list()
smoothplotlistnest2<-list()
smoothplotlistnest3<-list()
for (j in 1:length(outcomeclasses)) {
  o<-getViz(modellist[[j]]) #this produces gamViz object for a given class, from which smooth can be extracted and plotted using sm(o,n)
  gamcoef<-coef(modellist[[j]])[1] #intercept
  o1 <- plot.mgcv.smooth.1D( sm(o, 1) ) #for resistance class j, extract smooth 
  o1<-o1 + geom_hline(yintercept = 0,linetype='solid',colour='light grey',size=0.3) + l_ciPoly() +l_fitLine(colour = brewerpal[j])
  o1<-o1 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(lowerlim,upperlim)) +
    scale_x_continuous(breaks = scales::pretty_breaks())
  o2 <- plot.mgcv.smooth.1D( sm(o, 1), trans = function(x) plogis(x+gamcoef))
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
  o1<-o1 + ggtitle(gsub('.','-',outcomeclasses[j],fixed=TRUE)) + 
    annotate('text',x=-Inf,y=Inf,label=paste('edf =',edf,collapse=''),hjust=-.1,vjust=1.8,colour="#525252",size=4)
  o2<-o2 + ggtitle(gsub('.','-',outcomeclasses[j],fixed=TRUE)) + 
    annotate('text',x=-Inf,y=Inf,label=paste('edf =',edf,collapse=''),hjust=-.1,vjust=1.8,colour="#525252",size=4)
  #further customise plots
  if (j == 1) {
    o1<-o1 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-12,-4),'pt'),axis.text=element_text(size=rel(1.1)))
    o2<-o2 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-12,-4),'pt'),axis.text=element_text(size=rel(1.1)))
  } else if (j == 6) {
    o1<-o1 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-4,-4),'pt'),axis.text=element_text(size=rel(1.1)))
    o2<-o2 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-4,-4),'pt'),axis.text=element_text(size=rel(1.1)))
  } else if (j == 5) {
    o1<-o1 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,12,-12,-14),'pt'),axis.text=element_text(size=rel(1.1)))
    o2<-o2 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,12,-12,-14),'pt'),axis.text=element_text(size=rel(1.1)))
  } else if (j == 10) {
    o1<-o1 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,12,-4,-14),'pt'),axis.text=element_text(size=rel(1.1)))
    o2<-o2 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,12,-4,-14),'pt'),axis.text=element_text(size=rel(1.1)))
  } else if (j > 6 && j < 10) {
    o1<-o1 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-4,-14),'pt'),axis.text=element_text(size=rel(1.1)))
    o2<-o2 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-4,-14),'pt'),axis.text=element_text(size=rel(1.1)))
  } else {
    o1<-o1 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-12,-14),'pt'),axis.text=element_text(size=rel(1.1)))
    o2<-o2 + xlab('') + ylab('') + theme(plot.title=element_text(hjust=0,size=13,colour='#525252'),plot.margin = unit(c(2,4,-12,-14),'pt'),axis.text=element_text(size=rel(1.1)))
  }
  #adjust probability y axis scale
  o3<-o2
  o3<-o3 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,1))
  o2<-o2 + scale_y_continuous(breaks = scales::pretty_breaks(n = 10),limits=c(0,upperlim_prob))
  #adjust x axis scales for plasmid size and collection date
  if (modelname=='log10PlasmidSize') {
    o1<-o1 + scale_x_continuous(labels = log10size_tokb, breaks = c(-0.522,0,0.5444,1,1.478))
    o2<-o2 + scale_x_continuous(labels = log10size_tokb, breaks = c(-0.522,0,0.5444,1,1.478))
    o3<-o3 + scale_x_continuous(labels = log10size_tokb, breaks = c(-0.522,0,0.5444,1,1.478))
  }
  if (modelname=='CollectionDate') {
    o1<-o1 + scale_x_continuous(labels = yearsince_tocollectionyear) + theme(axis.text.x = element_text(angle = 30, hjust=1,vjust=1))
    o2<-o2 + scale_x_continuous(labels = yearsince_tocollectionyear) + theme(axis.text.x = element_text(angle = 30, hjust=1,vjust=1))
    o3<-o3 + scale_x_continuous(labels = yearsince_tocollectionyear) + theme(axis.text.x = element_text(angle = 30, hjust=1,vjust=1))
  }
  #assign plots to list
  smoothplotlistnest1[[outcomeclasses[j]]]<-o1 #logodds
  smoothplotlistnest2[[outcomeclasses[j]]]<-o2 #prob
  smoothplotlistnest3[[outcomeclasses[j]]]<-o3
}

smoothplotlist1[[modelname]]<-smoothplotlistnest1
smoothplotlist2[[modelname]]<-smoothplotlistnest2
smoothplotlist3[[modelname]]<-smoothplotlistnest3



# plotting; save plot lists to file; save plots to pdf

# save lists to file
saveRDS(smoothplotlist1,file = gsub('%s',modelname,'output_unadjusted/%s/coefficientplots/logoddsscale/smoothplotlist1.rds'))
saveRDS(smoothplotlist2,file = gsub('%s',modelname,'output_unadjusted/%s/coefficientplots/probscale/smoothplotlist2.rds'))
saveRDS(smoothplotlist3,file = gsub('%s',modelname,'output_unadjusted/%s/coefficientplots/probscale/smoothplotlist3.rds'))

# smooth plots
width=12
height=5.8
probwidth=12.7  # need to adjust for longer ylab text


pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=modelname),'output_unadjusted/%1/coefficientplots/logoddsscale/%2.pdf'),width=width,height=height)
gridPrint(grobs=(smoothplotlist1[[modelname]]),nrow=2,bottom=textGrob(myxlab,gp=gpar(fontsize=15), vjust=-0.1),left=textGrob('Effect on log odds',gp=gpar(fontsize=15),rot=90, vjust=1),top = grid::textGrob(myggtitle, x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=15,lineheight=1)))
dev.off()
pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=modelname),'output_unadjusted/%1/coefficientplots/probscale/%2.pdf'),width=probwidth,height=height)
gridPrint(grobs=(smoothplotlist2[[modelname]]),nrow=2,bottom=textGrob(myxlab,gp=gpar(fontsize=15), vjust=-0.1),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90, vjust=1),top = grid::textGrob(myggtitle, x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=15,lineheight=1)))
dev.off()
pdf(gsubfn('%1|%2',list('%1'=modelname,'%2'=modelname),'output_unadjusted/%1/coefficientplots/probscale/ylim0to1_%2.pdf'),width=probwidth,height=height)
gridPrint(grobs=(smoothplotlist3[[modelname]]),nrow=2,bottom=textGrob(myxlab,gp=gpar(fontsize=15), vjust=-0.1),left=textGrob('Effect on predicted probability',gp=gpar(fontsize=15),rot=90, vjust=1),top = grid::textGrob(myggtitle, x = 0, hjust = 0, vjust=0.75, gp=gpar(fontsize=15,lineheight=1)))
dev.off()


