library(ggplot2)
theme_set(theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
library(ggpubr)
library(cowplot)
library(gridExtra)

###Additional exploratory analyses, following on from assocaition stats, investiating specific confounding associations###
#https://people.ok.ubc.ca/jpither/modules/Visualizing_association_two_variables.html#Visualizing_association_between_a_numeric_and_a_categorical_variable

finaldftrunc<-read.table('data/Appendix_D_Table_6_splitbetalactam_truncated.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")
outcomeclasses<-c('aminoglycoside','betalactam_carbapenem','betalactam_ESBL','betalactam_other','macrolide','phenicol','quinolone','sulphonamide','tetracycline','trimethoprim')
numotherresgenesvec<-vector()
for (outcomeclass in outcomeclasses) {
  numotherresgenes<-gsub('%s',outcomeclass,'numotherresgenes%s')
  numotherresgenesvec<-c(numotherresgenesvec,numotherresgenes)
}

###Host taxonomy vs reptype / numothreresgenes (grouped barplots)
table(finaldftrunc$taxa)
table(finaldftrunc$reptypes)


#host taxonomy vs replicon
finaldftrunc$reptypes<-factor(finaldftrunc$reptypes,levels=c('untyped','single','multi'))
tax_rep_contingency<-as.data.frame(table(finaldftrunc$taxa,finaldftrunc$reptypes))
colnames(tax_rep_contingency)<-c('taxa','reptypes','Freq')
p<-ggplot(tax_rep_contingency, aes(x=taxa,y=Freq,fill=reptypes)) + geom_bar(position = 'dodge',stat='identity') + labs(x='Host taxonomy',y='Frequency') + theme(axis.text.x=element_text(angle=30,hjust=1),legend.position = 'top',plot.title = element_text(hjust=0.5,size=11)) + scale_x_discrete(limits=c("Enterobacteriaceae", "Proteobacteria", "Firmicutes","other")) + scale_fill_discrete(labels=c('untyped','single-replicon','multi-replicon')) + ggtitle('Replicon type multiplicity:') + guides(title=NULL,fill=guide_legend(title=NULL,nrow=1))
pdf('output_exploratory_confounding/taxa_reptypes.pdf',4.6,3.3)
p
dev.off()


#host taxonomy vs numotherresgenes
# #boxplots
# for (numotherresgenes in numotherresgenesvec) {
#   finaldftrunc[,'numotherresgenes']<-finaldftrunc[,numotherresgenes]
#   p<-ggplot(finaldftrunc,aes(x=taxa,y=numotherresgenes)) + geom_boxplot()
#   print(p)
# }

#grouped barplots
plotlist<-list()
for (i in 1:length(outcomeclasses)) {
  outcomeclass<-outcomeclasses[i]
  numotherresgenes<-numotherresgenesvec[i]
  tax_numotherresgenes_contingency<-as.data.frame(table(finaldftrunc$taxa,finaldftrunc[,numotherresgenes]))
  colnames(tax_numotherresgenes_contingency)<-c('taxa','numotherresgenes','Freq')
  p<-ggplot(tax_numotherresgenes_contingency, aes(x=taxa,y=Freq,fill=numotherresgenes)) + geom_bar(position = 'dodge',stat='identity') + theme(axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),legend.position = 'top',legend.key.size=unit(0.3,'cm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-8,0),plot.title = element_text(hjust=0.5,size=9),plot.margin = unit(c(0.1,0.1,0.1,1),'cm')) + scale_x_discrete(limits=c("Enterobacteriaceae", "Proteobacteria", "Firmicutes","other")) + ggtitle(gsub('%s',outcomeclass,'Resistance class: %s\nOther resistance gene classes:')) + guides(fill=guide_legend(title=NULL,nrow=1))
  plotlist[[i]]<-p
}

pdf('output_exploratory_confounding/taxa_numotherresgenes.pdf',10,12)
#plot_grid(plotlist=plotlist,nrow=2)
grid.arrange(grobs=plotlist,bottom='Host taxonomy',left='Frequency',ncol=3)
dev.off()




###integron vs bacmet vs numotherresgenes (grouped barplots)
#integron vs bacmet
integron_bacmet_contingency<-as.data.frame(table(finaldftrunc$integronbinary,finaldftrunc$bacmetbinary))
colnames(integron_bacmet_contingency)<-c('integronbinary','bacmetbinary','Freq')
p<-ggplot(integron_bacmet_contingency, aes(x=integronbinary,y=Freq,fill=bacmetbinary)) + geom_bar(position = 'dodge',stat='identity') + labs(x='Integron presence',y='Frequency') + theme(axis.text.x=element_text(angle=30,hjust=1),legend.position = 'top',plot.title = element_text(hjust=0.5,size=11)) + scale_x_discrete(labels=c("absence", "presence")) + scale_fill_discrete(labels=c('absence','presence')) + ggtitle('Biocide/metal resistance gene presence:') + guides(fill=guide_legend(title=NULL,nrow=1))
pdf('output_exploratory_confounding/integron_bacmet.pdf',4,4)
p
dev.off()

#integron vs numotherresgenes
plotlist<-list()
for (i in 1:length(outcomeclasses)) {
  outcomeclass<-outcomeclasses[i]
  numotherresgenes<-numotherresgenesvec[i]
  integron_numotherresgenes_contingency<-as.data.frame(table(finaldftrunc$integronbinary,finaldftrunc[,numotherresgenes]))
  colnames(integron_numotherresgenes_contingency)<-c('integronbinary','numotherresgenes','Freq')
  p<-ggplot(integron_numotherresgenes_contingency, aes(x=integronbinary,y=Freq,fill=numotherresgenes)) + geom_bar(position = 'dodge',stat='identity') + theme(axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),legend.position = 'top',legend.key.size=unit(0.3,'cm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-8,0),plot.title = element_text(hjust=0.5,size=9),plot.margin = unit(c(0.5,0.2,0.5,1),'cm')) + scale_x_discrete(labels=c('absence','presence')) + ggtitle(gsub('%s',outcomeclass,'Resistance class: %s\nOther resistance gene classes:')) + guides(fill=guide_legend(title=NULL,nrow=1))
  plotlist[[i]]<-p
}

pdf('output_exploratory_confounding/integron_numotherresgenes.pdf',10,12)
grid.arrange(grobs=plotlist,bottom='Integron presence',left='Frequency',ncol=3)
dev.off()

#bacmet vs numotherresgenes
plotlist<-list()
for (i in 1:length(outcomeclasses)) {
  outcomeclass<-outcomeclasses[i]
  numotherresgenes<-numotherresgenesvec[i]
  bacmet_numotherresgenes_contingency<-as.data.frame(table(finaldftrunc$bacmetbinary,finaldftrunc[,numotherresgenes]))
  colnames(bacmet_numotherresgenes_contingency)<-c('bacmetbinary','numotherresgenes','Freq')
  p<-ggplot(bacmet_numotherresgenes_contingency, aes(x=bacmetbinary,y=Freq,fill=numotherresgenes)) + geom_bar(position = 'dodge',stat='identity') + theme(axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),legend.position = 'top',legend.key.size=unit(0.3,'cm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-8,0),plot.title = element_text(hjust=0.5,size=9),plot.margin = unit(c(0.5,0.2,0.5,1),'cm')) + scale_x_discrete(labels=c('absence','presence')) + ggtitle(gsub('%s',outcomeclass,'Resistance class: %s\nOther resistance gene classes:')) + guides(fill=guide_legend(title=NULL,nrow=1))
  plotlist[[i]]<-p
}

pdf('output_exploratory_confounding/bacmet_numotherresgenes.pdf',10,12)
grid.arrange(grobs=plotlist,bottom='Biocide/metal resistance gene presence',left='Frequency',ncol=3)
dev.off()



###OLD
# plotlist<-list()
# for (i in 1:length(outcomeclasses)) {
#   outcomeclass<-outcomeclasses[i]
#   numotherresgenes<-numotherresgenesvec[i]
#   bacmet_numotherresgenes_contingency<-as.data.frame(table(finaldftrunc$bacmetbinary,finaldftrunc[,numotherresgenes]))
#   colnames(bacmet_numotherresgenes_contingency)<-c('bacmetbinary','numotherresgenes','Freq')
#   p<-ggplot(bacmet_numotherresgenes_contingency, aes(x=bacmetbinary,y=Freq,fill=numotherresgenes)) + geom_bar(position = 'dodge',stat='identity') + theme(axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),legend.position = 'top',legend.key.size=unit(0.3,'cm'),plot.title = element_text(hjust=0.5,size=9),plot.margin = unit(c(0.5,0.2,0.5,1),'cm')) + scale_x_discrete(labels=c('absence','presence')) + ggtitle(gsub('%s',outcomeclass,'Resistance class: %s\nOther resistance gene classes:')) + guides(fill=guide_legend(title=NULL,nrow=1))
#   plotlist[[i]]<-p
# }
# 
# pdf('output_exploratory_confounding/bacmet_numotherresgenes.pdf',12,10)
# grid.arrange(grobs=plotlist,bottom='Biocide/metal resistance gene presence',left='Frequency',ncol=4)
# dev.off()



