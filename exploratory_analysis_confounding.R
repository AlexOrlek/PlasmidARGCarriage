pacman::p_load(tidyverse,ggpubr,cowplot,gridExtra)
theme_set(theme_bw() + theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()))
dir.create('output_exploratory_confounding', showWarnings = FALSE)

# ---------------------------
# Additional exploratory analyses, following on from assocaition stats, investiating specific confounding associations
finaldftrunc<-read.table('data/plasmiddf_transformed.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")
finaldftrunc <- finaldftrunc %>% rename(`outcomeTEM-1` = outcomeTEM.1, `NumOtherResistanceClassesTEM-1` = NumOtherResistanceClassesTEM.1)
outcomeclasses<-c('aminoglycoside','phenicol','sulphonamide','tetracycline','macrolide','TEM-1','trimethoprim','ESBL', 'carbapenem','quinolone','colistin')

NumOtherResistanceClassesvec<-vector()
for (outcomeclass in outcomeclasses) {
  NumOtherResistanceClasses<-gsub('%s',outcomeclass,'NumOtherResistanceClasses%s')
  NumOtherResistanceClassesvec<-c(NumOtherResistanceClassesvec,NumOtherResistanceClasses)
}


levels(finaldftrunc$HostTaxonomy)<-c('Enterobacteriaceae','Firmicutes','other','Proteobacteria (non-Enterobacteriaceae)')
finaldftrunc$HostTaxonomy<-factor(finaldftrunc$HostTaxonomy,levels=c('Enterobacteriaceae','Proteobacteria (non-Enterobacteriaceae)','Firmicutes','other'))


# ---------------------------
# HostTaxonomy vs RepliconCarriage / NumOtherResistanceClasses (grouped barplots)
#table(finaldftrunc$HostTaxonomy)
#table(finaldftrunc$RepliconCarriage)

# HostTaxonomy vs RepliconCarriage
finaldftrunc$RepliconCarriage<-factor(finaldftrunc$RepliconCarriage,levels=c('untyped','single-replicon','multi-replicon'))
tax_rep_contingency<-as.data.frame(table(finaldftrunc$HostTaxonomy,finaldftrunc$RepliconCarriage))
colnames(tax_rep_contingency)<-c('HostTaxonomy','RepliconCarriage','Freq')
p<-ggplot(tax_rep_contingency, aes(x=HostTaxonomy,y=Freq,fill=RepliconCarriage)) + geom_bar(position = 'dodge',stat='identity') + labs(x='Host taxonomy',y='Frequency') + theme(axis.text.x=element_text(angle=30,hjust=1),legend.position = 'top',plot.title = element_text(hjust=0.5,size=11)) + scale_x_discrete(limits=c("Enterobacteriaceae", "Proteobacteria (non-Enterobacteriaceae)", "Firmicutes","other"), labels=c("Enterobacteriaceae", "Proteobacteria\n(non-Enterobacteriaceae)", "Firmicutes","other")) + scale_fill_discrete(labels=c('untyped','single-replicon','multi-replicon')) + ggtitle('Replicon carriage:') + guides(title=NULL,fill=guide_legend(title=NULL,nrow=1))

pdf('output_exploratory_confounding/HostTaxonomy_RepliconCarriage.pdf',4.6,3.3)
p
dev.off()

# HostTaxonomy vs NumOtherResistanceClasses (grouped barplots)
plotlist<-list()
for (i in 1:length(outcomeclasses)) {
  outcomeclass<-outcomeclasses[i]
  NumOtherResistanceClasses<-NumOtherResistanceClassesvec[i]
  tax_NumOtherResistanceClasses_contingency<-as.data.frame(table(finaldftrunc$HostTaxonomy,finaldftrunc[,NumOtherResistanceClasses]))
  colnames(tax_NumOtherResistanceClasses_contingency)<-c('HostTaxonomy','NumOtherResistanceClasses','Freq')
  p<-ggplot(tax_NumOtherResistanceClasses_contingency, aes(x=HostTaxonomy,y=Freq,fill=NumOtherResistanceClasses)) + geom_bar(position = 'dodge',stat='identity') + theme(axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),legend.position = 'top',legend.key.size=unit(0.3,'cm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-8,0),plot.title = element_text(hjust=0.5,size=9),plot.margin = unit(c(0.1,0.1,0.1,1),'cm')) + scale_x_discrete(limits=c("Enterobacteriaceae", "Proteobacteria (non-Enterobacteriaceae)", "Firmicutes","other"),labels=c("Enterobacteriaceae", "Proteobacteria\n(non-Enterobacteriaceae)", "Firmicutes","other")) + ggtitle(gsub('%s',outcomeclass,'Resistance class: %s\nOther resistance gene classes:')) + guides(fill=guide_legend(title=NULL,nrow=1))
  plotlist[[i]]<-p
}

pdf('output_exploratory_confounding/HostTaxonomy_NumOtherResistanceClasses.pdf',10,12)
grid.arrange(grobs=plotlist,bottom='Host taxonomy',left='Frequency',ncol=3)
dev.off()


# ---------------------------
# Integron vs BiocideMetalResistance vs NumOtherResistanceClasses (grouped barplots)
#Integron vs BiocideMetalResistance
Integron_BiocideMetalResistance_contingency<-as.data.frame(table(finaldftrunc$Integron,finaldftrunc$BiocideMetalResistance))
colnames(Integron_BiocideMetalResistance_contingency)<-c('Integron','BiocideMetalResistance','Freq')
p<-ggplot(Integron_BiocideMetalResistance_contingency, aes(x=Integron,y=Freq,fill=BiocideMetalResistance)) + geom_bar(position = 'dodge',stat='identity') + labs(x='Integron presence',y='Frequency') + theme(axis.text.x=element_text(angle=30,hjust=1),legend.position = 'top',plot.title = element_text(hjust=0.5,size=11)) + scale_x_discrete(labels=c("absence", "presence")) + scale_fill_discrete(labels=c('absence','presence')) + ggtitle('Biocide/metal resistance gene presence:') + guides(fill=guide_legend(title=NULL,nrow=1))
pdf('output_exploratory_confounding/Integron_BiocideMetalResistance.pdf',4,4)
p
dev.off()

#Integron vs NumOtherResistanceClasses
plotlist<-list()
for (i in 1:length(outcomeclasses)) {
  outcomeclass<-outcomeclasses[i]
  NumOtherResistanceClasses<-NumOtherResistanceClassesvec[i]
  Integron_NumOtherResistanceClasses_contingency<-as.data.frame(table(finaldftrunc$Integron,finaldftrunc[,NumOtherResistanceClasses]))
  colnames(Integron_NumOtherResistanceClasses_contingency)<-c('Integron','NumOtherResistanceClasses','Freq')
  p<-ggplot(Integron_NumOtherResistanceClasses_contingency, aes(x=Integron,y=Freq,fill=NumOtherResistanceClasses)) + geom_bar(position = 'dodge',stat='identity') + theme(axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),legend.position = 'top',legend.key.size=unit(0.3,'cm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-8,0),plot.title = element_text(hjust=0.5,size=9),plot.margin = unit(c(0.5,0.2,0.5,1),'cm')) + scale_x_discrete(labels=c('absence','presence')) + ggtitle(gsub('%s',outcomeclass,'Resistance class: %s\nOther resistance gene classes:')) + guides(fill=guide_legend(title=NULL,nrow=1))
  plotlist[[i]]<-p
}

pdf('output_exploratory_confounding/Integron_NumOtherResistanceClasses.pdf',10,12)
grid.arrange(grobs=plotlist,bottom='Integron presence',left='Frequency',ncol=3)
dev.off()

#BiocideMetalResistance vs NumOtherResistanceClasses
plotlist<-list()
for (i in 1:length(outcomeclasses)) {
  outcomeclass<-outcomeclasses[i]
  NumOtherResistanceClasses<-NumOtherResistanceClassesvec[i]
  BiocideMetalResistance_NumOtherResistanceClasses_contingency<-as.data.frame(table(finaldftrunc$BiocideMetalResistance,finaldftrunc[,NumOtherResistanceClasses]))
  colnames(BiocideMetalResistance_NumOtherResistanceClasses_contingency)<-c('BiocideMetalResistance','NumOtherResistanceClasses','Freq')
  p<-ggplot(BiocideMetalResistance_NumOtherResistanceClasses_contingency, aes(x=BiocideMetalResistance,y=Freq,fill=NumOtherResistanceClasses)) + geom_bar(position = 'dodge',stat='identity') + theme(axis.text.x=element_text(angle=30,hjust=1),axis.title=element_blank(),legend.position = 'top',legend.key.size=unit(0.3,'cm'),legend.margin=margin(0,0,0,0),legend.box.margin=margin(0,0,-8,0),plot.title = element_text(hjust=0.5,size=9),plot.margin = unit(c(0.5,0.2,0.5,1),'cm')) + scale_x_discrete(labels=c('absence','presence')) + ggtitle(gsub('%s',outcomeclass,'Resistance class: %s\nOther resistance gene classes:')) + guides(fill=guide_legend(title=NULL,nrow=1))
  plotlist[[i]]<-p
}

pdf('output_exploratory_confounding/BiocideMetalResistance_NumOtherResistanceClasses.pdf',10,12)
grid.arrange(grobs=plotlist,bottom='Biocide/metal resistance gene presence',left='Frequency',ncol=3)
dev.off()


# Notes
#https://people.ok.ubc.ca/jpither/modules/Visualizing_association_two_variables.html#Visualizing_association_between_a_numeric_and_a_categorical_variable
