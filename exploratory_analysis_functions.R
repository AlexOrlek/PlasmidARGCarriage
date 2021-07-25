###exploratory analysis functions

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  rabs <- abs(cor(x, y, use = "complete.obs",method="spearman"))
  r <- cor(x, y, use = "complete.obs",method="spearman")
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste(prefix, txt, sep = "")
  if (missing(cex.cor)) cex.cor <- 0.8/strwidth(txt)
  text(0.5, 0.5, txt, cex =  cex.cor * (1 + rabs) / 2)
}

panel.hist <- function(x, ...) {
  usr <- par("usr")
  on.exit(par(usr))
  par(usr = c(usr[1:2], 0, 1.5) )
  h <- hist(x, plot = FALSE)
  breaks <- h$breaks
  nB <- length(breaks)
  y <- h$counts
  y <- y/max(y)
  rect(breaks[-nB], 0, breaks[-1], y, col = "white", ...)
}


crosstabulate<-function(outcomeclass,factorvar) {
  outcomeclassbinary<-finaldf[,gsub('%s',outcomeclass,'outcome%s')]
  crosstable<-table(outcomeclassbinary,as.factor(finaldf[,factorvar]))
  FactorLevel<-colnames(crosstable)
  Count_NonResistancePlasmids<-crosstable[1,]
  Count_ResistancePlasmids<-crosstable[2,]
  ResistanceClass<-rep(outcomeclass,length(FactorLevel))
  FactorVariable<-rep(factorvar,length(FactorLevel))
  crosstable<-cbind(ResistanceClass,FactorVariable,FactorLevel,Count_NonResistancePlasmids,Count_ResistancePlasmids)
  rownames(crosstable)<-NULL
  return(crosstable)
}


myproptable<-function(x,y,z,style='dataframelist',oddstable=TRUE) {  #kable doesn't produce odds table
  tablelistout<-list()
  library(dplyr)
  library(kableExtra)
  options(knitr.table.format = "markdown") 
  mytable<-with(finaldf, table(x, y))
  prop <- prop.table(mytable)
  prop2<-prop.table(mytable,margin=2)
  out<-cbind(mytable,prop,prop2)
  colnum<-as.integer(ncol(out)/3)
  tablelistout[[1]]<-z
  #print(colnames(prop))
  if (style=='dataframelist') {
    tablelistout[[2]]<-kable(out)
    if (oddstable==TRUE) {
      #get odds and odds ratio table + confidence interval of odds ratio
      oddsoriginal<-(out[2,1]/sum(out[,1]))/(out[1,1]/sum(out[,1]))
      olist<-list()
      orlist<-list()
      orCIlowerlist<-list()
      orCIupperlist<-list()
      olist[[1]]<-oddsoriginal
      orlist[[1]]<-NA
      orCIlowerlist[[1]]<-NA
      orCIupperlist[[1]]<-NA
      for (i in 1:ncol(prop)) {
        if (i==1) {
          next
        }
        odds<-(out[2,i]/sum(out[,i]))/(out[1,i]/sum(out[,i]))
        oddsratio<-odds/oddsoriginal
        olist[[i]]<-odds
        orlist[[i]]<-oddsratio
        oddsratioSE<-sqrt(sum(1/as.vector(mytable)))
        oddsratioCIlower<-exp(log(oddsratio) - (1.96*oddsratioSE))
        oddsratioCIupper<-exp(log(oddsratio) + (1.96*oddsratioSE))
        orCIlowerlist[[i]]<-oddsratioCIlower
        orCIupperlist[[i]]<-oddsratioCIupper
      }
      oddstable<-rbind(unlist(olist),unlist(orlist),unlist(orCIlowerlist),unlist(orCIupperlist))
      rownames(oddstable)<-c('odds','odds ratio','lower 95% CI','upper 95% CI')
      colnames(oddstable)<-colnames(prop)
      tablelistout[[3]]<-kable(oddstable)
      return(tablelistout)
    }
  }
  else if (style=='kable') {
    library(kableExtra)
    kable(out) %>%
      kable_styling("striped") %>%
      add_header_above(c(" " = 1,"Counts" = colnum, "Proportion of total" = colnum, "Columnwise proportion per category" = colnum))
  }
  else {
    print('unknown style')
  }
}

writeproptablelist<-function(proptablelist,filename,appendfile=TRUE) {
  for (i in 1:length(proptablelist)) {
    if (i==1) {
      cat(proptablelist[[i]],file=filename,append=appendfile,sep='\n')
    }
    if (i>1) {
      cat(gsub('|','\t',proptablelist[[i]],fixed=TRUE),sep='\n',file = filename,append=appendfile)
    }
  } 
}



#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization #used correlation heat map code
#https://stackoverflow.com/questions/34659782/creation-of-a-contingency-table-based-on-two-columns-of-a-binary-matrix
resclassheatmap<-function(mymatrixinput,resclasses,twoway=TRUE,onewayplotlower=TRUE,hmtext='counts') {
  library(ggplot2)
  library(reshape2)
  #convert binary matrix to proportions
  resgenedfbinarymat<-t(apply(mymatrixinput[,resclasses], 2, function(x) apply(mymatrixinput[,resclasses], 2, function(y) sum(x==1 & y ==1))))
  if (twoway==TRUE) {
    diag(resgenedfbinarymat)<-0
  } else {
    if (onewayplotlower==TRUE) {
      resgenedfbinarymat[lower.tri(resgenedfbinarymat)] <- 0
      resgenedfbinarymat[lower.tri(resgenedfbinarymat,diag=TRUE)] <- 0
    } else {
      resgenedfbinarymat[upper.tri(resgenedfbinarymat)] <- 0
      resgenedfbinarymat[upper.tri(resgenedfbinarymat,diag=TRUE)] <- 0
    }
  }
  mypropmatrix<-resgenedfbinarymat/colSums(mymatrixinput[,resclasses])
  #convert matrix to long dataframe  
  melted_hm <- reshape2::melt(mypropmatrix, na.rm = TRUE)
  # Heatmap
  if (hmtext=='counts') {
    mygeomtext<-as.character(resgenedfbinarymat)
    mygeomtext[mygeomtext=='0']<-''
    #add totals to diagonal
    mygeomtext[seq(1,length(mygeomtext),nrow(resgenedfbinarymat)+1)]<-colSums(mymatrixinput[,resclasses])
    #plot
    ggplot(data = melted_hm, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
      scale_fill_gradientn(colours=c("white", "light blue", "red"), limit = c(0,max(mypropmatrix)), space = "Lab", name="Proportion") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),plot.title = element_text(size=11)) + coord_fixed() +ylab("Resistance class") + xlab("") + geom_text(aes(label = mygeomtext)) + ggtitle('Resistance class intersect heatmap\n(counts indicated)')
  }
  else if (hmtext=='prop') {
    mygeomtext<-as.character(round(mypropmatrix,2))
    mygeomtext[mygeomtext=='0']<-''
    ggplot(data = melted_hm, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
      scale_fill_gradientn(colours=c("white", "light blue", "red"), limit = c(0,max(mypropmatrix)), space = "Lab", name="Proportion") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),plot.title = element_text(size=11)) + coord_fixed() +ylab("Resistance class") + xlab("") + geom_text(aes(label = mygeomtext)) + ggtitle('Resistance class intersect heatmap')
  }
  else {
    print('hmtext type unrecognised')
  }
}






#relationship between continuous predictors and outcome
contpredictorplot<-function(contpredictor,outcomeclassbinary,outcomeclass,contpredictorname) {
  library(ggplot2)
  overplottingdf<-data.frame(contpredictor,outcomeclassbinary)
  overplottingdf$outcomeclassbinary<-as.numeric(as.character(overplottingdf$outcomeclassbinary))
  sp<-ggplot(overplottingdf, aes(x=contpredictor, y=outcomeclassbinary))
  sp<-sp + stat_bin2d(bins=50)
  spdat<-ggplot_build(sp)
  spdat<-spdat[[1]][[1]][,'count']
  sp<-sp + scale_fill_gradientn(colours=c('white','cyan','green','yellow','orange','red'), limits=c(0, max(spdat))) + 
    theme(plot.title = element_text(size=10),panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + xlab(contpredictorname) + ylab(gsub('%s',outcomeclass,'%s binary')) + scale_y_continuous(breaks=c(0,1))
  return(sp)
}

