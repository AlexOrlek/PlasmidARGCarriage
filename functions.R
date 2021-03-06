pacman::p_load(ggplot2, reshape2)

# ---------------------------
# exploratory analysis functions

# calculate antibiotic resistance gene totals at gene and plasmid level
get_resgene_data <- function(plasmiddf, outcomeclasses) {
  resgenedf<-plasmiddf[,c('TotalResGenes',outcomeclasses)]
  genestotal<-as.data.frame(rev(sort(colSums(resgenedf))))
  # total resistance genes by class
  colnames(genestotal)<-'Total resistance genes by class'
  # binary totals (i.e. per-plasmid presence/absence) of resistance genes by class
  resgenedfbinary<-resgenedf
  resgenedfbinary[resgenedfbinary>0]<-1
  plasmidstotal<-as.data.frame(rev(sort(colSums(resgenedfbinary))))
  colnames(plasmidstotal)<-'Total plasmids with a resistance gene by class'
  resgene_data<-list(genestotal, plasmidstotal, resgenedf, resgenedfbinary)
  names(resgene_data)<-c('total_genes', 'total_plasmids', 'resgenedf', 'resgenedfbinary')
  return(resgene_data)
}

# cross-tabulate categorical variables
crosstabulate<-function(plasmiddf,outcomeclass,factorvar) {
  plasmidf_crosstabulate<-plasmiddf
  outcomeclassbinary<-plasmidf_crosstabulate[,gsub('%s',outcomeclass,'outcome%s')]
  crosstable<-table(outcomeclassbinary,as.factor(plasmidf_crosstabulate[,factorvar]))
  FactorLevel<-colnames(crosstable)
  Count_NonResistancePlasmids<-crosstable[1,]
  Count_ResistancePlasmids<-crosstable[2,]
  ResistanceClass<-rep(outcomeclass,length(FactorLevel))
  FactorVariable<-rep(factorvar,length(FactorLevel))
  crosstable<-cbind(ResistanceClass,FactorVariable,FactorLevel,Count_NonResistancePlasmids,Count_ResistancePlasmids)
  rownames(crosstable)<-NULL
  return(crosstable)
}

# heatmap of resistance gene class intersections
#http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization #used correlation heat map code
#https://stackoverflow.com/questions/34659782/creation-of-a-contingency-table-based-on-two-columns-of-a-binary-matrix
resclassheatmap<-function(mymatrixinput,resclasses,twoway=TRUE,onewayplotlower=TRUE,hmtype='counts',hmtitle='') {
  # convert binary matrix to proportions
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
  # convert matrix to long dataframe  
  melted_hm <- reshape2::melt(mypropmatrix, na.rm = TRUE)
  # plot heatmap
  if (hmtype=='counts') {
    mygeomtext<-as.character(resgenedfbinarymat)
    mygeomtext[mygeomtext=='0']<-''
    #add totals to diagonal
    mygeomtext[seq(1,length(mygeomtext),nrow(resgenedfbinarymat)+1)]<-colSums(mymatrixinput[,resclasses])
    # plot
    ggplot(data = melted_hm, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
      scale_fill_gradientn(colours=c("white", "light blue", "red"), limit = c(0,max(mypropmatrix)), space = "Lab", name="Proportion") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),plot.title = element_text(size=11)) + coord_fixed() +ylab("") + xlab("") + geom_text(aes(label = mygeomtext)) + ggtitle(hmtitle)
  }
  else if (hmtype=='prop') {
    mygeomtext<-as.character(round(mypropmatrix,2))
    mygeomtext[mygeomtext=='0']<-''
    ggplot(data = melted_hm, aes(Var2, Var1, fill = value)) + geom_tile(color = "white") +
      scale_fill_gradientn(colours=c("white", "light blue", "red"), limit = c(0,max(mypropmatrix)), space = "Lab", name="Proportion") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1),plot.title = element_text(size=11)) + coord_fixed() +ylab("") + xlab("") + geom_text(aes(label = mygeomtext)) + ggtitle(hmtitle)
  }
  else {
    print('hmtype unrecognised')
  }
}


# ---------------------------
# exploratory adjusted modelling functions

writegamcheck<-function(model,outcomeclass,gamcheckfile,gamcheckplotfile) {
  if (is.null(gamcheckfile)) {
    #print gam.check output to terminal
    gcouts<-capture.output(gam.check(model[[outcomeclass]]))
    for (gcout in gcouts) {
      cat(gcout,sep='\n')
    }
  } else {
    #write gam.check output to file
    pdf(gamcheckplotfile)
    gcouts<-capture.output(gam.check(model[[outcomeclass]]))
    dev.off()
    for (gcout in gcouts) {
      if (grepl('^s\\(',gcout)==TRUE || grepl('^\\s+k',gcout)==TRUE) {
        gcout<-gsub("*\\s+", "\t", gcout)
      }
      cat(gcout,sep='\n',append = TRUE,file = gamcheckfile)
    }
  }
}

writeconcurvity<-function(model,outcomeclass,concurvityfile) {
  couts<-capture.output(concurvity(model[[outcomeclass]],full=TRUE))
  for (cout in couts) {
    cout<-gsub("*\\s+", "\t", cout)
    cat(cout,sep='\n',append = TRUE,file = concurvityfile)
  }
  cat('',sep='\n',append = TRUE,file = concurvityfile)
  couts<-capture.output(concurvity(model[[outcomeclass]],full=FALSE))
  for (cout in couts) {
    cout<-gsub("*\\s+", "\t", cout)
    cat(cout,sep='\n',append = TRUE,file = concurvityfile)
  }
}


# ---------------------------
# unadjusted and adjusted analysis functions

# formatter functions for axis tick labels
log10size_tokb <- function(x) {
  x <- x + 4
  x <- 10 ** x
  x <- round(floor(x/1000))
  return(x)
}

yearsince_tocollectionyear <- function(x, referenceyr = 1994) {
  x <- round(as.integer(x) + as.integer(referenceyr))
  return(x)
} 



