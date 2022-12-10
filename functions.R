pacman::p_load(ggplot2, reshape2, stringr)

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

# heatmap functions for visualising resistance gene class intersections
simililarity_index = function(x, y) {
  intersection <- sum(x == 1 & y == 1)
  union <- sum(x == 1 | y == 1)
  min_xy <- min(sum(x == 1), sum(y == 1))
  jaccard_similarity <- intersection / union
  overlap_coefficient <- intersection / min_xy
  return(list("jaccard_similarity" = jaccard_similarity, "overlap_coefficient" = overlap_coefficient, "intersection" = intersection, "union" = union, "min_xy" = min_xy))
}

resclassheatmap <- function(mymatrixinput, resclasses) {
  mymatrixinput <- mymatrixinput[,resclasses]
  # calculate similarity indices (jaccard and overlap coefficient) and count intersections
  m_jaccard <- matrix(data = NA, nrow = length(resclasses), ncol = length(resclasses))  # jaccard similarity matrix
  m_overlap_coefficient <- matrix(data = NA, nrow = length(resclasses), ncol = length(resclasses))  # overlap coefficient matrix
  m_intersection <- matrix(data = NA, nrow = length(resclasses), ncol = length(resclasses))  # count matrix of intersections
  for (i in 1:length(resclasses)) {  # rows
    for (j in 1:length(resclasses)) {  # columns
      if (i > j) {  # lower left triangle
        m_jaccard[i,j] <- 0
        m_overlap_coefficient[i,j] <- 0
        m_intersection[i,j] <- NA
      }
      else if (i == j) {
        m_jaccard[i,j] <- 0
        m_overlap_coefficient[i,j] <- 0
        m_intersection[i,j] <- sum(mymatrixinput[,i])
      } else {  # j > i : column index > row index (fill out upper right triangle)
        output <- simililarity_index(mymatrixinput[,i], mymatrixinput[,j])
        m_jaccard[i,j] <- output$jaccard_similarity
        m_overlap_coefficient[i,j] <- output$overlap_coefficient
        m_intersection[i,j] <- output$intersection
      }
    }
  }
  # plot heatmps
  hm_text <- str_replace_na(as.character(m_intersection),'')
  
  # jaccard heatmap
  colnames(m_jaccard) <- outcomeclasses
  rownames(m_jaccard) <- outcomeclasses
  melted_m_jaccard <- reshape2::melt(m_jaccard, na.rm = TRUE)
  
  p_jaccard <- ggplot(data = melted_m_jaccard, aes(Var1, Var2, fill = value)) + geom_tile(color = "white") + scale_fill_gradientn(colours=c("white", "light blue", "red"), limit = c(0,max(m_jaccard)), name="Jaccard\nsimilarity coefficient") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(size=11)) + coord_fixed() + ylab("") + xlab("") + geom_text(aes(label = hm_text)) + ggtitle('')
  
  # overlap coefficient heatmap
  colnames(m_overlap_coefficient) <- outcomeclasses
  rownames(m_overlap_coefficient) <- outcomeclasses
  melted_m_sim <- reshape2::melt(m_overlap_coefficient, na.rm = TRUE)
  
  p_overlap_coefficient <- ggplot(data = melted_m_sim, aes(Var1, Var2, fill = value)) + geom_tile(color = "white") + scale_fill_gradientn(colours=c("white", "light blue", "red"), limit = c(0,max(m_overlap_coefficient)), name="Overlap coefficient") + theme_minimal() + theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1), plot.title = element_text(size=11)) + coord_fixed() + ylab("") + xlab("") + geom_text(aes(label = hm_text)) + ggtitle('')
  
  return(list('jaccard heatmap' = p_jaccard, 'overlap coefficient heatmap' = p_overlap_coefficient))
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



