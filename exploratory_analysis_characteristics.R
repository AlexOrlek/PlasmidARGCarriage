setwd("C:/Users/alexo/Dropbox/Oxford_2015/thesis/CORRECTIONS/Chapter_5/resistance_prediction")


finaldftrunc<-read.table('data/Appendix_D_Table_6_splitbetalactam_truncated.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")
finaldftrunc$logplasmidsize_uncentred<-log10(finaldftrunc$plasmidsize)
colnames(finaldftrunc)
finaldftrunc<-finaldftrunc[,c('logplasmidsize_uncentred','isscore','coldateimputed','integronbinary','bacmetbinary','vfdbbinary','conj','reptypes','taxa','geographies','isolationsources')]  # can't include numotherresgenes
colnames(finaldftrunc)<-c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)', 'Integron presence', 'Biocide/metal resistance gene presence', 'Virulence gene presence', 'Conjugative system', 'Replicon type multiplicity', 'Host taxonomy', 'Geographic location', 'Isolation source')

str(finaldftrunc)
summary(finaldftrunc)

#re-name and re-order factor levels; change binary variables from integer to factor
finaldftrunc$`Integron presence`<-factor(finaldftrunc$`Integron presence`)
levels(finaldftrunc$`Integron presence`)<-c('Absence','Presence')

finaldftrunc$`Biocide/metal resistance gene presence`<-factor(finaldftrunc$`Biocide/metal resistance gene presence`)
levels(finaldftrunc$`Biocide/metal resistance gene presence`)<-c('Absence','Presence')

finaldftrunc$`Virulence gene presence`<-factor(finaldftrunc$`Virulence gene presence`)
levels(finaldftrunc$`Virulence gene presence`)<-c('Absence','Presence')


levels(finaldftrunc$`Conjugative system`)<-c('Conjugative','Mobilisable','Non-mobilisable')
finaldftrunc$`Conjugative system`<-factor(finaldftrunc$`Conjugative system`,levels=c('Non-mobilisable','Mobilisable','Conjugative'))

levels(finaldftrunc$`Replicon type multiplicity`)<-c('Multi-replicon','Single-replicon','Untyped')
finaldftrunc$`Replicon type multiplicity`<-factor(finaldftrunc$`Replicon type multiplicity`,levels=c('Untyped','Single-replicon','Multi-replicon'))

levels(finaldftrunc$`Host taxonomy`)<-c('Enterobacteriaceae','Firmicutes','Other','Proteobacteria')
finaldftrunc$`Host taxonomy`<-factor(finaldftrunc$`Host taxonomy`,levels=c('Enterobacteriaceae','Proteobacteria','Firmicutes','Other'))

levels(finaldftrunc$`Geographic location`)<-c('China','EU','High-income','Middle-income','Other','United States')
finaldftrunc$`Geographic location`<-factor(finaldftrunc$`Geographic location`,levels=c('High-income','Middle-income','EU','China','United States','Other'))

levels(finaldftrunc$`Isolation source`)<-c('Human','Livestock','Other')
finaldftrunc$`Isolation source`<-factor(finaldftrunc$`Isolation source`,levels=c('Human','Livestock','Other'))



#make table of variable characteristics
summary(finaldftrunc)

library(tableone)
?tableone
?CreateTableOne
#Create a variable list which we want in Table 1
listVars<-c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)', 'Integron presence', 'Biocide/metal resistance gene presence', 'Virulence gene presence', 'Conjugative system', 'Replicon type multiplicity', 'Host taxonomy', 'Geographic location', 'Isolation source')

#Define categorical variables
catVars <- c('Integron presence', 'Biocide/metal resistance gene presence', 'Virulence gene presence', 'Conjugative system', 'Replicon type multiplicity', 'Host taxonomy', 'Geographic location', 'Isolation source')

table1 <- CreateTableOne(vars = listVars, data = finaldftrunc, factorVars = catVars)
print(table1,nonnormal = c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)'), quote = TRUE, noSpaces = TRUE)


#export table
tabMat <- print(table1, nonnormal = c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)'), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)
tabMat
## Save to a CSV file
write.table(tabMat, file = 'output_exploratory/characteristics.tsv',sep='\t')


#Notes
#https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html
#https://www.r-bloggers.com/2016/02/table-1-and-the-characteristics-of-study-population/

