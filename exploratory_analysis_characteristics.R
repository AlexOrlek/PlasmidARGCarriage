pacman::p_load(tableone)

# load data
finaldftrunc<-read.table('data/plasmiddf_transformed.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")
finaldftrunc$log10PlasmidSize_uncentred<-log10(finaldftrunc$PlasmidSize)
#colnames(finaldftrunc)
finaldftrunc<-finaldftrunc[,c('log10PlasmidSize_uncentred','InsertionSequenceDensity','CollectionDate','Integron','BiocideMetalResistance','Virulence','ConjugativeSystem','RepliconCarriage','HostTaxonomy','GeographicLocation','IsolationSource')]  # can't include numotherresgenes
colnames(finaldftrunc)<-c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)', 'Integron presence', 'Biocide/metal resistance gene presence', 'Virulence gene presence', 'Conjugative system', 'Replicon carriage', 'Host taxonomy', 'Geographic location', 'Isolation source')


# re-name and re-order factor levels; change binary variables from 0/1 to absence/presence

# binary variables
finaldftrunc$`Integron presence`<-factor(finaldftrunc$`Integron presence`)
levels(finaldftrunc$`Integron presence`)<-c('Absence','Presence')

finaldftrunc$`Biocide/metal resistance gene presence`<-factor(finaldftrunc$`Biocide/metal resistance gene presence`)
levels(finaldftrunc$`Biocide/metal resistance gene presence`)<-c('Absence','Presence')

finaldftrunc$`Virulence gene presence`<-factor(finaldftrunc$`Virulence gene presence`)
levels(finaldftrunc$`Virulence gene presence`)<-c('Absence','Presence')


# factor variables
levels(finaldftrunc$`Conjugative system`)<-c('Conjugative','Mobilisable','Non-mobilisable')
finaldftrunc$`Conjugative system`<-factor(finaldftrunc$`Conjugative system`,levels=c('Non-mobilisable','Mobilisable','Conjugative'))

levels(finaldftrunc$`Replicon carriage`)<-c('Multi-replicon','Single-replicon','Untyped')
finaldftrunc$`Replicon carriage`<-factor(finaldftrunc$`Replicon carriage`,levels=c('Untyped','Single-replicon','Multi-replicon'))

levels(finaldftrunc$`Host taxonomy`)<-c('Enterobacteriaceae','Firmicutes','Other','Proteobacteria_other')
finaldftrunc$`Host taxonomy`<-factor(finaldftrunc$`Host taxonomy`,levels=c('Enterobacteriaceae','Proteobacteria_other','Firmicutes','Other'))

levels(finaldftrunc$`Geographic location`)<-c('China','EU','High-income','Middle-income','Other','United States')
finaldftrunc$`Geographic location`<-factor(finaldftrunc$`Geographic location`,levels=c('High-income','Middle-income','EU','China','United States','Other'))

levels(finaldftrunc$`Isolation source`)<-c('Human','Livestock','Other')
finaldftrunc$`Isolation source`<-factor(finaldftrunc$`Isolation source`,levels=c('Human','Livestock','Other'))


# summarise variable characteristics
#summary(finaldftrunc)

# define all categorical variables
listVars<-c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)', 'Integron presence', 'Biocide/metal resistance gene presence', 'Virulence gene presence', 'Conjugative system', 'Replicon carriage', 'Host taxonomy', 'Geographic location', 'Isolation source')

# define factor variables
catVars <- c('Integron presence', 'Biocide/metal resistance gene presence', 'Virulence gene presence', 'Conjugative system', 'Replicon carriage', 'Host taxonomy', 'Geographic location', 'Isolation source')

table1 <- CreateTableOne(vars = listVars, data = finaldftrunc, factorVars = catVars)
#print(table1,nonnormal = c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)'), quote = TRUE, noSpaces = TRUE)

# export table
tabMat <- print(table1, nonnormal = c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)'), quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

# save to file
write.table(tabMat, file = 'output_exploratory/characteristics.tsv',sep='\t')


# Notes
#https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html
#https://www.r-bloggers.com/2016/02/table-1-and-the-characteristics-of-study-population/

