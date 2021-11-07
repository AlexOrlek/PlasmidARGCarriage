pacman::p_load(tidyverse,tableone)

# load data
finaldftrunc<-read.table('data/plasmiddf_transformed.tsv',header=TRUE,sep='\t',stringsAsFactors = TRUE,quote = "",comment.char = "")
finaldftrunc$log10PlasmidSize_uncentred<-log10(finaldftrunc$PlasmidSize)
finaldftrunc$log10PlasmidSize <- NULL
#colnames(finaldftrunc)
finaldftrunc <- finaldftrunc %>% rename(`log10 Plasmid size` = log10PlasmidSize_uncentred, `Insertion sequence density (frequency per 10 kb)` = InsertionSequenceDensity, `Collection date (years since initial collection year)` = CollectionDate, `Integron presence` = Integron, `Biocide/metal resistance gene presence` = BiocideMetalResistance, `Virulence gene presence` = Virulence, `Conjugative system` = ConjugativeSystem, `Replicon carriage` = RepliconCarriage, `Host taxonomy` = HostTaxonomy, `Geographic location` = GeographicLocation, `Isolation source` = IsolationSource)


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

levels(finaldftrunc$`Host taxonomy`)<-c('Enterobacteriaceae','Firmicutes','Other','Proteobacteria (non-Enterobacteriaceae)')
finaldftrunc$`Host taxonomy`<-factor(finaldftrunc$`Host taxonomy`,levels=c('Enterobacteriaceae','Proteobacteria (non-Enterobacteriaceae)','Firmicutes','Other'))

levels(finaldftrunc$`Geographic location`)<-c('China','EU','High-income','Middle-income','Other','United States')
finaldftrunc$`Geographic location`<-factor(finaldftrunc$`Geographic location`,levels=c('High-income','Middle-income','EU','China','United States','Other'))

levels(finaldftrunc$`Isolation source`)<-c('Human','Livestock','Other')
finaldftrunc$`Isolation source`<-factor(finaldftrunc$`Isolation source`,levels=c('Human','Livestock','Other'))


# split data into resistant / non-resistant plasmids
finaldftrunc_resistant <- finaldftrunc[finaldftrunc %>% select(starts_with('outcome')) %>% rowSums() > 0,]
finaldftrunc_nonresistant <- finaldftrunc[finaldftrunc %>% select(starts_with('outcome')) %>% rowSums() == 0,]
nrow(finaldftrunc_resistant)  # 3692


# summarise variable characteristics
#summary(finaldftrunc)

# define all categorical variables
listVars<-c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)', 'Integron presence', 'Biocide/metal resistance gene presence', 'Virulence gene presence', 'Conjugative system', 'Replicon carriage', 'Host taxonomy', 'Geographic location', 'Isolation source')

# define factor variables
catVars <- c('Integron presence', 'Biocide/metal resistance gene presence', 'Virulence gene presence', 'Conjugative system', 'Replicon carriage', 'Host taxonomy', 'Geographic location', 'Isolation source')

# define nonnormal variables
nonnormal = c('log10 Plasmid size', 'Insertion sequence density (frequency per 10 kb)', 'Collection date (years since initial collection year)')

# ------- all plasmids
table1 <- CreateTableOne(vars = listVars, data = finaldftrunc, factorVars = catVars)

# export table
tabMat <- print(table1, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

# save to file
write.table(tabMat, file = 'output_exploratory/characteristics_allplasmids.tsv',sep='\t')

# ------- resistant
table1_resistant <- CreateTableOne(vars = listVars, data = finaldftrunc_resistant, factorVars = catVars)

# export table
tabMat_resistant <- print(table1_resistant, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

# save to file
write.table(tabMat_resistant, file = 'output_exploratory/characteristics_resistant.tsv',sep='\t')

# ------ non-resistant
table1_nonresistant <- CreateTableOne(vars = listVars, data = finaldftrunc_nonresistant, factorVars = catVars)

# export table
tabMat_nonresistant <- print(table1_nonresistant, nonnormal = nonnormal, quote = FALSE, noSpaces = TRUE, printToggle = FALSE)

# save to file
write.table(tabMat_nonresistant, file = 'output_exploratory/characteristics_nonresistant.tsv',sep='\t')


# Notes
#https://cran.r-project.org/web/packages/tableone/vignettes/introduction.html
#https://www.r-bloggers.com/2016/02/table-1-and-the-characteristics-of-study-population/

