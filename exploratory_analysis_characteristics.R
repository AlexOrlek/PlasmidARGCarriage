pacman::p_load(tidyverse, tableone, rio, readxl)

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

levels(finaldftrunc$`Geographic location`)<-c('China','EU & UK','High-income','Middle-income','Other','United States')
finaldftrunc$`Geographic location`<-factor(finaldftrunc$`Geographic location`,levels=c('High-income','Middle-income','EU & UK','China','United States','Other'))

levels(finaldftrunc$`Isolation source`)<-c('Human','Livestock','Other')
finaldftrunc$`Isolation source`<-factor(finaldftrunc$`Isolation source`,levels=c('Human','Livestock','Other'))

# load and join additional data
resdf <- read_excel("data/Data_S1.xlsx", sheet = 'I')
resdf <- resdf %>% mutate(ARGName = str_replace_all(ARGName, "; ", ",")) %>% as.data.frame()
plasmiddf <- read_excel("data/Data_S1.xlsx", sheet = 'G')
plasmiddf <- plasmiddf %>% mutate(across(everything(), ~ str_replace_all(.x, "^\'|\'$", ""))) %>% as.data.frame()
colnames(plasmiddf) <- colnames(plasmiddf) %>% str_replace_all("^\'|\'$", "")

# check if all EU & UK countries are high-income
all_data <- finaldftrunc %>% left_join(plasmiddf, by = 'Accession') %>% left_join(resdf, by = 'Accession')
eu_uk_countries <- all_data %>% filter(`Geographic location` == 'EU & UK') %>% distinct(Country)
incomedf<-read.table('data/GlobalHealthObservatoryMetadata_edit.tsv',header=TRUE,sep='\t',as.is = TRUE,quote = "",comment.char = "")
eu_uk_countries %>% left_join(incomedf, by = c('Country' = 'DisplayString')) %>% select(Country, WORLD_BANK_INCOME_GROUP_CODE)  # all except Romania are high-income (in WB 2018 categorisation)

# split data into resistant / non-resistant plasmids and output resistant plasmid data for microreact
finaldftrunc_resistant <- finaldftrunc[finaldftrunc %>% select(starts_with('outcome')) %>% rowSums() > 0,]
finaldftrunc_nonresistant <- finaldftrunc[finaldftrunc %>% select(starts_with('outcome')) %>% rowSums() == 0,]
nrow(finaldftrunc_resistant)

microreact_data <- finaldftrunc_resistant %>% left_join(plasmiddf, by = 'Accession') %>% left_join(resdf, by = 'Accession')
microreact_data <- microreact_data %>% mutate(CollectionDate2 = CollectionDate) %>% separate(col = CollectionDate, into = c("year", "month", "day"), sep = "-", fill = "right") %>% rename(CollectionDate = CollectionDate2)
microreact_data <- microreact_data %>% mutate(LatitudeLongitude = ifelse(LatitudeLongitude == '-', '', LatitudeLongitude)) %>% separate(col = LatitudeLongitude, into = c("latitude", "longitude"), sep = ",") %>% mutate(latitude = str_replace_all(latitude, ' ', ''), longitude = str_replace_all(longitude, ' ', ''))

microreact_data <- microreact_data %>% select(id = Accession, Description, CreateDate, SequenceLength, Phylum, Class, Order, Family, Genus, Species, BiosampleAccession, year, month, day, LocationDescription_uncurated, latitude, longitude, Country, IsolationSource__autocolour = IsolationSource_maincategories, BacMetGenes, NumBacMetGenes, VFDBgenes, NumVFDBgenes, MobType, ConjType, NumIntegron, NumIn0, NumCALIN, NumIS, RepliconType, RepliconFamily, ResistanceGenes = ARGName, TotalResGenes, aminoglycoside, sulphonamide, tetracycline, phenicol, macrolide, trimethoprim, ESBL = betalactam_ESBL, carbapenem = betalactam_carbapenem, quinolone, colistin)  # problem also need to import resistance genes; try just joining pre-created microreact data

write.table(microreact_data, file = 'data/microreact-data-ncbi-plasmid-antibiotic-resistance.tsv', sep = '\t', col.names = TRUE, row.names = FALSE, na = "")


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

