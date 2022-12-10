pacman::p_load(tidyr, tidyverse, readxl, openxlsx, stringr, igraph, RColorBrewer)
dir.create('output_descriptive_analysis', showWarnings = FALSE)

# load data
accessions <- read_excel("data/Data_S1.xlsx", sheet = 'J')
restyping <- read_excel("data/Data_S1.xlsx", sheet = 'I')
plasmiddf <- read_excel("data/Data_S1.xlsx", sheet = 'G')
plasmiddf <- plasmiddf %>% mutate(across(everything(), ~ str_replace_all(.x, "^\'|\'$", "")))
colnames(plasmiddf) <- colnames(plasmiddf) %>% str_replace_all("^\'|\'$", "")

# ---------------------------
# AMINOGLYCOSIDE/SULPHONAMIDE CO-OCCURRENCE
#restyping %>% filter(str_detect(ARGlabel, ';')) %>% filter(str_detect(ARGlabel, 'sul1') & str_detect(ARGlabel, 'aadA2b')) %>% select(ARGlabel) %>% nrow()  # 90 edges between aadA2b and sul1

restyping_ami_sul <- restyping %>% filter(str_detect(ARGType, 'aminoglycoside') | str_detect(ARGType, 'sulphonamide'))

restyping_list <- list()
for (i in 1:nrow(restyping_ami_sul)) {
  Accession <- restyping_ami_sul[i, 'Accession', drop = TRUE]
  ARGlabels <- restyping_ami_sul[i, 'ARGlabel', drop = TRUE]
  ARGlabels <- sort(unique(unlist(str_split(ARGlabels, '; '))))
  if (length(ARGlabels) > 1) {
    ARGlabels <- as.data.frame(t(as.data.frame(combn(ARGlabels, m = 2), stringsAsFactors = FALSE)))  # get all pairwise combinations of ARGs per plasmid
    ARGlabels$Accession <- Accession
    restyping_list[[Accession]] <- ARGlabels
  }
}

restyping_ami_sul_edges <- do.call(rbind, restyping_list)
#restyping_ami_sul_edges %>% mutate(edge = str_c(V1, V2)) %>% filter(str_detect(edge, 'sul1') & str_detect(edge, 'aadA2b')) %>% distinct(Accession) %>% nrow()
#restyping_ami_sul_edges %>% mutate(edge = str_c(V1, V2)) %>% count(edge) %>% arrange(desc(n))

# exclude same ARG type pairs
restyping_ami_sul_edges <- restyping_ami_sul_edges %>% filter((str_detect(V1, 'aminoglycoside') & str_detect(V2, 'sulphonamide')) | (str_detect(V1, 'sulphonamide') & str_detect(V2, 'aminoglycoside'))) %>% mutate(edge = str_c(V1, V2, sep = '--'))


# plot graph
set.seed(42)
nodes <- unique(c(as.character(restyping_ami_sul_edges$V1), as.character(restyping_ami_sul_edges$V2))) %>% as.data.frame()
colnames(nodes) <- 'id'
nodes <- nodes %>% separate(col = id, into = c('id', 'type'), sep = '\\|') %>% arrange(id)
nodes <- nodes %>% mutate(type_recoded = type) %>% mutate(type_recoded = recode(type_recoded, 'aminoglycoside_quinolone' = 'aminoglycoside'))
nodes$color <- c("#FBB4AE", "#B3CDE3")[as.numeric(as.factor(nodes$type_recoded))]

links <- restyping_ami_sul_edges %>% count(edge) %>% arrange(desc(n)) %>% separate(col = edge, into = c('from', 'to'), sep = '--') %>% rename('weight' = n) %>% mutate(from = str_replace(from, '\\|.+', ''), to = str_replace(to, '\\|.+', ''))

g <- graph_from_data_frame(links,vertices = nodes,directed = F)
#plot(g)
g_pruned <- delete_edges(g, which(E(g)$weight < 50))
g_pruned <- delete_vertices(g_pruned, which(degree(g_pruned) < 1))
#plot(g_pruned)
E(g_pruned)$width <- E(g_pruned)$weight/50
layout <- layout.fruchterman.reingold(g_pruned)

pdf('output_descriptive_analysis/ami_sul_network.pdf')
plot(g_pruned, vertex.label.color = 'black', vertex.frame.color = "light grey", layout = layout)
legend("bottomright", legend = c("aminoglycoside", "sulphonamide"), col = c("#FBB4AE", "#B3CDE3"), pch = 19, pt.cex = 1.5, bty = "n", title = "Antibiotic resistance type")
dev.off()



# count edges and link to host taxonomy info from accessions data
restyping_ami_sul_count_all <- restyping_ami_sul_edges %>% count(edge) %>% arrange(desc(n)) %>% separate(col = edge, into = c('gene1', 'gene2'), sep = '--')
restyping_ami_sul_count_host_taxonomy <- restyping_ami_sul_edges %>% left_join((accessions %>% select(Accession, HostTaxonomy)), by = 'Accession') %>% group_by(HostTaxonomy) %>% count(edge) %>% arrange(desc(n), .by_group = TRUE) %>% separate(col = edge, into = c('gene1', 'gene2'), sep = '--')

restyping_ami_sul_count_host_taxonomy_split <- split(restyping_ami_sul_count_host_taxonomy, restyping_ami_sul_count_host_taxonomy$HostTaxonomy)
restyping_ami_sul_count_Enterobacteriaceae  <- restyping_ami_sul_count_host_taxonomy_split$Enterobacteriaceae
restyping_ami_sul_count_Proteobacteria_non_Enterobacteriaceae <- restyping_ami_sul_count_host_taxonomy_split$`Proteobacteria_non-Enterobacteriaceae`
restyping_ami_sul_count_other <- restyping_ami_sul_count_host_taxonomy_split$other

# save
write.xlsx(list("all_taxa" = restyping_ami_sul_count_all, "Enterobacteriaceae" = restyping_ami_sul_count_Enterobacteriaceae, "Proteobacteria_nonEnterobac" = restyping_ami_sul_count_Proteobacteria_non_Enterobacteriaceae, "restyping_ami_sul_count_other" = restyping_ami_sul_count_other), file = "output_descriptive_analysis/ami_sul_cooccurrence.xlsx")


# ---------------------------
# QUINOLONE RESISTACE IN SMALL PLASMIDS
# filter for small quinolone plasmids - what quinolone genes, and what rep types?

restyping_qnl <- restyping %>% filter(str_detect(ARGType, 'quinolone')) %>% select(Accession, ARGName, ARGType, ARGlabel)

# link to replicon typing info from plasmiddf
restyping_qnl <- restyping_qnl %>% left_join((plasmiddf %>% transmute(Accession, RepliconType, RepliconFamily, PlasmidSize = as.numeric(SequenceLength))), by = 'Accession')
restyping_qnl_small <- restyping_qnl %>% filter(PlasmidSize < 10000)

ARGlabels_all <- vector()
for (i in 1:nrow(restyping_qnl_small)) {
  ARGlabels <- restyping_qnl_small[i, 'ARGlabel', drop = TRUE]
  ARGlabels <- unique(unlist(str_split(ARGlabels, '; ')))
  for (ARGlabel in ARGlabels) {
    if (str_detect(ARGlabel, 'quinolone')) {
      ARGlabels_all <- c(ARGlabels_all, ARGlabel)
    }
  }
}

restyping_qnl_small_count <- table(ARGlabels_all) %>% as.data.frame() %>% separate(col = ARGlabels_all, into = c('ARGName','ARGClass'), sep = '\\|') %>% arrange(desc(Freq))

qnrD1_reptypes <- restyping_qnl_small %>% filter(str_detect(ARGlabel, 'qnrD1')) %>% count(RepliconType) %>% arrange(desc(n))
qnrS2_reptypes <- restyping_qnl_small %>% filter(str_detect(ARGlabel, 'qnrS2')) %>% count(RepliconType) %>% arrange(desc(n))
qnrB19_reptypes <- restyping_qnl_small %>% filter(str_detect(ARGlabel, 'qnrB19')) %>% count(RepliconType) %>% arrange(desc(n))


# save
write.xlsx(list("qnl_ARG_count" = restyping_qnl_small_count, "qnrD1_reptypes" = qnrD1_reptypes, "qnrS2_reptypes" = qnrS2_reptypes, "qnrB19_reptypes" = qnrB19_reptypes), file = "output_descriptive_analysis/qnl_small_plasmids.xlsx")




#################### OLD

hist(log10(restyping_qnl$PlasmidSize))
test <- plasmiddf %>% mutate(PlasmidSize = as.numeric(SequenceLength)) %>% filter(PlasmidSize < 500000)
hist(log10(test$PlasmidSize))


restyping_qnl_small %>% count(RepliconType) %>% arrange(desc(n))
restyping_qnl_small %>% count(ARGName) %>% arrange(desc(n))  # need to split apart and tally individual genes (max 1 per plasmid)
# need to ask what plasmids carry the major qnr genes. First ask what are the most common quinolone ARGs. For top quinolone ARGs, what are the most common replicon type haplotpes and types.




##################### DEBUGGING
# problem KP975077.1 has 3 edges that are the same
ARGlabels <- restyping %>% as.data.frame() %>% filter(Accession == 'KP975077.1') %>% select(ARGlabel)
#ARGlabels <- sort(unique(str_split(ARGlabels, '; ', simplify = TRUE)))  # fail
is.vector(unlist(str_split(ARGlabels, '; ')))








############################### OLD CODE
# separate collapsed column data into rows
restyping <- restyping %>% separate_rows(ARGProbe, ARGName, ARGClass, ARGType, ARGlabel, sep = '; ')

# exclude duplicate ARGs per plasmid
#restyping %>% nrow()  # 27237
#restyping %>% group_by(Accession, ARGlabel) %>% slice_head(n = 1) %>% nrow()  # 26304
restyping <- restyping %>% group_by(Accession, ARGlabel) %>% slice_head(n = 1) %>% ungroup()

# ---------------------------
# ANALYSE CO-OCCURRENCE OF AMINOGLYCOSIDE/SULPHONAMIDE
# create edgelist from ARGlabel column
restyping_ami_sul <- restyping %>% filter(str_detect(ARGType, 'aminoglycoside') | str_detect(ARGType, 'sulphonamide'))
restyping_ami_sul <- restyping_ami_sul %>% separate_rows(ARGProbe, ARGName, ARGClass, ARGType, ARGlabel, sep = '; ')

restyping_ami_sul <- restyping_ami_sul %>% as.data.frame() %>% select(Accession, ARGlabel) %>% group_by(Accession) %>%
  filter(n()>=2) %>% group_by(Accession) %>%
  do(data.frame(t(combn(.$ARGlabel, 2)), stringsAsFactors=FALSE))

restyping_ami_sul_edges <- restyping_ami_sul %>% filter((str_detect(X1, 'aminoglycoside') & str_detect(X2, 'sulphonamide')) | (str_detect(X1, 'sulphonamide') & str_detect(X2, 'aminoglycoside'))) %>% mutate(edge = str_c(X1, X2, sep = '--'))


restyping_ami_sul_edges %>% count(edge) %>% arrange(desc(n)) %>% separate(col = edge, into = c('gene1', 'gene2'), sep = '--')

