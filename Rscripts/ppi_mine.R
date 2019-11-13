###############
# PPI network ##########
# from STRING database #
########################
dir.string <- "/Users/schischlikf2/datasets/STRING/"

string.actions <- read.csv(
  paste(dir.string, "9606.protein.actions.v11.0.txt", sep=""),
  sep="\t")

string.links.full <- read.csv(
  paste(dir.string, "9606.protein.links.full.v11.0.txt", sep=""),
  sep="\n")

string.links <- read.csv(
  paste(dir.string, "9606.protein.links.v11.0.tab-separated.txt", sep=""),
  sep="\t")

string.info <- read.csv(
  paste(dir.string, "9606.protein.info.v11.0_cut.txt", sep=""),
  sep="\t")
  
#################################
# Extract MYCN physical binding #
#################################
string.info[string.info$preferred_name=="MYCN", ]
#protein_external_id preferred_name protein_size
#3864 9606.ENSP00000281043           MYCN          464

mycn <- string.actions[string.actions$item_id_a=="9606.ENSP00000281043" &
                         string.actions$mode=="binding", ]

mycn.merge.name <- merge(mycn, string.info, by.x="item_id_b", by.y="protein_external_id")

write.table(mycn.merge.name, 
          file="~/datasets/SLpediatric/Analysis/ppi_MYCN_partner_genes/MYCN_partner_genes.tsv",
          sep="\t", row.names=FALSE, col.names=TRUE, quote=FALSE)

##############################################
# Extract all proteins with physical binding #
##############################################

string_binding <- string.actions[string.actions$mode=="binding", ]
string_binding <- string_binding[c("item_id_a", "item_id_b", "score")]

m.string_binding <- merge(string_binding, string.info, by.x="item_id_a", by.y="protein_external_id")
colnames(m.string_binding) <- c("ensembl_a", "ensembl_b", "score", "gene_symbol_a", "potein_size")

mm.sb <- merge(m.string_binding, string.info, by.x="ensembl_b", by.y="protein_external_id")
colnames(mm.sb) <- c("ensembl_b", "ensembl_a", "score", "gene_symbol_a", "potein_size", 
                     "gene_symbol_b", "protein_size")

ppi.string <- mm.sb[c("gene_symbol_a", "gene_symbol_b", "ensembl_a", "ensembl_b", "score")]
ppi.string <- ppi.string[order(ppi.string$ensembl_a), ]

ppi.string$gene_symbol_a <- as.character(ppi.string$gene_symbol_a)
ppi.string$gene_symbol_b <- as.character(ppi.string$gene_symbol_b)
ppi.string$ensembl_a <- as.character(ppi.string$ensembl_a)
ppi.string$ensembl_b <- as.character(ppi.string$ensembl_b)

save(ppi.string, file="~/datasets/STRING/ppi.string.RData")
      