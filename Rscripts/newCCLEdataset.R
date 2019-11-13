# Create new CCLE dataset
# Author Fiorella Schischlik, some functions were provided by Joo Sang Lee

library(Rcpp)
library(parallel)

dir.code <- "~/datasets/SLpediatric/R/"
dir.rdata <- "~/datasets/SLpediatric/Rdata_object/"
load(paste(dir.code, "prob.avana.extended.RData", sep=""))
load(paste(dir.code, "prob.TCGA.extended.qnorm.RData", sep=""));prob0=prob

dir <- "/Users/schischlikf2/datasets/DepMap/CCLE_2019_02_09/"

# achilles eq avana
achilles <- read.csv(paste(dir, "Achilles_gene_effect.csv", sep=""),
                     sep=",", header=TRUE, row.names = 1)
copynumber <- read.csv(paste(dir, "CCLE_gene_cn.csv", sep=""),
                     sep=",", header=TRUE, row.names = 1)
expression <- read.csv(paste(dir, "CCLE_expression.csv", sep=""),
                       sep=",", header=TRUE, row.names = 1)
mutations <- read.csv(paste(dir, "CCLE_mutations.csv", sep=""),
                       sep=",", header=TRUE)
meta <- read.csv(paste(dir, "sample_info.csv", sep=""),
                      sep=",", header=TRUE, row.names = 1)
meta$DepMap_ID <- row.names(meta)

# only take overlapping genes between 
# achilles, copynumber, expression
genes.achilles <- colnames(achilles)
genes.copynumber <- colnames(copynumber)
genes.expression <- colnames(expression)

intersect.genes.ach.cn <- intersect(genes.achilles, genes.copynumber)
intersect.genes.all <- intersect(union.genes.ach.cn, genes.expression)

achilles.g <- achilles[intersect.genes.all]
copynumber.g <- copynumber[intersect.genes.all]
expression.g <- expression[intersect.genes.all]

# only take cell lines with overlapping copynumber
# and expression information

cell.lines.achilles <- row.names(achilles)
cell.lines.copynumber <- row.names(copynumber)
cell.lines.expression <- row.names(expression)

intersect.cl.ach.cn <- intersect(cell.lines.achilles, cell.lines.copynumber)
intersect.cell.lines.all <- intersect(intersect.cl.ach.cn, cell.lines.expression)

# order rownmames
achilles.g <- achilles.g[match(intersect.cell.lines.all, rownames(achilles.g)),]
copynumber.g <- copynumber.g[match(intersect.cell.lines.all, rownames(copynumber.g)),]
expression.g <- expression.g[match(intersect.cell.lines.all, rownames(expression.g)),]

# replace gene names 
gene.symbols <- sapply(strsplit(colnames(achilles.g), split="\\.\\."), head, 1)
colnames(achilles.g) <- gene.symbols
colnames(copynumber.g) <- gene.symbols
colnames(expression.g) <- gene.symbols

# replace column names with CCLE.Name
m.cl <- meta[row.names(achilles.g),]
m.cl <- m.cl[match(row.names(achilles.g), row.names(m.cl)),]

row.names(achilles.g) <- m.cl$CCLE.Name
row.names(copynumber.g) <- m.cl$CCLE.Name
row.names(expression.g) <- m.cl$CCLE.Name

# transpose all datasets
achilles.g <- t(achilles.g)
copynumber.g <- t(copynumber.g)
expression.g <- t(expression.g)

# cell lines & types
cell.lines <- as.character(m.cl$stripped_cell_line_name)

prob.avana.09.19 <- list(
  genes=row.names(achilles.g),
  samples=colnames(achilles.g),
  mat=achilles.g,
  scna=copynumber.g,
  celllines=as.character(m.cl$stripped_cell_line_name),
  types=as.character(toupper(m.cl$sample_collection_site)), 
  mRNA=expression.g
)

# preprocessing of data...
preprocess <- function(mod.prob, prob, micro.rna=FALSE){
  
  #m mod.prob --> prob for modification
  
  if(micro.rna){
    ix=match(prob$genes, mod.prob$genes)
    mod.prob$genes=mod.prob$genes[ix]
    mod.prob$mat=mod.prob$mat[ix,]
    mod.prob$mRNA=mod.prob$mRNA[ix,]
    mod.prob$scna=mod.prob$scna[ix,]
    
    sorted <- sort(mod.prob$genes)
    id = match(sorted, mod.prob$genes)
    mod.prob$genes=mod.prob$genes[id]
    mod.prob$mat=mod.prob$mat[id,]
    mod.prob$mRNA=mod.prob$mRNA[id,]
    mod.prob$scna=mod.prob$scna[id,]
  }
  
  else{
    ix=match(prob$genes, mod.prob$genes)
    mod.prob$genes=mod.prob$genes[ix]
    mod.prob$mat=mod.prob$mat[ix,]
    mod.prob$mRNA=mod.prob$mRNA[ix,]
    mod.prob$scna=mod.prob$scna[ix,]
  }
  
  mRNA.norm = mclapply(1:length(mod.prob$samples), 
                       function(tt) qnorm.array(mod.prob$mRNA[,tt]),mc.cores=5)
  mRNA.norm = do.call(cbind, mRNA.norm)
  mRNA.norm = mclapply(1:length(mod.prob$genes), 
                       function(tt) qnorm.array(mRNA.norm[tt,]),mc.cores=5)
  mod.prob$mRNA.norm = t(do.call(cbind, mRNA.norm))
  scna.norm = mclapply(1:length(mod.prob$samples), 
                       function(tt) qnorm.array(mod.prob$scna[,tt]),mc.cores=5)
  scna.norm = do.call(cbind, scna.norm)
  scna.norm = mclapply(1:length(mod.prob$genes), 
                       function(tt) qnorm.array(scna.norm[tt,]),mc.cores=5)
  mod.prob$scna.norm = t(do.call(cbind, scna.norm))
  
  mRNA.rank = apply(mod.prob$mRNA, 1, rank.array) 
  mod.prob$mRNAq2 = t(1*(mRNA.rank > 1/3) + 1*(mRNA.rank > 2/3))
  scna.rank = apply(mod.prob$scna, 1, rank.array) 
  mod.prob$scnaq2 = t(1*(scna.rank > 1/3) + 1*(scna.rank > 2/3))
  return(mod.prob)
}

rank.array <- function(mat){
  mat.back = mat 
  mat = mat[!is.na(mat)]
  mat = rank(mat, ties.method = "average")/length(mat);
  mat.back[!is.na(mat.back)] = mat
  mat.back
}

qnorm.array <- function(mat){
  mat.back = mat 
  mat = mat[!is.na(mat)]
  mat = rank(mat, ties.method = "average");
  mat = qnorm(mat / (length(mat)+1));
  mat.back[!is.na(mat.back)] = mat 
  mat.back
}

prob.mod.avana.09.19 <- preprocess(prob.avana.09.19, prob0, micro.rna=FALSE)

save(prob.mod.avana.09.19, 
     file=paste(dir.rdata, "prob.mod.avana.09.19.Rdata", sep=""))

# select meta & mutations overlapping with prob.mod.avana.09.19
ccle.mutations.09.19 <- mutations[mutations$DepMap_ID %in% row.names(m.cl), ]
ccle.meta.09.19 <- m.cl

write.table(ccle.mutations.09.19, 
            file=paste(dir.code, "ccle.mutations.09.19.txt", sep=""),
            quote=FALSE, 
            sep="\t", 
            row.names=FALSE)

write.table(ccle.meta.09.19, 
            file=paste(dir.code, "ccle.meta.09.19.txt", sep=""),
            quote=FALSE, 
            sep="\t", 
            row.names=FALSE)









