# Library -----------------------------------------------------------------
library(survcomp)
library(mppa)
library(RobustRankAggreg)
library(tidyverse)
library(data.table)
library(org.Hs.eg.db)
library(tidyverse)


# Note --------------------------------------------------------------------
batch = "2023-12-21" #to change for different version of data
#download gtex
#GTEx v8 specific gene expression		https://fuma.ctglab.nl/tutorial/download_variants
#download gtex_v8_ts_avg_log2TPM.txt from the above link and save it to data folder
#it has been downloaded into the data folder in the github repo

#save impc_mouse_human_homlogs.RData to data/ folder. it has been downloaded into the data folder in the github repo
#gwas.out.RData this file is larger than 100 mb. Download it from box https://thejacksonlaboratory.ent.box.com/folder/164757838622
# and save to data/ folder.

# KOMP AT results ---------------------------------------------------------
meta.data.file <- paste0("data/KOMP2018_Meta_Table_v1.2.csv")
meta.data <- read.csv(meta.data.file, header=TRUE)
meta.data <- meta.data[meta.data$Domain != "EKG",]
pheno.domain.list <- unique(meta.data$Domain)
#read
AT.summary.stat.path <- paste0("output/", batch, "_AT/")
all.AT.summary.stat <- NULL
for(pd in pheno.domain.list){
  AT.summary.stat.file <- paste0(AT.summary.stat.path, pd, "_stats_LASSO_", batch, "_wing1months.csv")
  AT.summary.stat <- as.data.frame(fread(AT.summary.stat.file, header=TRUE, sep=","))
  all.AT.summary.stat <- rbind(all.AT.summary.stat, AT.summary.stat)
}
dim(all.AT.summary.stat)
pheno.list <- as.character(meta.data$Phenotype)
all.AT.summary.stat <- subset(all.AT.summary.stat, pheno %in% pheno.list)
all.AT.summary.stat$geno.fdr <- p.adjust(all.AT.summary.stat$geno.pval, method="BH", n=length(na.omit(all.AT.summary.stat$geno.pval)))
## add short domain and phenotype names
dim(all.AT.summary.stat)
all.AT.summary.stat <- merge(droplevels(all.AT.summary.stat), droplevels(meta.data), by.x=c("domain","pheno"), by.y=c("Domain","Phenotype"), all.x=TRUE)
dim(all.AT.summary.stat)
# add short gene name
all.AT.summary.stat$mutant.short <- str_replace(all.AT.summary.stat$mutant, "n-TAagc16", "n_TAagc16")
all.AT.summary.stat$mutant.short[all.AT.summary.stat$mutant.short == "C57BL/6J-Tg(tetO-NEK1)5Lutzy/J"] = "Tg_tetO_NEK1"
all.AT.summary.stat$mutant.short <- gsub(".*-\\s*|<.*", "",all.AT.summary.stat$mutant.short)


# Ranking 1 ---------------------------------------------------------------
#Ranking 1: Combining the pvalues of the behavior phenotype using simes.test. Only use the behavior phenotype pvalue.
beha_gene_pheno <- all.AT.summary.stat[all.AT.summary.stat$Category == "Behavior", ]
beha_gene_pheno_split <- split(beha_gene_pheno, beha_gene_pheno$mutant.short)
beha_gene_pheno_simes <- lapply(beha_gene_pheno_split, function(x){
  simes.test(na.omit(x$geno.pval))
})
beha_gene_pheno_simesP <- data.frame(gene = names(unlist(beha_gene_pheno_simes)),
                                     simes.p = unlist(beha_gene_pheno_simes))
# Ranking 2 ---------------------------------------------------------------
#Ranking 2: Rank genes for psychiatric phewas. In a similar way as 1, combine pvalues for each geneusing PsyGenenetics GWAS phenotypes.
#phewas results
load("data/impc_mouse_human_homlogs.RData")
load("data/gwas.out.RData")
phewas_gene <- mouse.human.homlogs[mouse.human.homlogs$mouseGene %in% beha_gene_pheno_simesP$gene,]
phewas_gene_simes <- gwas[phewas_gene[phewas_gene$humanGene %in% names(gwas),"humanGene"]]
phewas_gene_simes <- lapply(phewas_gene_simes, function(x){
  simes.test(x[x$Domain == "Psychiatric",]$P.value)
})
phewas_gene_simesP <- data.frame(gene = names(unlist(phewas_gene_simes)),
                                 simes.p = unlist(phewas_gene_simes))
phewas_gene_simesP <- phewas_gene_simesP %>%
  dplyr::arrange(simes.p) %>%
  dplyr::distinct(gene, .keep_all = TRUE)

# Ranking 3 ---------------------------------------------------------------
#3,Rank genes based on brain specific expression – use Gtex or another resource (does not have to be mouse) to determine how brain enriched each gene is.
gtex_gene_log2TPM <- read.table("data/gtex_v8_ts_avg_log2TPM.txt", header = TRUE, sep = "\t")
annots <- AnnotationDbi::select(org.Hs.eg.db, keys=gtex_gene_log2TPM$GENE,
                 columns="SYMBOL", keytype="ENSEMBL") %>%
  as.data.frame() %>%
  dplyr::distinct(ENSEMBL, .keep_all = T)
gtex_gene_log2TPM <- dplyr::left_join(gtex_gene_log2TPM, annots, by = c("GENE" = "ENSEMBL")) %>%
  dplyr::relocate(SYMBOL, .after = GENE)
gtex_gene_log2TPM$ratio1 <- gtex_gene_log2TPM$ratio2 <- NULL
for(i in 1:nrow(gtex_gene_log2TPM)){
  gtex_gene_log2TPM$ratio1[i] <- mean(as.numeric(gtex_gene_log2TPM[i,c(10,12,15,16,17,19,20,22)]))/mean(as.numeric(gtex_gene_log2TPM[i,3:56]))
  gtex_gene_log2TPM$ratio2[i] <- mean(as.numeric(gtex_gene_log2TPM[i,c(10,12,15,16,17,19,20,22)]))/mean(as.numeric(gtex_gene_log2TPM[i,10:22]))
}
#merge gene
beha_genes_merge <- merge(beha_gene_pheno_simesP,
                          phewas_gene,
                          by.x = "gene",
                          by.y = "mouseGene",
                          all.x = TRUE)[,c(-4,-5)]
beha_genes_merge <- merge(beha_genes_merge,
                          phewas_gene_simesP,
                          by.x = "humanGene",
                          by.y = "gene",
                          all.x = TRUE)
beha_genes_merge <- merge(beha_genes_merge,
                          gtex_gene_log2TPM[,c(1,2,58)],
                          by.x = "humanGene",
                          by.y = "SYMBOL",
                          all.x = TRUE)[,-5]
beha_genes_merge_nona <- na.omit(beha_genes_merge)
colnames(beha_genes_merge_nona)[3:5] <- c("simes.p.beha",
                                          "simes.p.psy",
                                          "ratio.gtex")

rank_beha_genes_merge_nona <- data.frame(matrix(0, nrow = nrow(beha_genes_merge_nona),ncol = 3))
rownames(rank_beha_genes_merge_nona) <- beha_genes_merge_nona$gene
rank_beha_genes_merge_nona[,1] <- rank(log10(beha_genes_merge_nona$simes.p.beha))
rank_beha_genes_merge_nona[,2] <- rank(log10(beha_genes_merge_nona$simes.p.psy))
rank_beha_genes_merge_nona[,3] <- rank(-beha_genes_merge_nona$ratio.gtex)
colnames(rank_beha_genes_merge_nona) <- c("rank.simes.beha",
                                          "rank.simes.psy",
                                          "rank.ratio.gtex")
glist <- list(rownames(rank_beha_genes_merge_nona[order(rank_beha_genes_merge_nona[,1]),]),
              rownames(rank_beha_genes_merge_nona[order(rank_beha_genes_merge_nona[,2]),]),
              rownames(rank_beha_genes_merge_nona[order(rank_beha_genes_merge_nona[,3]),]))

# Aggregate the inputs by RRA
aggrg_rank_RRA <- aggregateRanks(glist = glist, N = length(glist[[1]]), method = "RRA")

rank_beha_genes_merge_nona_RRA <- merge(aggrg_rank_RRA,
                                    rank_beha_genes_merge_nona,
                                    by.x = "Name",
                                    by.y = "row.names",
                                    all.x = TRUE)
rank_beha_genes_merge_nona_RRA <- rank_beha_genes_merge_nona_RRA[order(rank_beha_genes_merge_nona_RRA$Score),]
write.csv(rank_beha_genes_merge_nona_RRA, file = "output/rank_beha_genes_merge_nona_RRA.csv", quote = F, row.names = F)
