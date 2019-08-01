#! /bin/env Rscript


library(Biostrings)
library(ggplot2)
library(getopt)
library(gdata)
library(idr)
library(tidyverse)
library(qvalue)
library(limma)
library(edgeR)

current_dir = "/home/jwen/projects/qi/elba/ELBA/"
setwd(current_dir)
data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")

#####################################################################################################
######################################### functions #################################################

do_limma_p_rep <-function (ddx, libsize, savef ) {
	
	# ddx = genelevel_counts
	
	libsize_sel = libsize[libsize$lib %in% colnames(ddx), ]
	ddx = ddx[, libsize_sel$lib]
	dge <- DGEList(counts=ddx, lib.size=libsize_sel$mapped_reads)
	dim(dge)

	
	dge <- calcNormFactors(dge)
	

	Treat <- factor(libsize_sel$genotype)
	design <- model.matrix(~0+Treat)
	colnames(design) <- levels(Treat)
	

	if (!do_quality_weights) {
		v <- voom(dge,design,plot=FALSE)
	} else {
		v <- voomWithQualityWeights(dge,design,plot=FALSE)
	}
	
	
	if (do_weights) {
		fit <- lmFit(v,design)
	} 
	
	
	cont.matrix <- makeContrasts(
			elba1_vs_yw = elba1 - yw,
			elba2_vs_yw = elba2 - yw,
			elba3_vs_yw = elba3 - yw,
			insv_vs_yw = insv - yw,
			levels=design)
	mypairs = c("elba1_vs_yw", "elba2_vs_yw", "elba3_vs_yw","insv_vs_yw")
	
	
	
	fit2  <- contrasts.fit(fit, cont.matrix)
	fit2  <- eBayes(fit2)
	colnames(fit2)
	
	
	
	res = list()
	for (pp in mypairs) {
		res[[pp]] = topTreat(fit2, coef=pp, number=nrow(dge$counts))
	}
   
	save(res,v, design, cont.matrix,fit2,  file=savef)
}

output_table <- function (res, v, savef) {
	
	sigGenes_list = list()
	
	vv = data.frame(v$E)
	vv = mutate(vv, gene=rownames(vv))
	
	
	for (mm in names(res)) {
		sigGenes = res[[mm]]
		sigGenes = sigGenes %>% mutate( gene=rownames(sigGenes), absLogFC = abs(logFC)) %>% select(gene,logFC,absLogFC, P.Value,adj.P.Val)
		sigGenes = sigGenes %>% inner_join(vv, by="gene") %>% select(-matches("logFC|P.Value|adj.P.Val"), matches("logFC|P.Value|adj.P.Val"))

		sigGenes_list[[mm]] = sigGenes
	}
	save(sigGenes_list,vv, file=gsub(".RData$", "_final.RData", savef) )
	
}

#######################################################################
#######################################################################


load(file=paste(data_dir, "elba_RNAseq_new.RData", sep=""))
#libsize
libsize = libsize %>% mutate(genotype=gsub("\\_(\\d+)","", lib)) %>% arrange(genotype)


load(file=paste(data_dir, "dm3_RNAseq_gene_elba_counts.RData",sep=""))
# mycounts_merge

genelevel_counts = mycounts_merge
rownames(genelevel_counts) = genelevel_counts$Gene
genelevel_counts = genelevel_counts[,-1]
genelevel_counts = genelevel_counts[rowSums(genelevel_counts) >=50,]

############## do limma ###############################################
do_quality_weights = FALSE 
do_splicing = TRUE 
do_robust_fit = FALSE
do_weights = TRUE

savef = paste(data_dir, "elba__gene.rep.sig.RData", sep="")
do_limma_p_rep(ddx = genelevel_counts, libsize,savef)

load(file=savef)
output_table(res,v,savef)



