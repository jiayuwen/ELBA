#! /bin/env Rscript


library(Biostrings)
library(getopt)
library(gdata)
library(ggplot2)
library(tidyverse)
library(readr)
library(limma)
library(Rsubread)
library(edgeR)
library(openxlsx)
library(GenomicFeatures)
library(GenomicRanges)


current_dir = "/home/jwen/projects/qi/elba/ELBA/"
setwd(current_dir)

data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")
peak_dir = paste(data_dir, "proseq_peaks/",sep="")
genome_f =  paste(data_dir, "dm3.chrom.sizes",sep="")

############################################################################################################
##################################### function #############################################################

get_neighborGeneDist <- function (dpeak,neighborGene, genes, direction) {
	#  dpeak = dpeak; neighborGene=peak_down; genes=dpro_starts
	
	dpeak$neighborGeneId = neighborGene
	dpeak  = dpeak[!is.na(dpeak$neighborGeneId)]
	neighborGeneIds = dpeak$neighborGeneId
	
	dist2up = distance(dpeak, genes[neighborGeneIds], ignore.strand=TRUE)
	dpeak$neighborGeneStrand = strand(genes[neighborGeneIds])
	dpeak$neighborGeneId = genes[neighborGeneIds]$id
	dpeak$neighborGeneStart = start(genes[neighborGeneIds]) 
#	dpeak$neighborGeneEnd = end(genes[neighborGeneIds]) 
	
	dpeak$neighborGeneDistance = dist2up
	dpeak = dpeak %>% data.frame() 
	dpeak$mergedPeakId = as.character(dpeak$mergedPeakId); dpeak$strand = as.character(dpeak$strand); dpeak$neighborGeneStrand=as.character(dpeak$neighborGeneStrand)
	dpeak$neighborGeneId = as.character(dpeak$neighborGeneId)
	dpeak2 = dpeak %>% mutate(coord=paste(seqnames,":", start ,  sep=""))
	
	
	load(file=paste(data_dir, "combined_pro__ProseqE_tssFold4_bodyFold3_annot_round2.RData",sep=""))
# dpeak,dPI, dPI_promoter ,dPI_coding, dPI_noncoding, 
	dPI = dPI %>% dplyr::filter(grepl("AlignToDmOnly",genotype)) %>% mutate(genotype=gsub("_sorted.AlignToDmOnly", "", genotype)) %>% 
			dplyr::select(PeakID,genotype, body_log2TPM, promoter_log2TPM,PI, PI.ntile, Gene.Name,Gene.Strand, Distance.to.TSS)
	colnames(dPI)[2:ncol(dPI)] = paste(direction,"streamGene.",colnames(dPI)[2:ncol(dPI)],sep="")
	
	dpeak2 = dpeak2 %>% left_join(dPI, by=c("neighborGeneId"="PeakID"))
	colnames(dpeak2) = gsub("neighborGene",paste(direction,"streamGene.",sep=""),colnames(dpeak2))
	
	dpeak2 = dpeak2 %>% dplyr::select(-strand,-gid) 
	
	dpeak2
	
}

get_paircate <- function (g_xx) {
	
	g_xx[g_xx$upstreamGene.Strand %in% "+" & g_xx$downstreamGene.Strand %in% "+", "paircate"] = "tandem"
	g_xx[g_xx$upstreamGene.Strand %in% "+" & g_xx$downstreamGene.Strand %in% "-", "paircate"] = "convergent"
	g_xx[g_xx$upstreamGene.Strand %in% "-" & g_xx$downstreamGene.Strand %in% "+", "paircate"] = "divergent"
	g_xx[g_xx$upstreamGene.Strand %in% "-" & g_xx$downstreamGene.Strand %in% "-", "paircate"] = "tandem"
	
	g_xx
}

dist2peak <- function (dd, inbed, tarbed) {
	
	#inbed =probedf ; tarbed = peakbedf 
	exprf = tempfile(pattern = paste("tmp_",Sys.getpid(),sep=""), tmpdir = tmp_dir, fileext ="txt") 
	
	system(paste('closestBed -D a -t all -g ', genome_f, ' -a ',inbed,' -b ',tarbed,'  > ',exprf,sep=''))
#	system(paste('closestBed -D a -io -t all  -a ',inbed,' -b ',inbed,'  > ',exprf,sep=''))
	
	if (file.info(exprf)$size == 0) {
		outi = NULL
	} else {
		outi = read.delim(file=exprf,header=F,stringsAsFactors=F)
	}
	
	outi = subset(outi, select=c(V4,V10, V13))
	colnames(outi) = c("peakid", "nearestPeak", "dist2NearestPeak")
	
	
	dt<-as.data.table(outi)
	setkey(dt,peakid)
	xx1 = dt[,list(nearestPeak=paste(nearestPeak, collapse=";",sep=""), dist2NearestPeak= paste(dist2NearestPeak, collapse=";",sep="")),by=peakid]
	outi  =data.frame(xx1)			
	
	unlink(exprf)
	outi
	
}

############################################################################################################


load(file=paste(data_dir, "combined_pro__ProseqE_tssFold4_bodyFold3_annot_round2.RData",sep=""))
# dpeak,dPI, dPI_promoter ,dPI_coding, dPI_noncoding, 
dPI_high = dPI %>% filter(promoter_log2TPM>=1)
probedf = paste(data_dir, "proseq_peaks/combined_pro_tssFold4_bodyFold3_round2_anchor.bed",sep="")
dpro = read.delim(probedf,header=F, stringsAsFactors=F)
colnames(dpro) = c("chr","start", "end", "id","score","strand")
dpro = dpro %>% dplyr::select(chr,start,end,strand,id)  %>% filter(id %in% dPI_high$PeakID) %>% 
		arrange(chr,start,end,id) %>% distinct() %>%  mutate(gid=1:n())
dpro_starts = makeGRangesFromDataFrame(dpro, keep.extra.columns=TRUE)


dss = read.delim(paste(data_dir,"/merge_chipseq_chipnexus_q10/strigent_merge2chip_q10_given.bed",sep=""),stringsAsFactors=F, header=F)
peakfile  = paste(data_dir,"/merge_chipseq_chipnexus_q10/merge2chip_q10_given.txt",sep="")
dm = read.delim(peakfile,stringsAsFactors=F)
colnames(dm)[1] = c("mergedPeakId")
colnames(dm) = gsub("_narrowPeak_q10.bed|_narrowPeak_q5.bed","", colnames(dm) )
dmm = dm %>% dplyr::select(chr,start,end,mergedPeakId) %>% 
		arrange(chr,start,end,mergedPeakId) %>% distinct() %>%  mutate(gid=1:n())
dpeak = makeGRangesFromDataFrame(dmm, keep.extra.columns=TRUE)


peak_down = precede(dpeak, dpro_starts, ignore.strand=TRUE)
dist_peak_down = get_neighborGeneDist(dpeak = dpeak, neighborGene=peak_down, genes=dpro_starts, direction="down")

peak_up = follow(dpeak, dpro_starts, ignore.strand=TRUE)
dist_peak_up = get_neighborGeneDist(neighborGene=peak_up, dpeak = dpeak, genes=dpro_starts, direction="up")

dist_peak = dist_peak_up %>% full_join(dist_peak_down)
dist_peak  = get_paircate(dist_peak )
dist_peak = dist_peak %>% mutate(pro.promoter_absFC=abs(upstreamGene.promoter_log2TPM-downstreamGene.promoter_log2TPM)) %>% 
		mutate(pro.body_absFC=abs(upstreamGene.body_log2TPM-downstreamGene.body_log2TPM)) %>%
		mutate(pro.PI_absFC=abs(upstreamGene.PI-downstreamGene.PI)) %>%
		dplyr::filter(upstreamGene.genotype == downstreamGene.genotype)

save(dist_peak, file=paste(data_dir, "merge_chipseq_chipnexus_q10__pro_round2.RData",sep=""))

