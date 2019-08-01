#! /bin/env Rscript


library(Biostrings)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(venneuler)
library(Rsubread)
library(getopt)

current_dir = "/home/jwen/projects/qi/elba/ELBA/"
setwd(current_dir)

data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")

## bam files can be downloaded from GEO
bam_dir = "/home/jwen/projects/qi/ATACseq_ChIPnexus/bam/"

##########################################################################
########################### functions ####################################
##########################################################################

featurecounts_bed <- function (bamfiles,dother, strandSpecific=2,countMultiMappingReads=F, fraction=F, primaryOnly=TRUE,minMQS=0) {
	
	dother = dother %>% arrange(GeneID,Chr,Start,End) %>% distinct()
	
	fc_PE <- featureCounts(bamfiles, 
			annot.ext = dother,
			nthreads=5,
			minMQS=minMQS,
			countMultiMappingReads = countMultiMappingReads,
			fraction = fraction, 
			strandSpecific = strandSpecific, 
			useMetaFeatures=TRUE, 
			allowMultiOverlap=TRUE,
			largestOverlap = TRUE,
			primaryOnly=primaryOnly,
			isPairedEnd=FALSE)
	
	counts=fc_PE$counts
	colnames(counts) = gsub("(.*)\\.\\.","", colnames(counts) )
	colnames(counts) = gsub("_sorted.mapped|.unique.mapped|.bam|_sorted.bam","", colnames(counts) )
	counts = counts %>% data.frame() %>% mutate(PeakID=rownames(counts)) %>% dplyr::select(PeakID, everything())
#	counts_sum = colSums(counts)
	
#	counts_sum
	
	counts
}

plot_cdf_per_motif <- function (vv, mytitle="") {
	
	vv =  vv %>% group_by(ss)  %>% mutate(nloci = n_distinct(coord)) %>% ungroup %>% select(-nloci) %>% data.frame
	
	p <- ggplot(vv, aes(orientation.index, ecdf, colour = ss))
	pg <- p + geom_step() + 
			xlab("Orientation Index (x)") + ylab("1-cummulative peak fraction") +
			scale_colour_manual(values = c("blue","green","magenta","black", "red", "darkgray", "orange","dodgerblue2","deeppink", "yellow4","brown4", "cyan2", "blue4","yellow","darkorchid2","darkslategray")) +
			theme_bw() +
			ggtitle(mytitle) +
			theme(legend.position="right")
	print(pg)
	
}

plot_orientation_index_per_motif <- function (dtop, pdff) {
	
	dmotif = dtop %>% select(ss,coord, motifcombine, orientation.index,Peak.Score) %>% distinct()
	dplot = dmotif %>% group_by(ss,motifcombine) %>% mutate(ecdf = 1-ecdf(orientation.index)(orientation.index)) %>% data.frame()
	
	dmotif %>% group_by(ss, motifcombine) %>% summarise(n())
	
	pdf(file=pdff, height=5, width=8)
	for (ii in c("insv_CCAATTGG","insv_CCAATAAG")){
		dplot_i = dplot %>% filter(motifcombine== ii)
		plot_cdf_per_motif(dplot_i, mytitle=ii)
	}
	dev.off()
	system(paste("pdfjam ",pdff,"  --frame false --nup 1x2 --suffix 2up --outfile ",pdff,sep=""))
}

########################################################################
##### Calculate ChIP-nexus orientation index ##################
########################################################################

bamfiles = list.files(path=bam_dir,pattern = "_ChIPnexus_sorted.bam$", full=T)
#Rfiles = list.files(path=data_dir,pattern = paste("_ChIPnexus_q1.3_macs2_anno.RData$",sep=""),full.names=T)
Rfiles = list.files(path=data_dir,pattern = paste("_ChIPnexus_q2_macs2_anno.RData$",sep=""),full.names=T)

Dcounts = NULL
for (Rfile in Rfiles) {
	bn = gsub("_q2_macs2_anno.RData","",basename(Rfile))
	load(Rfile)

	dother = danno %>% dplyr::select(PeakID,Chr,Start,End,Strand) %>% distinct() %>% rename(GeneID=PeakID)
	
	bnn = gsub("_ChIPnexus_q2_macs2_anno.RData", "", basename(Rfile))
	bamfile = bamfiles[grepl(bnn, bamfiles)]
	counts_fw = featurecounts_bed(bamfile, dother=dother, strandSpecific=1, countMultiMappingReads=FALSE)
	counts_rv = featurecounts_bed(bamfile, dother=dother, strandSpecific=2, countMultiMappingReads=FALSE)
	
	counts_fw_rv = counts_fw %>% full_join(counts_rv,by=c("PeakID"="PeakID")) %>%
			`colnames<-`(c("PeakID","fw","rv"))  %>%
			mutate(orientation.index=pmax(fw,rv)/(fw+rv))  %>% mutate(ss = bnn) %>% 
			left_join(dplyr::select(danno, PeakID,coord, Gene.Name, genomic_region,Distance.to.TSS, Peak.Score,peaksumit,motifDistToPeakSummit, insv_CCAATTGG,insv_TCCAATTGGA, insv_CCAATAAG,insv_TCCAATAAGA, CP190,GAGA))
	
	Dcounts = rbind(Dcounts,counts_fw_rv)
}

Dcounts[is.na(Dcounts$insv_TCCAATTGGA),"insv_TCCAATTGGA"] = ""
Dcounts[is.na(Dcounts$insv_TCCAATAAGA),"insv_TCCAATAAGA"] = ""

Dcounts = Dcounts %>% mutate(motif="")
Dcounts[!Dcounts$insv_TCCAATTGGA %in% "" & !Dcounts$insv_TCCAATAAGA %in% "" & Dcounts$motif %in% "", "motif"] = "insv_TCCAATTGGA_TCCAATAAGA"
Dcounts[!Dcounts$insv_CCAATTGG %in% "" & !Dcounts$insv_CCAATAAG %in% "" & Dcounts$motif %in% "", "motif"] = "insv_CCAATTGG_CCAATAAG"
Dcounts[!Dcounts$insv_TCCAATTGGA %in% "" & Dcounts$motif %in% "", "motif"] = "insv_TCCAATTGGA"
Dcounts[!Dcounts$insv_CCAATTGG %in% "" & Dcounts$motif %in% "", "motif"] = "insv_CCAATTGG"
Dcounts[!Dcounts$insv_TCCAATAAGA %in% "" & Dcounts$motif %in% "", "motif"] = "insv_TCCAATAAGA"
Dcounts[!Dcounts$insv_CCAATAAG %in% "" & Dcounts$motif %in% "", "motif"] = "insv_CCAATAAG"
Dcounts[Dcounts$motif %in% "", "motif"] = "others"

Dcounts = Dcounts %>% mutate(motifcombine="")
Dcounts[(!Dcounts$insv_CCAATTGG %in% "" | !Dcounts$insv_TCCAATTGGA %in% "") & (!Dcounts$insv_CCAATAAG %in% "" | !Dcounts$insv_TCCAATAAGA %in% "") & Dcounts$motifcombine %in% "", "motifcombine"] = "insv_CCAATTGG_CCAATAAG"
Dcounts[(!Dcounts$insv_CCAATTGG %in% "" | !Dcounts$insv_TCCAATTGGA %in% "") & Dcounts$motifcombine %in% "", "motifcombine"] = "insv_CCAATTGG"
Dcounts[(!Dcounts$insv_CCAATAAG %in% "" | !Dcounts$insv_TCCAATAAGA %in% "") & Dcounts$motifcombine %in% "", "motifcombine"] = "insv_CCAATAAG"
Dcounts[Dcounts$motifcombine %in% "", "motifcombine"] = "others"


table(Dcounts$motif)
table(Dcounts$motifcombine)

save(Dcounts, file=paste(data_dir, "ChIPnexus_orientation_index.RData",sep=""))


########################################################################
##### figure 3.A motif coverage chipseq vs chip-nexus ##################
########################################################################




########################################################################
##### figure 3.D strand orientation index ##############################
########################################################################

load(file=paste(data_dir, "ChIPnexus_orientation_index.RData",sep=""))
# Dcounts, Dcounts_merge

###top 500 motifs from chipnexus for each factor
topn = 500
dtop = Dcounts %>% dplyr::select(ss,coord, motifcombine,Peak.Score,  orientation.index) %>% distinct() %>% 
		group_by(ss,coord,  motifcombine, orientation.index) %>% summarise(Peak.Score = max(Peak.Score))  %>% data.frame() %>%
		arrange(-Peak.Score) %>% group_by(ss,motifcombine)  %>%  top_n(topn, Peak.Score) %>% data.frame() #Peak.Score
dtop %>% group_by(ss, motifcombine) %>% summarise(n=n_distinct(coord))
plot_orientation_index_per_motif(dtop, pdff=paste(res_dir, "fig3d.ChIPNexus_OI_top",topn,"_motifs.pdf", sep=""))



