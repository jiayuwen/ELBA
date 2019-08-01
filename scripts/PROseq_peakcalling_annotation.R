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


current_dir = "/home/jwen/projects/qi/elba/ELBA/"
setwd(current_dir)

data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")
peak_dir = paste(data_dir, "proseq_peaks/",sep="")

genome_f =  paste(data_dir, "dm3.chrom.sizes",sep="")

############################################################################################################

get_anchor_pos <- function (ddx, mypos,outfile="") {
	colnames(ddx) = c("V1","V2", "V3","V4","V5","V6")
	
	ddx_plus = ddx[ddx$V6 == "+",]
	ddx_minus = ddx[ddx$V6 == "-",]
	
	if (mypos == "start") {
		ddx_plus$V3 = ddx_plus$V2
		ddx_minus$V2 = ddx_minus$V3
	} else if (mypos == "end") {
		ddx_plus$V2 = ddx_plus$V3
		ddx_minus$V3 = ddx_minus$V2
	}
	
	
	ddx =  rbind(ddx_plus,ddx_minus )
	ddx = ddx[with(ddx, order(V1, V2, V3)), ]
	ddx = ddx[complete.cases(ddx),]
	ddx = ddx %>% arrange(V1, V2, V3)
	
	if (outfile != "") {
		unlink(outfile)
		write.table(ddx,file=outfile,row.names=F,col.names=F,quote = F,sep="\t")
	}
	ddx
}

make_expr <- function (bedf, nm) {
	dd = read.table(file=bedf, header=F, stringsAsFactors=F)
	if (ncol(dd) < 6) {
		dd = rbind(mutate(dd, V6="+"), mutate(dd, V6="-"))
		dd = dd %>% select(V4, V1, V2, V3, V6)  %>% mutate(nm=nm)
		
	} else {
		dd = dd %>% select(V4, V1, V2, V3, V6)  %>% mutate(nm=nm)
	}
	colnames(dd) = c("GeneID","Chr","Start","End","Strand", "nm")
	dd
}

fpkmToTpm <- function(fpkm, dolog=T)
{
#	exp(log(fpkm) - log(sum(fpkm)) + log(1e6))
	#	tpm =  log2(fpkm) - log2(sum(fpkm)) + log2(1e6)
	#	tpm = (fpkm / sum(fpkm)) * 10^6		
	
	J = t(data.frame(colSums(fpkm)))
	J = J[rep(1,nrow(fpkm)), ]
	
	tpm = log2(fpkm) - log2(J) + log2(1e6)			
	if (!dolog) {
		tpm = 2^tpm
	}
	tpm	
}

do_annotation_bed <- function (bamfiles,dother, savef,strandSpecific=2,countMultiMappingReads=F, fraction=F, primaryOnly=TRUE,minMQS=0, region) {
	
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
	
	
	
	x_PE <- DGEList(counts=fc_PE$counts, genes=fc_PE$annotation) 
	x_PE_rpkm = rpkm(x_PE,x_PE$genes$Length, log=T)
	x_PE_tpm = fpkmToTpm(2^x_PE_rpkm, dolog=T)	
	colnames(x_PE_tpm) = gsub("_q20.bam|_sorted.bam|_primary.bam|.swapstrand.bam|(.*)P10817\\.","", colnames(x_PE_tpm), perl=T)
	x_tpm = x_PE_tpm
	
	
	x_tpm = data.frame(x_tpm)
	x_tpm = cbind(fc_PE$annotation,x_tpm )
	
	save(fc_PE,  dother, x_tpm, file=savef)
	
}


######################################################
######## peak calling parameters #####################

# findPeaks -style groseq -tssFold 4 -bodyFold 3

######################################################
######## peaks annotation ############################

bn = "combined_pro"
dpeak = read.delim(file=paste(peak_dir,  bn,"_tssFold4_bodyFold3_annot.txt",sep=""), header=T, stringsAsFactors=F)
colnames(dpeak)[1] = "PeakID" 
colnames(dpeak)[20] = "normalization.factor" 
colnames(dpeak)[21] = "Distance.to.nearest.Peak..Peak.ID" 
dpeak = dpeak %>% mutate(coord=paste(Chr,":",Start,"-",End,sep="")) %>% mutate(Annotation.category=gsub("\\s+(.*)", "", Annotation)) %>% 
		dplyr::select(PeakID,Gene.Name,coord, Chr,Start,End,Strand,Peak.Score,Focus.Ratio.Region.Size,Annotation.category, Annotation,Detailed.Annotation,Distance.to.TSS,Nearest.PromoterID,Gene.Type,Distance.to.nearest.Peak..Peak.ID)
    
############ produce bed for promoter/body TPM ############
outbed = paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_round2.bed",sep="")
dpeak_bed = dpeak %>% dplyr::select(Chr,Start,End,PeakID, Peak.Score,Strand ) %>% mutate(Peak.Score=round(Peak.Score,digits=0)) %>% mutate(Start=Start-1) %>% arrange(Chr,Start,End)  %>% distinct()
write.table(dpeak_bed,file=outbed,row.names = F,col.names=F, quote = F,sep="\t",append = F)	

outbed_startAnchor = paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_round2_anchor.bed",sep="")
xx = get_anchor_pos(ddx=dpeak_bed , mypos="start", outfile = outbed_startAnchor)
pflank = 250 
outbed_promoter = paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_promoter", pflank,"_round2.bed",sep="")
system(paste("slopBed -i ",outbed_startAnchor," -g ", genome_f," -s -l 0 -r ",pflank," | sort  -k1,1 -k2,2n > ",outbed_promoter,sep=""))

pp = read.delim(file=outbed_promoter, header=F, stringsAsFactors=F)
colnames(pp) = c("Chr", "pStart",   "pEnd", "PeakID", "Peak.Score", "Strand")
pp2 = pp %>% inner_join(dplyr::select( dpeak_bed, PeakID, Start, End))
pp_plus = pp2 %>% filter(Strand == "+")  %>% mutate(blen=End-pEnd) %>% 
		mutate(End = ifelse(blen <200, End+200, End)) %>% 
		mutate(bStart=Start+400, bEnd=End) %>%
		mutate(blen=bEnd-bStart) %>% filter(blen >10)  %>% 
        select(Chr, bStart, bEnd, PeakID, Peak.Score, Strand)  %>% rename(Start=bStart, End=bEnd)
pp_minus = pp2 %>% filter(Strand == "-")  %>% mutate(blen=pStart-Start) %>% 
		mutate(Start = ifelse(blen <200, Start-200, Start)) %>% 
		mutate(bStart=Start, bEnd=pEnd-400) %>%
		mutate(blen=bEnd-bStart) %>% filter(blen >10)  %>% 
		select(Chr,  bStart, bEnd,  PeakID, Peak.Score, Strand)  %>% rename(Start=bStart, End=bEnd)
BB = pp_plus %>% bind_rows(pp_minus) %>%  mutate(PeakID= paste(PeakID, ".body",sep="")) %>% arrange(Chr,Start, End )
PP = pp %>% rename(Start=pStart, End=pEnd) %>% mutate(PeakID= paste(PeakID, ".promoter",sep="")) %>% arrange(Chr,Start, End )
PPBB = PP %>% bind_rows(BB) %>% arrange(Chr,Start, End )
write.table(PPBB,file=paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_promoter_body_", pflank,"_round2.bed",sep=""),row.names = F,col.names=F, quote = F,sep="\t",append = F)	


system(paste("annotatePeaks.pl ",outbed_promoter," dm3 -d ", peak_dir,"/  -p ", outbed_promoter," -strand + > ",peak_dir,"/",bn,"_tssFold4_bodyFold3_promoter_annot_round2.txt", sep=""))
dpeak = read.delim(file=paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_promoter_annot_round2.txt",sep=""), header=T, stringsAsFactors=F)
colnames(dpeak)[1] = "PeakID" 
colnames(dpeak)[20] = "totalReads" 
colnames(dpeak)[21] = "Distance.to.nearest.Peak..Peak.ID" 
dpeak = dpeak %>% mutate(coord=paste(Chr,":",Start,"-",End,sep="")) %>% mutate(Annotation.category=gsub("\\s+(.*)", "", Annotation)) %>% 
		select(PeakID,Gene.Name,coord, Chr,Start,End,Strand,Peak.Score,Focus.Ratio.Region.Size,Annotation.category, Annotation,Detailed.Annotation,Distance.to.TSS,Nearest.PromoterID,Gene.Type,Distance.to.nearest.Peak..Peak.ID)


outbed_body = paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_body_TSSdownstream400_round2.bed",sep="")
write.table(pp_plus %>% bind_rows(pp_minus) %>% arrange(Chr,Start, End ),file=outbed_body,row.names = F,col.names=F, quote = F,sep="\t",append = F)	
system(paste("annotatePeaks.pl ",outbed_body," dm3 -d ", peak_dir,"/  -strand + > ",peak_dir,"/",bn,"_tssFold4_bodyFold3_body_annot_round2.txt", sep=""))
dpeak_body = read.delim(file=paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_body_annot_round2.txt",sep=""), header=T, stringsAsFactors=F)
colnames(dpeak_body)[1] = "PeakID" 
colnames(dpeak_body)[20] = "totReads" 
dpeak_body = dpeak_body %>% mutate(coord=paste(Chr,":",Start,"-",End,sep="")) %>% mutate(Annotation.category=gsub("\\s+(.*)", "", Annotation)) %>% 
		select(PeakID,Gene.Name,coord, Chr,Start,End,Strand,Peak.Score,Focus.Ratio.Region.Size,Annotation.category, Annotation,Detailed.Annotation,Distance.to.TSS,Nearest.PromoterID,Gene.Type,totReads)

save(dpeak,dpeak_body, file=paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_promoter_annot_round2.RData",sep="") )

############################## calculate TPM ###################

bamfiles = list.files(path=peak_dir,pattern = ".swapstrand.bam$", full=TRUE)
ddx = make_expr(bedf = paste(peak_dir, "/",bn,"_tssFold4_bodyFold3_promoter_body_",pflank,"_round2.bed",sep=""), nm = bn)
ddx = ddx %>% filter(End-Start>0)

do_annotation_bed(bamfiles, dother=ddx, savef=paste(peak_dir, bn, "__ProseqE_tssFold4_bodyFold3_promoter_body_",pflank,"_round2.RData",sep=""), strandSpecific=1, countMultiMappingReads=TRUE, fraction=FALSE, primaryOnly=TRUE, minMQS=0)


############ Pausing Index ######################################

dgene = read.delim(paste(data_dir, "dm6flybase_gene.bed",sep=""),stringsAsFactors=F, header=F)
dgene = dgene %>% dplyr::select(V4, V6) %>% rename(Gene.Name=V4, Gene.Strand=V6)

load(paste(peak_dir, "/", bn,"_tssFold4_bodyFold3_promoter_annot_round2.RData",sep=""))
#dpeak
load(paste(peak_dir, bn, "__ProseqE_tssFold4_bodyFold3_promoter_body_",pflank,"_round2.RData",sep=""))

dPI = x_tpm %>% mutate(PeakID=gsub("(.*)\\.(body|promoter)", "\\1", GeneID,perl=T), region=gsub("(.*)\\.(body|promoter)", "\\2", GeneID,perl=T)) %>% 
		 select(-GeneID, -Chr,-Start,-End,-Strand,-Length) %>% 
		 gather(genotype, tpm, -c(PeakID,region)) %>%  spread(region, tpm) %>%
         rename(body_log2TPM=body, promoter_log2TPM=promoter) %>%
		 mutate(PI=2^promoter_log2TPM/(2^promoter_log2TPM+2^body_log2TPM)) %>% 
		 mutate(PI.ntile= ntile( PI,4)) 

dPI =  dPI %>% inner_join(dpeak) %>% left_join(dgene) %>% arrange(desc(PI))

dPI_coding = dPI %>% filter(Gene.Type %in% "protein-coding")  %>% filter(Gene.Strand == Strand)
dPI_promoter = dPI_coding %>% filter(abs(Distance.to.TSS) <= 1000 | Annotation.category %in% "promoter-TSS") %>% mutate(PI.promoter.ntile= ntile( PI,4)) 
dPI_coding = dPI_coding %>% filter(!PeakID %in% dPI_promoter$PeakID)
dPI_noncoding =dPI %>% filter(!Gene.Type %in% "protein-coding")
 
save(dpeak,dPI, dPI_promoter ,dPI_coding, dPI_noncoding, file=paste(data_dir, bn,"__ProseqE_tssFold4_bodyFold3_annot_round2.RData",sep="") )


