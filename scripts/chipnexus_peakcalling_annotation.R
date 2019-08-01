#! /bin/env Rscript


library(Biostrings)
library(ggplot2)
library(getopt)
library(gdata)
library(idr)
library(tidyverse)


current_dir = "/home/jwen/projects/qi/elba/code_NC/"
setwd(current_dir)
data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")
out_dir = paste(data_dir, "chipnexus_peaks/",sep="")
## bam files can be downloaded from GEO
bam_dir = "/home/jwen/projects/qi/ATACseq_ChIPnexus/bam/"
genome_f =  paste(data_dir, "dm3.chrom.sizes",sep="")
genome_fasta = paste(data_dir, "dmel_noUextra.fasta",sep="") 

motif_f = paste(data_dir,"dm_TF_motif_qi.txt",sep="")

##########################################################################
########################### functions ####################################
##########################################################################

create_beds <- function(bn_dir, bn, dnp,qcutoff, nm, dofasta=F, do_fimo=F) {
	
	dnp = dnp %>% arrange(chr,chr_start, chr_end,strand)
	outf = paste(out_dir,bn_dir,"/",bn, "_",nm,".bed", sep="")
	write.table(dnp[,1:6],file=outf,row.names = F,col.names=F, quote = F,sep="\t",append = F)
	
	if (dofasta) {
		outfasta = paste(out_dir,bn_dir,"/",bn, "_",nm,".fasta", sep="")
		system(paste("fastaFromBed -fi ",genome_fasta," -name -s -bed ",outf," -fo ",outfasta,sep="" ))
		
		ffseq  = readDNAStringSet(outfasta)
		names(ffseq) = gsub("::(.*)|\\((.*)", "",names(ffseq), perl=T)
		writeXStringSet(ffseq, file=outfasta, width=20000)
	}
	
	if (do_fimo) {
		infasta_f = paste(out_dir,bn_dir,"/",bn, "_",nm,".fasta", sep="")
		system(paste("fimo --parse-genomic-coord --oc ", out_dir,bn_dir,"/",qcutoff,"/ ", motif_f, " ",infasta_f,  sep=""))
	}
}

merge_motif <- function (danno, dmotif) {
	

	mmm = readLines(motif_f )
	mmm = mmm[mmm != ""]
	mmm = mmm[regexpr("MOTIF", mmm, perl=T) !=-1]
	mmm = gsub("MOTIF ", "", mmm, perl=T)
	
	colnames(dmotif) = c("motif.pattern", "xx", "peakid", "motif.start", "motif.end", "motif.strand","motif.score","motif.pval","motif.qval","motif.matchedSequence")		
	
	dmotif = dmotif %>% 
			filter((motif.pattern %in% "insv_CCAATTGG" & motif.matchedSequence %in% c("CCTATTGG","CCAATTAG","CTAATTGG","CCAATTGG")) | (motif.pattern %in% "insv_CCAATAAG" & motif.matchedSequence %in% c("CCGATAAG","CCAATAAG")) | !(motif.pattern %in% c("insv_CCAATTGG","insv_CCAATAAG")))  %>% 
	        mutate(motif = gsub("_(.*)","",motif.pattern)) %>% dplyr::select(-xx)
	
	
	dmotif2 = data.frame(peakid = unique(dmotif$peakid), stringsAsFactors=F)
	for (oo in mmm) {					
		hh = dmotif  %>% filter(motif.pattern== oo) %>% select(peakid,motif.pattern ) %>% distinct()
		colnames(hh) = c("peakid", oo)
		dmotif2 = dmotif2 %>% left_join(hh) 
	}
	
	dmotif2[is.na(dmotif2$insv_CCAATTGG), "insv_CCAATTGG"] = "" 
	dmotif2[is.na(dmotif2$insv_CCAATAAG), "insv_CCAATAAG"] = "" 
	
	dmotif2 = mutate(dmotif2, CP190="", GAGA="")
	dmotif2[(!is.na(dmotif2$CP190_NTGGCMACACTR)), "CP190"] = "CP190" 
	dmotif2[(!is.na(dmotif2$GAGA_AGAGAGMGAGAG)), "GAGA"] = "GAGA" 
	
	dmotif = dmotif %>% left_join(dmotif2 )
	
	dmotif 
	
}

merge_motif_bak <- function (dnp,dmotif) {

	colnames(dmotif) = c("motif.pattern", "xx", "peakid", "motif.start", "motif.end", "motif.strand","motif.score","motif.pval","motif.qval","motif.matchedSequence")		

	dmotif = dmotif %>% 
			filter((motif.pattern %in% "insv_CCAATTGG" & motif.matchedSequence %in% c("CCTATTGG","CCAATTAG","CTAATTGG","CCAATTGG")) | (motif.pattern %in% "insv_CCAATAAG" & motif.matchedSequence %in% c("CCGATAAG","CCAATAAG")) | !(motif.pattern %in% c("insv_CCAATTGG","insv_CCAATAAG")))
	
	dmotif2 = data.frame(peakid = unique(dnp$PeakID), stringsAsFactors=F)
	
	for (oo in mmm) {					
		hh = dmotif  %>% filter(motif.pattern== oo) %>% select(peakid,motif.pattern ) %>% distinct()
		colnames(hh) = c("peakid", oo)
		dmotif2 = dmotif2 %>% left_join(hh) 
	}
	
	dmotif2[is.na(dmotif2$insv_CCAATTGG), "insv_CCAATTGG"] = "" 
	dmotif2[is.na(dmotif2$insv_CCAATAAG), "insv_CCAATAAG"] = "" 
	
	dmotif2 = mutate(dmotif2, CP190="", GAGA="")
	dmotif2[(!is.na(dmotif2$CP190_NTGGCMACACTR)), "CP190"] = "CP190" 
	dmotif2[(!is.na(dmotif2$GAGA_AGAGAGMGAGAG)), "GAGA"] = "GAGA" 
	
	dmotif2
	
}

############################################################################
#############################################################################

## peak calling using macs2 ##
inbams = list.files(path=bam_dir,pattern = "(.*)ChIPnexus_sorted.bam$", full=T)
for (inbam in inbams) {
	nm = gsub("_sorted.bam","", basename(inbam), perl=T)
	system(paste("macs2 callpeak -t ", inbam, " --format=BAM  --gsize=dm  --qvalue=0.05   --cutoff-analysis --call-summits  --outdir=",out_dir,nm," --name ", nm, sep=""))
}	

## peak annotation ##
setwd(out_dir)
macs_dirs = list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
macs_dirs = macs_dirs[grepl("ChIPnexus",macs_dirs)]

for ( rr in macs_dirs ) {
	cat(rr, "\n")
	bn_dir = basename(rr)
	bn = gsub("_spikein","", bn_dir, perl=T)
	cutoff_f = paste(bn_dir,"/", bn, "_cutoff_analysis.txt",sep="")
	narrowPeak_f = paste(bn_dir,"/", bn, "_peaks.narrowPeak",sep="")
	summit_f = paste(bn_dir,"/", bn, "_summits.bed",sep="")
	peak_xls = paste(bn_dir,"/", bn, "_peaks.xls",sep="")
	
	if (file.exists(narrowPeak_f) && file.info(narrowPeak_f)$size != 0)		{
		
		for (qcutoff in c(-log10(0.05), -log10(0.01), -log10(1e-5), -log10(1e-10))) {
			qcutoff=round(qcutoff,digits=1)
			
			dcutoff = read.delim(cutoff_f, header=T, stringsAsFactors=F) %>% 
					mutate( pval=10^(-pscore), qval=10^(-qscore) )
			
			dnp = read.delim(narrowPeak_f, header=F, stringsAsFactors=F) %>% 
					`colnames<-`(c("chr", "chr_start", "chr_end","peakid","peakscore","strand","foldchange", "pscore","qscore","sumitpos2peak"))  %>%
					filter(qscore >= qcutoff)
			
			dsumit = read.delim(summit_f, header=F, stringsAsFactors=F) %>% 
					`colnames<-`(c("chr", "chr_start", "chr_end","peakid", "peakscore"))  %>%
					filter(peakscore >= qcutoff) %>% mutate(strand=".")
				
			create_beds(bn_dir, bn, dnp, qcutoff=qcutoff, nm=paste("q",qcutoff,"_narrowPeak",sep=""), dofasta=T, do_fimo=T)
			create_beds(bn_dir, bn, qcutoff=round(qcutoff,digits=1), dsumit, nm=paste("q",qcutoff,"_sumit",sep=""))	
			
	
			peakfile= paste(out_dir,bn_dir,"/",bn, "_q",qcutoff,"_narrowPeak.bed", sep="")
			peakannofile = gsub(".bed", "_annot.txt", peakfile)
			statfile = gsub(".bed", "_stat.txt", peakfile)
	        system(paste("annotatePeaks.pl ", peakfile," dm3  -genomeOntology ",out_dir,bn_dir,"/ -annStats ", statfile, " > ", peakannofile,sep=""))
			
			danno = read.delim(peakannofile, header=T, stringsAsFactors=F) 
			colnames(danno)[1] = "PeakID"
			
			motiff = paste(out_dir, bn, "/", qcutoff,"/fimo.tsv",sep="")
			dmotif  = read.delim(motiff,header=T, stringsAsFactors=F)
			dmotif = merge_motif(danno, dmotif)
		
			
			danno = danno  %>% 
					mutate(coord=paste(Chr,":",Start,"-", End,sep="")) %>% 
					mutate(Annotation=gsub(" UTR", "UTR", Annotation)) %>% 
					mutate(genomic_region=gsub("\\s+(.*)", "", Annotation)) %>% 
					rename(chipsets=Focus.Ratio.Region.Size) %>% 
					left_join( dsumit %>% dplyr::select(peakid, chr_end) %>% rename(PeakID=peakid, peaksumit=chr_end))  %>% 
					left_join( dmotif,by=c("PeakID"="peakid"))  %>% filter(!is.na(motif.pattern))  %>% 
					mutate(motif.chr.start=Start+motif.start-1, motif.chr.end=Start+motif.end-1) %>% 
					mutate(motif_chr.mid=motif.chr.start+round((motif.end-motif.start)/2)) %>% 
					mutate(motifDistToPeakSummit = motif_chr.mid-peaksumit) %>%
					mutate(uid = paste(PeakID, ".", motif.pattern, sep=""))  
			
			
			load(file=paste(data_dir,"elba__gene.rep.sig_final.RData",sep=""))
			# sigGenes_list
			gg = tolower(gsub("_ChIPnexus", "_vs_yw", bn))
			rnaseq = sigGenes_list[[gg]]
			rnaseq = rnaseq %>% mutate(elba1=1/3*(elba1_1+elba1_2+elba1_3),elba2=1/3*(elba2_1+elba2_2+elba2_3), elba3=1/3*(elba3_1+elba3_2+elba3_3),insv=1/3*(insv_1+insv_2+insv_3), yw=1/3*(yw_1+yw_2+yw_3) )
			if (gg == "elba1_vs_yw") {
				rnaseq = rnaseq %>% rename(elba1_yw.logFC = logFC, elba1_yw.adj.P.Val=adj.P.Val )  %>% 
						select(gene, yw, elba1, elba1_yw.logFC, elba1_yw.adj.P.Val)
			} else if (gg == "elba2_vs_yw") {
				rnaseq = rnaseq %>% rename(elba2_yw.logFC = logFC, elba2_yw.adj.P.Val=adj.P.Val )  %>% 
						select(gene, yw, elba2, elba2_yw.logFC, elba2_yw.adj.P.Val)
			} else if (gg == "elba3_vs_yw") {
				rnaseq = rnaseq %>% rename(elba3_yw.logFC = logFC, elba3_yw.adj.P.Val=adj.P.Val )  %>% 
						select(gene, yw, elba3, elba3_yw.logFC, elba3_yw.adj.P.Val)
			} else if (gg == "insv_vs_yw") {
				rnaseq = rnaseq %>% rename(insv_yw.logFC = logFC, insv_yw.adj.P.Val=adj.P.Val )  %>% 
						select(gene, yw, insv, insv_yw.logFC, insv_yw.adj.P.Val)
			}
			danno= left_join(danno, rnaseq, by=c("Gene.Name"="gene"))	
			
			save(dcutoff, dnp, dsumit, danno, file=paste(data_dir, bn, "_q",qcutoff, "_macs2_anno.RData",sep=""))
			
			
		}
	}
}


