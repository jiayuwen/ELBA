#! /bin/env Rscript


library(Biostrings)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(venneuler)
library(getopt)

current_dir = "/home/jwen/projects/qi/elba/code_NC/"


data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")

genome_f =  paste(data_dir, "dm3.chrom.sizes",sep="")
genome_fasta = paste(data_dir, "dmel_noUextra.fasta",sep="") 

motif_f = paste(data_dir, "dm_TF_motif_qi.txt",sep="")
mmm = readLines(motif_f )
mmm = mmm[mmm != ""]
mmm = mmm[regexpr("MOTIF", mmm, perl=T) !=-1]
mmm = gsub("MOTIF ", "", mmm, perl=T)

############################################################################################################
##############################################    function           #######################################

mergefiles <- function (dist, lnfiles, nm) {
	
	lnfiles_str = paste(sort(lnfiles,decreasing=T), sep="", collapse=" ")
	outf = paste(nm, "_",dist,".txt", sep="")
	vennf = paste(nm, "_",dist,"_venn.txt", sep="")
	vennxls = paste(nm, "_",dist,"_venn.xls", sep="")
	
	system(paste("mergePeaks -d ", dist, " ", lnfiles_str ," -venn ",vennf, " > ", outf,sep=""))

}

create_beds <- function(bn_dir, dnp, nm, dist, dofasta=F) {
	
	
	outf = paste(bn_dir,"/", nm,"_",dist,".bed", sep="")
	unlink(outf)
	write.table(dnp[,1:6],file=outf,row.names = F,col.names=F, quote = F,sep="\t",append = F)
	
	outf_bb = paste(bn_dir,"/", nm,"_",dist,".bb", sep="")
	system(paste("bedToBigBed ",outf," ",genome_f," ",outf_bb,sep=""))
	
	if (dofasta) {
		outfasta = paste(bn_dir,"/", nm,"_",dist,".fasta", sep="")
		system(paste("fastaFromBed -fi ",genome_fasta," -name -s -bed ",outf," -fo ",outfasta,sep="" ))
		
		dseq = readDNAStringSet(outfasta)
		dseq = toupper(dseq)
		dseq2 = DNAStringSet(dseq)
		names(dseq2) = gsub("::(.*)", "", names(dseq2), perl=T)
		writeXStringSet(dseq2, outfasta, format="fasta", width=20000)
	}
	
}

call_fimo <- function (merge_dir, infile, motif_f, nm) {
	
	outdir = paste(merge_dir, sep="")
	
	system(paste("fimo --parse-genomic-coord --oc ", outdir," ", motif_f, " ",infile,  sep=""))

}

merge_motif <- function (dd, dmotif) {
	
	colnames(dmotif) = c("motif.pattern", "xx", "peakid", "motif.start", "motif.end", "motif.strand","motif.score","motif.pval","motif.qval","motif.matchedSequence")		
	dmotif = dmotif %>% mutate(peakid = gsub("\\((.*)", "", peakid))
	
	tt = dmotif %>% filter((motif.pattern %in% "insv_CCAATAAG"))
	table(tt$motif.matchedSequence)
	
	tt = dmotif %>% filter((motif.pattern %in% "insv_CCAATTGG"))
	table(tt$motif.matchedSequence)
	
	dmotif = dmotif %>% 
			filter((motif.pattern %in% "insv_CCAATTGG" & motif.matchedSequence %in% c("CCTATTGG","CCAATTAG","CTAATTGG","CCAATTGG")) | (motif.pattern %in% "insv_CCAATAAG" & motif.matchedSequence %in% c("CCGATAAG","CCAATAAG")) | !(motif.pattern %in% c("insv_CCAATTGG","insv_CCAATAAG")))
	
	dmotif2 = data.frame(peakid = unique(dd$PeakID), stringsAsFactors=F)
	for (oo in mmm) {
		
		hh = dmotif  %>% filter(motif.pattern== oo) %>% select(peakid,motif.pattern ) %>% distinct()
		colnames(hh) = c("peakid", oo)
		dmotif2 = dmotif2 %>% left_join(hh) 
#						
	}
	
	dmotif2[is.na(dmotif2$insv_CCAATTGG), "insv_CCAATTGG"] = "" 
	dmotif2[is.na(dmotif2$insv_CCAATAAG), "insv_CCAATAAG"] = "" 
	
	dmotif2 = mutate(dmotif2, CP190="", GAGA="")
	dmotif2[(!is.na(dmotif2$CP190_NTGGCMACACTR)), "CP190"] = "CP190" 
	dmotif2[(!is.na(dmotif2$GAGA_AGAGAGMGAGAG)), "GAGA"] = "GAGA" 
	
	dmotif2
	
}

get_maxScore <- function (merge_dir,dd) {
	
	dm = dd 
	dscore = dd %>% select(PeakID) %>% data.frame()
	sss =c("Elba1_ChIPnexus", "Elba2_ChIPnexus", "Elba3_ChIPnexus", "Insv_ChIPnexus","wtElba1__elba1Elba1", "wtElba2__elba2Elba2", "wtElba3__elba3Elba3","wtInsv__insvInsv")
	for (bn2 in sss) {
		cat(bn2,"\n")
		dm2_i = dm[, regexpr(paste("PeakID|",bn2,sep="",collapse=""), colnames(dm),perl=T) !=-1]
		colnames(dm2_i) = c("mergedPeakId", "peakid")
		npeaks = sapply(str_extract_all(paste(dm2_i$peakid,",",sep=""), ","), length)
		dm2_i = mutate(dm2_i, npeaks=npeaks)
		dm2_i = data.frame(mergedPeakId=dm2_i[rep(1:nrow(dm2_i), dm2_i$npeaks),"mergedPeakId"], peakid = unlist(strsplit(paste(dm2_i$peakid,",",sep=""),split=",")),stringsAsFactors=F )
		
		rfile = paste(merge_dir,"/", bn2,"_macs2.RData",sep="")
		load(rfile)
		dm2_i = dm2_i %>% inner_join(select(dnp, peakid, peakscore), by=c("peakid"="peakid"))  %>% distinct()   %>%
				group_by(mergedPeakId) %>% summarise(maxPeakscore=max(peakscore,na.rm = TRUE))
		dscore = dscore %>% left_join(dm2_i, by=c("PeakID" = "mergedPeakId"))
	}
	
	colnames(dscore) = c("PeakID", paste(sss,".maxScore",sep=""))
	dscore[is.na(dscore)] = 0
	dscore = dscore  %>%  
			mutate(ChIPnexus.maxScore=apply(cbind(Elba1_ChIPnexus.maxScore, Elba2_ChIPnexus.maxScore,Elba3_ChIPnexus.maxScore, Insv_ChIPnexus.maxScore), 1,  max, na.rm = TRUE)) %>%
			mutate(ChIPseq.maxScore=apply(cbind(wtElba1__elba1Elba1.maxScore,wtElba2__elba2Elba2.maxScore,wtElba3__elba3Elba3.maxScore,wtInsv__insvInsv.maxScore), 1,  max, na.rm = TRUE))
	dscore = dscore  %>%  mutate(maxScore=apply(cbind(ChIPnexus.maxScore,ChIPseq.maxScore), 1,  max, na.rm = TRUE)) %>% 
			mutate(Elba1.maxScore=apply(cbind(Elba1_ChIPnexus.maxScore,wtElba1__elba1Elba1.maxScore), 1,  max, na.rm = TRUE))  %>% 
			mutate(Elba2.maxScore=apply(cbind(Elba2_ChIPnexus.maxScore,wtElba2__elba2Elba2.maxScore), 1,  max, na.rm = TRUE))  %>% 
			mutate(Elba3.maxScore=apply(cbind(Elba3_ChIPnexus.maxScore,wtElba3__elba3Elba3.maxScore), 1,  max, na.rm = TRUE))  %>% 
			mutate(Insv.maxScore=apply(cbind(Insv_ChIPnexus.maxScore,wtInsv__insvInsv.maxScore), 1,  max, na.rm = TRUE))  
	dscore
}	

merge2chip_anno <- function (merge_dir, nm) {
	
	mergef = paste(merge_dir,"/", nm,"_given.txt",  sep="")
	dd0 = read.delim(mergef, stringsAsFactors=F, header=T)
	colnames(dd0)[1] = "mergeid"
	colnames(dd0) =gsub("_narrowPeak_q(.*).bed", "",colnames(dd0) )
	
	dl=list()
	cc = colnames(dd0)[(which(colnames(dd0) == "Total.subpeaks")+1):ncol(dd0)]
	for (ii in cc) {
		dl[[ii]] = unique(dd0[dd0[,ii] != "","mergeid"])
	}
	lapply(dl, length)
	dbed = dd0 %>% dplyr::select( chr, start,end, mergeid, Total.subpeaks, strand) %>% mutate(start = start-1) %>% arrange( chr, start,end) 
	create_beds(bn_dir=merge_dir, dbed, nm=nm,dist="given", dofasta=T)
	
	for (bn in c("Elba1","Elba2","Elba3","Insv")) {
		dbed_ii = dd0 %>% filter_at(vars(contains(bn)), all_vars(. !="")) %>%
				dplyr::select( chr, start,end, mergeid, Total.subpeaks, strand) %>% 
				mutate(start = start-1) %>% arrange( chr, start,end) 
		create_beds(bn_dir=merge_dir, dnp=dbed_ii, nm=paste(bn, "_",nm,sep=""),dist="given", dofasta=F)
	}
	
	
	myset = paste(nm,"_given",sep="")
	dsumit = dbed %>% mutate(peak=start+round((end-start)/2,digits=0)) %>% mutate(start=peak-1, end=peak)%>%select(-peak)
	bedf=paste(merge_dir,"/",myset,"_sumit.bed",sep="")
	write.table(dsumit,file=bedf,row.names = F,col.names=F, quote = F,sep="\t",append = F)
	
	
	for (flank in c(100, 200, 500)) {
		bedf_flank = paste(merge_dir,"/",nm,"_given_sumit_flank", flank,".bed",sep="")
		fastaf_flank =  paste(merge_dir,"/merge2chip_given_sumit_flank", flank,".fasta",sep="")
		system(paste("slopBed -i ",bedf," -g ",genome_f," -s -l ",flank," -r ",flank," > ",bedf_flank,sep=""))
		system(paste("fastaFromBed -fi ",genome_fasta," -name -s -bed ",bedf_flank," -fo ",fastaf_flank,sep="" ))
		
		dseq = readDNAStringSet(fastaf_flank)
		dseq = toupper(dseq)
		dseq2 = DNAStringSet(dseq)
		names(dseq2) = gsub("::(.*)", "", names(dseq2), perl=T)
		writeXStringSet(dseq2, fastaf_flank, format="fasta", width=20000)
	}
	
	infasta_f = paste(merge_dir,"/",myset,".fasta", sep="")
	call_fimo(merge_dir, infile=infasta_f, motif_f, nm=nm)
	
	peakfile = paste(merge_dir,"/", myset, ".txt", sep="")
	peakannofile = gsub(".txt", "_annot.txt", peakfile, perl=T)
	statfile = gsub(".txt", "_stat.txt", peakfile, perl=T)
	system(paste("annotatePeaks.pl ", peakfile," dm3  -go ", merge_dir,"/GO/", myset,"/ -genomeOntology ",merge_dir,"/GenomeOntology/", myset,"/ -annStats ", statfile, " > ", peakannofile,sep=""))
	
	
	dd = read.delim(peakannofile, stringsAsFactors=F, header=T)
	colnames(dd)[1] = "PeakID"
	dd = dd %>% 
			mutate(coord=paste(Chr,":",Start,"-", End,sep="")) %>% 
			mutate(Annotation=gsub(" UTR", "UTR", Annotation)) %>% 
			mutate(genomic_region=gsub("\\s+(.*)", "", Annotation)) %>% 
			rename(chipsets=Focus.Ratio.Region.Size) %>% 
			left_join(select(dd0, mergeid, matches("ChIPnexus|wtElba|wtInsv")), by=c("PeakID"="mergeid"))
	
	dscore = get_maxScore(merge_dir,dd)
	
	dd = dd  %>% 
			left_join(dscore) %>% arrange(desc( maxScore))  %>%
			select(PeakID, coord, Gene.Name, genomic_region, Elba1_ChIPnexus.maxScore, Elba2_ChIPnexus.maxScore, maxScore, chipsets, Distance.to.TSS,Annotation, Detailed.Annotation, Nearest.PromoterID, everything()) 
	
	
	motiff = paste(merge_dir,"/fimo.tsv",sep="")
	dmotif  = read.delim(motiff,header=T, stringsAsFactors=F)
	dmotif2 = merge_motif(dd, dmotif )
	dd = left_join(dd, dmotif2, by=c("PeakID"="peakid"))
	
	
	load(file=paste(data_dir,"elba__gene.rep.sig_final.RData",sep=""))
	for (gg in names(sigGenes_list)) {
		
		rnaseq = sigGenes_list[[gg]]
		rnaseq = rnaseq %>% mutate(elba1=1/3*(elba1_1+elba1_2+elba1_3),elba2=1/3*(elba2_1+elba2_2+elba2_3), elba3=1/3*(elba3_1+elba3_2+elba3_3),insv=1/3*(insv_1+insv_2+insv_3), yw=1/3*(yw_1+yw_2+yw_3) )
		if (gg == "elba1_vs_yw") {
			rnaseq = rnaseq %>% rename(elba1_yw.logFC = logFC, elba1_yw.adj.P.Val=adj.P.Val )  %>% 
					select(gene, yw, elba1, elba1_yw.logFC, elba1_yw.adj.P.Val)
		} else if (gg == "elba2_vs_yw") {
			rnaseq = rnaseq %>% rename(elba2_yw.logFC = logFC, elba2_yw.adj.P.Val=adj.P.Val )  %>% 
					select(gene, elba2, elba2_yw.logFC, elba2_yw.adj.P.Val)
		} else if (gg == "elba3_vs_yw") {
			rnaseq = rnaseq %>% rename(elba3_yw.logFC = logFC, elba3_yw.adj.P.Val=adj.P.Val )  %>% 
					select(gene, elba3, elba3_yw.logFC, elba3_yw.adj.P.Val)
		} else if (gg == "insv_vs_yw") {
			rnaseq = rnaseq %>% rename(insv_yw.logFC = logFC, insv_yw.adj.P.Val=adj.P.Val )  %>% 
					select(gene, insv, insv_yw.logFC, insv_yw.adj.P.Val)
		}
		dd = left_join(dd, rnaseq, by=c("Gene.Name"="gene"))	
		
	}
	
	save(dd, dsumit, file=paste(data_dir, myset,  "_anno.RData",sep=""))
}


############################################################################################################
################################################# merge diff peak ##########################################
setwd(data_dir)
merge_dir = "merge_chipseq_chipnexus_q10"
dir.create(merge_dir)
setwd(merge_dir)
lnfiles = list.files(path=".",pattern = "(.*)", full=F)
lnfiles = lnfiles[grepl("ChIPnexus|wtElba|wtInsv",lnfiles) & !grepl("sumit|merge2chip|RData",lnfiles) ]
mergefiles(lnfiles, dist="given", nm="merge2chip_q10")
setwd(data_dir)
merge2chip_anno(merge_dir, nm="merge2chip_q10")














