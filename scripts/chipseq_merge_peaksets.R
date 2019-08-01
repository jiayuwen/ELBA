#! /bin/env Rscript


library(Biostrings)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(venneuler)
library(getopt)

current_dir = "/home/jwen/projects/qi/elba/code_NC/"
setwd(current_dir)

data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")
peak_dir  = paste(data_dir, "chipseq_peaks/",sep="")


############################################################################################################
##############################################    function           #######################################

create_beds <- function(bn_dir, dnp, nm, dist, dofasta=F) {
	
	
	outf = paste(bn_dir,"/", nm,"_",dist,".bed", sep="")
	unlink(outf)
	write.table(dbed[,1:6],file=outf,row.names = F,col.names=F, quote = F,sep="\t",append = F)
	
	outf_bb = paste(bn_dir,"/", nm,"_",dist,".bb", sep="")
	system(paste("bedToBigBed ",outf," ",genome_f," ",outf_bb,sep=""))
	
	if (dofasta) {
		outfasta = paste(bn_dir,"/", nm,"_",dist,".fasta", sep="")
		system(paste("fastaFromBed -fi ",genome_fasta," -name -s -bed ",outf," -fo ",outfasta,sep="" ))
		
		ffseq  = readDNAStringSet(outfasta)
		names(ffseq) = gsub("\\((.*)", "",names(ffseq), perl=T)
		writeXStringSet(ffseq, file=outfasta, width=20000)
	}
	
}

call_fimo <- function (infile, motif_f, nm) {
	
	outdir = paste(fimo_dir,  nm, sep="")
	
	system(paste("fimo --parse-genomic-coord --oc ", outdir," ", motif_f, " ",infile,  sep=""))
	
#	mast [options] <motif file> <sequence file>
	
}


mergeset <- function (indirs, dist, nm) {

	lnfiles = NULL
	for (jj in indirs) {
		cat(jj,"\n")
		infile = list.files(path=paste("../",jj,sep=""),pattern = "(.*)_narrowPeak.bed", full=T)
		if (length(infile) == 1) {
			lnfile = gsub("diff__|_q20|_narrowPeak.bed", "", basename(infile), perl=T)
			system(paste("ln -sf ", infile, " ",lnfile, sep=""))
			lnfiles = c(lnfiles, lnfile)
		}
	}
	
	lnfiles_str = paste(sort(lnfiles,decreasing=T), sep="", collapse=" ")
	
	outf = paste(nm, "_",dist,".txt", sep="")
	vennf = paste(nm, "_",dist,"_venn.txt", sep="")
	vennxls = paste(nm, "_",dist,"_venn.xls", sep="")
	
	system(paste("mergePeaks -d ", dist, " ", lnfiles_str ," -venn ",vennf, " > ", outf,sep=""))
		
}


mergeset_anno <- function (indirs, outdir,dist="given", nm, savef) {
	
	peakfile  = paste(nm, "_given.txt",sep="")
	dm = read.delim(peakfile,stringsAsFactors=F)
	colnames(dm)[1] = c("mergedPeakId")
	write.table(dm,file=peakfile,row.names = F,col.names=F, quote = F,sep="\t",append = T)
	dm = dm %>% mutate(coord=paste(chr,":",start,"-", end,sep="")) 
	
	
	peakannofile = gsub(".txt", "_annot.txt", peakfile, perl=T)
	statfile = gsub(".txt", "_stat.txt", peakfile, perl=T)
	GO_dir = paste("GO/",sep="")
	GenomeOntology_dir = paste("GenomeOntology/",sep="")
	system(paste("annotatePeaks.pl ", peakfile," dm3  -go ", GO_dir, " -genomeOntology ",GenomeOntology_dir, " -annStats ", statfile, " > ", peakannofile,sep=""))
	
	load(file = paste(data_dir,"diffpeaks_anno_q20TRUE_spikeinFALSE_correctedMotif.RData",sep=""))
	selected = intersect(names(peak_l_uni), indirs)
	
	Dsuper1 = NULL
	Dsuper2 = select(dm, mergedPeakId) %>% distinct()
	for (bn in selected) { 
		
		bn2 = gsub("_q20|diff__","",bn,perl=T)
		
		ddi = peak_l_uni[[bn]]
		ddi = ddi %>% dplyr::select(peakid, nearestTSS,peakscore, insv_CCAATTGG,insv_CCAATAAG, CP190, GAGA, BEAF32,EBox, SuH,  matches("dist2NearestTSS|TSS_up2k|five_prime_UTR|three_prime_UTR|CDS|intron")) %>% distinct()		
		ddi[is.na(ddi)] = ""
		
		dm2_i = dm[, regexpr(paste("mergedPeakId|",bn2,sep="",collapse=""), colnames(dm),perl=T) !=-1]
		colnames(dm2_i) = c("mergedPeakId", "peakid")
		npeaks = sapply(str_extract_all(paste(dm2_i$peakid,",",sep=""), ","), length)
		dm2_i = mutate(dm2_i, npeaks=npeaks)
		dm2_i = data.frame(mergedPeakId=dm2_i[rep(1:nrow(dm2_i), dm2_i$npeaks),"mergedPeakId"], peakid = unlist(strsplit(paste(dm2_i$peakid,",",sep=""),split=",")),stringsAsFactors=F )
		dm2_i = dm2_i %>% inner_join(ddi)
		
		
		dm2_i1 = dm2_i %>% dplyr::select( mergedPeakId,nearestTSS, insv_CCAATTGG,  insv_CCAATAAG, CP190,GAGA, BEAF32,EBox, SuH,  matches("dist2NearestTSS|TSS_up2k|five_prime_UTR|three_prime_UTR|CDS|intron"))
		dm2_i2 = dm2_i %>% dplyr::select(mergedPeakId,peakid, peakscore) %>% distinct()
		dm2_i2 = dm2_i2 %>% group_by(mergedPeakId) %>% summarise(peakid=paste(unique(peakid), collapse=";",sep=""), maxPeakscore=max(peakscore,na.rm = TRUE))
		colnames(dm2_i2 )[2:ncol(dm2_i2 )] = paste(bn2,".",colnames(dm2_i2 )[2:ncol(dm2_i2 )],sep="" )
		
		if (is.null (Dsuper1) ) {
			Dsuper1 = dm2_i1
		} else {	
			Dsuper1 = bind_rows(Dsuper1, dm2_i1)
		}
		Dsuper2 = Dsuper2 %>% left_join(dm2_i2)
	}
	
	Dsuper1 = Dsuper1 %>% distinct() %>% group_by(mergedPeakId) %>% summarise(nearestTSS=first(nearestTSS), 
			insv_CCAATTGG=paste(unique(insv_CCAATTGG), collapse=";",sep=""), insv_CCAATAAG=paste(unique(insv_CCAATAAG), collapse=";",sep=""), 
			CP190=paste(unique(CP190), collapse=";",sep=""), GAGA=paste(unique(GAGA), collapse=";",sep=""), 
			BEAF32=paste(unique(BEAF32), collapse=";",sep=""), EBox=paste(unique(EBox), collapse=";",sep=""),  SuH=paste(unique(SuH), collapse=";",sep=""), 
			TSS_up2k=paste(unique(TSS_up2k), collapse=";",sep=""),five_prime_UTR=paste(unique(five_prime_UTR), collapse=";",sep=""),
			three_prime_UTR=paste(unique(three_prime_UTR), collapse=";",sep=""),CDS=paste(unique(CDS), collapse=";",sep=""),
			intron=paste(unique(intron), collapse=";",sep="")) 
	
	Dsuper1 = Dsuper1 %>% mutate(insv_CCAATTGG = gsub("^;|;$","",insv_CCAATTGG,perl=T))  %>% mutate(insv_CCAATAAG=gsub("^;|;$","",insv_CCAATAAG,perl=T))  %>% 
			mutate(CP190=gsub("^;|;$","",CP190,perl=T)) %>%  mutate(GAGA=gsub("^;|;$","",GAGA,perl=T)) %>% 
			mutate(BEAF32=gsub("^;|;$","",BEAF32,perl=T))  %>%  mutate(EBox=gsub("^;|;$","",EBox,perl=T)) %>%  mutate(SuH=gsub("^;|;$","",SuH,perl=T)) %>%
			mutate(TSS_up2k=gsub("^;|;$","",TSS_up2k,perl=T), five_prime_UTR=gsub("^;|;$","",five_prime_UTR,perl=T),three_prime_UTR=gsub("^;|;$","",three_prime_UTR,perl=T)) %>% 
			mutate(CDS=gsub("^;|;$","",CDS,perl=T), intron=gsub("^;|;$","",intron,perl=T))
	
	Dsuper = select(dm, mergedPeakId, coord, Total.subpeaks) %>% distinct()  %>% inner_join(Dsuper1)  %>% inner_join(Dsuper2) %>%
			mutate(maxScore=apply(cbind(wtElba1__elba1Elba1.maxPeakscore, wtElba2__elba2Elba2.maxPeakscore, wtElba3__elba3Elba3.maxPeakscore, wtInsv__insvInsv.maxPeakscore), 1,  max, na.rm = TRUE)) %>% 
			arrange(desc(maxScore))
	

	dd = read.delim(peakannofile, stringsAsFactors=F, header=T)
	colnames(dd)[1] = "PeakID"
	dd$Focus.Ratio.Region.Size = gsub("_peaks.txt","", dd$Focus.Ratio.Region.Size)
	dd = dd %>% filter(PeakID %in% Dsuper$mergedPeakId)  %>%
	        mutate(Annotation=gsub(" UTR", "UTR", Annotation)) %>% 
			mutate(genomic_region=gsub("\\s+(.*)", "", Annotation)) %>% 
			arrange(desc(Peak.Score))  %>%
			select(PeakID, Gene.Name, genomic_region, Distance.to.TSS) 
	
	Dsuper2 = Dsuper %>% select(mergedPeakId, coord, insv_CCAATTGG, insv_CCAATAAG, CP190, GAGA,BEAF32,EBox, SuH, matches("peakid"),maxScore)%>%
			left_join(dd, by=c("mergedPeakId"="PeakID")) %>% 
			select(mergedPeakId, coord, genomic_region, Gene.Name, Distance.to.TSS,  everything())
	
	Dsuper2 = Dsuper2 %>% mutate(motifs = paste(insv_CCAATTGG,insv_CCAATAAG,CP190,GAGA,sep=";"))
	Dsuper2 = Dsuper2 %>% mutate(motifs=gsub("\\;{2,}",";", motifs)) %>% mutate(motifs=gsub("^\\;{1,}","", motifs)) %>% mutate(motifs=gsub("\\;{1,}$","", motifs))
	Dsuper2[Dsuper2$motifs %in% "", "motifs"] = "nomotifs"
	
	Dsuper2 = Dsuper2 %>% mutate(motifsCombine = "")
	Dsuper2[grepl("insv_CCAATTGG;insv_CCAATAAG", Dsuper2$motifs) &  (Dsuper2$motifsCombine %in% ""), "motifsCombine"] = "insv_CCAATTGG;insv_CCAATAAG"
	Dsuper2[grepl("insv_CCAATTGG", Dsuper2$motifs) &  (Dsuper2$motifsCombine %in% ""), "motifsCombine"] = "insv_CCAATTGG"
	Dsuper2[grepl("insv_CCAATAAG", Dsuper2$motifs) &  (Dsuper2$motifsCombine %in% ""), "motifsCombine"] = "insv_CCAATAAG"
	table(Dsuper2$motifsCombine)
	
	
	dm = dm %>% select(matches("mergedPeakId|Elba|Insv"))
	dm_peak = NULL
	flist = list()
	for (ii in 2:ncol(dm)) {
		bn = colnames(dm)[ii]
		dm_i = dm[, c("mergedPeakId", bn)]
		colnames(dm_i) = c("mergedPeakId", "peak")
		dm_i = dm_i %>% filter(peak != "")
		
		flist[[bn]] = unique(dm_i$mergedPeakId)
		dm_peak_i = data.frame(mergedPeakId=unique(dm_i$mergedPeakId), factor_motif=bn, stringsAsFactors=F)
		dm_peak = rbind(dm_peak, dm_peak_i)
	}
#	save(flist,dm_peak,file=paste(RData_dir,"Elba_merge_all_mutant_peaks.RData",sep="") )
	
	
	save(Dsuper,Dsuper2,dd,flist,dm_peak,file=savef)		
}


############################################################################################################
################################################# merge diff peak ##########################################

setwd(peak_dir)
macs_dirs = list.dirs(path = ".", full.names = FALSE, recursive = FALSE)

######## 4 major binding sites: wt vs cognate mutant ########
indirs =  macs_dirs[grepl("wtElba1_q20__elba1Elba1|wtElba2_q20__elba2Elba2|wtElba3_q20__elba3Elba3|wtInsv_q20__insvInsv",macs_dirs)]

outdir = paste("merge__major", sep="")
dir.create(outdir)
setwd(outdir)
mergeset(indirs, dist="given", nm="merge_major")
mergeset_anno(indirs, outdir,dist="given", nm="merge_major", savef = paste(data_dir, "chip_supertable_with_RNAseq_correctedMotif.RData",sep=""))


########  wt+mutant vs cognate mutant: 16 ########
setwd(peak_dir)
indirs = macs_dirs[grepl("Elba1_q20__elba1Elba1|Elba3_q20__elba3Elba3|Elba2_q20__elba2Elba2|Insv_q20__insvInsv",macs_dirs)]
outdir = paste("merge__major_ElbaAntibody", sep="")
dir.create(outdir)
setwd(outdir)
mergeset(indirs, dist="given", nm="mergemajor_ElbaAntibody")
mergeset_anno(indirs, outdir,dist="given", nm="mergemajor_ElbaAntibody", savef = paste(data_dir, "chip_ElbaAntibody_supertable_with_RNAseq_correctedMotif.RData",sep=""))



############################################################################################################
################################################# chipseq subsets ##########################################


load(file=paste(data_dir, "chip_ElbaAntibody_supertable_with_RNAseq_correctedMotif.RData",sep=""))	
Dsuper2 = Dsuper2[!is.infinite(Dsuper2$maxScore), ]

chip_subsets = list()

########################### Wt major ###########################
chip_subsets[["wtElba1"]] = Dsuper2 %>% filter(!is.na(wtElba1__elba1Elba1.peakid))
chip_subsets[["wtElba2"]] = Dsuper2 %>% filter(!is.na(wtElba2__elba2Elba2.peakid))
chip_subsets[["wtElba3"]] = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid))
chip_subsets[["wtInsv"]] = Dsuper2 %>% filter(!is.na(wtInsv__insvInsv.peakid))


########################### Elba3 #############################
elba3_elba12dep = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid)  & is.na(elba1Elba3__elba3Elba3.peakid) & is.na(elba2Elba3__elba3Elba3.peakid)) 
elba3_elba12indep = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid)  & (!is.na(elba1Elba3__elba3Elba3.peakid) | !is.na(elba2Elba3__elba3Elba3.peakid))) 
elba12_gained = Dsuper2 %>% filter(is.na(wtElba3__elba3Elba3.peakid)  & (!is.na(elba1Elba3__elba3Elba3.peakid) | !is.na(elba2Elba3__elba3Elba3.peakid))) 		
elba1_gained = Dsuper2 %>% filter(is.na(wtElba3__elba3Elba3.peakid)  & (!is.na(elba1Elba3__elba3Elba3.peakid))) 	
elba1_indep = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid)  & (!is.na(elba1Elba3__elba3Elba3.peakid))) 	

chip_subsets[["elba3_elba12indep"]] = elba3_elba12indep
chip_subsets[["elba3_elba12dep"]] = elba3_elba12dep
chip_subsets[["elba12_gained"]] =elba12_gained
chip_subsets[["elba1_gained"]] = elba1_gained
chip_subsets[["elba1_indep"]] = elba1_indep

########################### Elba1 #############################
elba1_wt = Dsuper2 %>% filter((!is.na(wtElba1__elba1Elba1.peakid)))
elba1_elba2mut = Dsuper2 %>% filter((!is.na(elba2Elba1__elba1Elba1.peakid)) )
elba3_elba2mut = Dsuper2 %>% filter((!is.na(elba2Elba3__elba3Elba3.peakid)))

chip_subsets[["elba1_wt"]] = elba1_wt
chip_subsets[["elba1_elba2mut"]] = elba1_elba2mut
chip_subsets[["elba3_elba2mut"]] = elba3_elba2mut

########################### Wt ################################

elba3_wt = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid) & is.na(wtInsv__insvInsv.peakid)) 
insv_wt = Dsuper2 %>% filter(!is.na(wtInsv__insvInsv.peakid) & is.na(wtElba3__elba3Elba3.peakid)) 
elba3_insv_wt = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid) & !is.na(wtInsv__insvInsv.peakid)) 

chip_subsets[["elba3_wt_noinsv"]] = elba3_wt
chip_subsets[["insv_wt_noElba3"]] = insv_wt
chip_subsets[["elba3_insv_wt"]] = elba3_insv_wt


########################### Elba3-only,  insv-only, elba1/2/3 overlap, 4 factor overlap ###########################
elba3_elba12dep = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid)  & is.na(elba1Elba3__elba3Elba3.peakid) & is.na(elba2Elba3__elba3Elba3.peakid)) 
elba3_only = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid) & is.na(wtElba1__elba1Elba1.peakid) & is.na(wtElba2__elba2Elba2.peakid)  & is.na(wtInsv__insvInsv.peakid) ) 
insv_only = Dsuper2 %>% filter(!is.na(wtInsv__insvInsv.peakid) & is.na(wtElba1__elba1Elba1.peakid) & is.na(wtElba2__elba2Elba2.peakid)  & is.na(wtElba3__elba3Elba3.peakid) ) 
insv_nonunique = Dsuper2 %>% filter(!is.na(wtInsv__insvInsv.peakid) &  !mergedPeakId %in% insv_only$mergedPeakId) 
elba123_ovlp =  Dsuper2 %>% filter(!is.na(wtElba1__elba1Elba1.peakid) & !is.na(wtElba2__elba2Elba2.peakid)  & !is.na(wtElba3__elba3Elba3.peakid) ) 
elba3.only_elba12dep_ovlp = elba3_elba12dep %>% filter(wtElba3__elba3Elba3.peakid %in% elba3_only$wtElba3__elba3Elba3.peakid)
elba123_ovlp.ovlp_insv = elba123_ovlp  %>% filter(!is.na(wtInsv__insvInsv.peakid) )
elba123_ovlp.notovlp_insv = elba123_ovlp  %>% filter( is.na(wtInsv__insvInsv.peakid) )
elba3insv_noelba12 =  Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid)  & !is.na(wtInsv__insvInsv.peakid) & is.na(wtElba1__elba1Elba1.peakid) & is.na(wtElba2__elba2Elba2.peakid) ) 
elba13_noelba2_noinsv = Dsuper2 %>% filter((!is.na(wtElba3__elba3Elba3.peakid) | is.na(wtElba1__elba1Elba1.peakid)) & is.na(wtElba2__elba2Elba2.peakid)  & is.na(wtInsv__insvInsv.peakid) ) 

chip_subsets[["elba3_only"]] = elba3_only
chip_subsets[["insv_only"]] = insv_only
chip_subsets[["insv_nonunique"]] = insv_nonunique
chip_subsets[["elba123_ovlp"]] = elba123_ovlp
chip_subsets[["elba123_ovlp.ovlp_insv"]] = elba123_ovlp.ovlp_insv
chip_subsets[["elba123_ovlp.notovlp_insv"]] = elba123_ovlp.notovlp_insv
chip_subsets[["elba3.only_elba12dep_ovlp"]] = elba3.only_elba12dep_ovlp
chip_subsets[["elba3insv_noelba12"]] = elba3insv_noelba12
chip_subsets[["elba13_noelba2_noinsv"]] = elba13_noelba2_noinsv

save(chip_subsets, file =paste(data_dir,"chip_elba_supertable_subsets.RData",sep=""))

########################################################################################################################################
######################### chipseq subsets: top 200 peaks with motif, peak with/without motifs ##########################################

get_genes <- function (genes) {
	mygenes = unlist(strsplit(genes, split=";"))
	mygenes = unique(mygenes)
	mygenes
}

get_myset <-function (peak_l_uni) {
	# peak_l_uni = diffpeak_peak_l_uni
	mysets = list()
	bns = names(peak_l_uni)[regexpr("q20", names(peak_l_uni), perl=T)!=-1 ]
	for (bn in bns) {
		dpeak0 = peak_l_uni[[bn]]
		colnames(dpeak0)[colnames(dpeak0) == "PeakID"] = "peakid"
		bn2 = gsub("_q20","", bn, perl=T)
		
		dpeak = dpeak0 %>% arrange(desc(peakscore)) %>% dplyr::select(peakid, nearestTSS,dist2NearestTSS, peakscore, starts_with("insv"),"GAGA", "CP190" )
		dpeak_insv = filter(dpeak, (insv_CCAATTGG != "" | insv_CCAATAAG != ""))
		dpeak_insv_CCAATTGG = filter(dpeak, (insv_CCAATTGG != "" & insv_CCAATAAG == ""))
		dpeak_insv_CCAATAAG = filter(dpeak, (insv_CCAATAAG != "" & insv_CCAATTGG == "" ))
		dpeak_insvNOT = filter(dpeak, (insv_CCAATTGG == "" & insv_CCAATAAG == ""))
		
		dpeak_CP190 = filter(dpeak, (CP190 != "" & (insv_CCAATTGG == "" & insv_CCAATAAG == "")))
		dpeak_GAGA = filter(dpeak, (GAGA != "" & (insv_CCAATTGG == "" & insv_CCAATAAG == "")))
		dpeak_insv_GAGA = filter(dpeak, (GAGA != "" & (insv_CCAATTGG != "" | insv_CCAATAAG != "")))
		
		
		dpeak_tss2k = dpeak[abs(dpeak$dist2NearestTSS) <= 2000, ]
		dyy0 = dpeak0 %>% filter(peakid %in% dpeak_tss2k$peakid) %>% dplyr::select(peakid,peakscore, nearestTSS, coord)
#		create_beds(bn2, dyy0, nm = "insv")

		dpeak_tss2k_insv = dpeak_insv %>% filter(peakid %in% dpeak_tss2k$peakid) %>% mutate(ntiles = ntile(peakscore,5))
		dpeak_tss2k_insv_top200 = dpeak_tss2k_insv[1:200,]
		dyy1 = dpeak0 %>% filter(peakid %in% dpeak_tss2k_insv_top200$peakid) %>% dplyr::select(peakid,peakscore, nearestTSS, coord)
#		create_beds(bn2, dyy1, nm = "insv_top200")
		
		dpeak_tss2k_insv_top100 = dpeak_tss2k_insv[1:100,]
		
		dpeak_tss2k_insvNOT = filter(dpeak_insvNOT, peakid %in% dpeak_tss2k$peakid)
		dyy2 = dpeak0 %>% filter(peakid %in% dpeak_tss2k_insvNOT$peakid) %>% dplyr::select(peakid,peakscore, nearestTSS, coord)
#		create_beds(bn2, dyy2, nm = "insvNOT")
		
		dpeak_tss2k_insv_CCAATTGG = filter(dpeak_insv_CCAATTGG, peakid %in% dpeak_tss2k$peakid)
		dpeak_tss2k_insv_CCAATAAG = filter(dpeak_insv_CCAATAAG, peakid %in% dpeak_tss2k$peakid)
		dpeak_tss2k_CP190 = filter(dpeak_CP190, peakid %in% dpeak_tss2k$peakid)
		dpeak_tss2k_GAGA = filter(dpeak_GAGA, peakid %in% dpeak_tss2k$peakid)
		dpeak_tss2k_insv_GAGA = filter(dpeak_insv_GAGA, peakid %in% dpeak_tss2k$peakid)
		
		mysets[[paste(bn2, "__insv_tss2k", sep="")]] = get_genes(dpeak_tss2k_insv$nearestTSS)
		mysets[[paste(bn2, "__insv_tss2k_top100", sep="")]] = get_genes(dpeak_tss2k_insv_top100$nearestTSS)
		mysets[[paste(bn2, "__insv_tss2k_top200", sep="")]] = get_genes(dpeak_tss2k_insv_top200$nearestTSS)
		mysets[[paste(bn2, "__insvNOT_tss2k", sep="")]] = get_genes(dpeak_tss2k_insvNOT$nearestTSS)
		
		for (jj in 1:5) {
			xx = dpeak_tss2k_insv %>% filter(ntiles == jj)
			mysets[[paste(bn2, "__insv_tss2k_",jj,"tile", sep="")]] = get_genes(xx$nearestTSS)
		}
			
	}
	print(lapply(mysets, length))
	mysets
}


load( file = paste(data_dir,"diffpeaks_anno_q20TRUE_spikeinFALSE_correctedMotif.RData",sep=""))	

diffpeak_peak_l_uni = peak_l_uni
ww = names(diffpeak_peak_l_uni)
ww = gsub("diff__|_q20", "", ww, perl=T)
ww = gsub("elba\\d+|insv|wt", "",ww, perl=T)
ww1 =  gsub("(.*)__(.*)", "\\1",ww, perl=T)
ww2 =  gsub("(.*)__(.*)", "\\2",ww, perl=T)
diffpeak_peak_l_uni = diffpeak_peak_l_uni[names(diffpeak_peak_l_uni)[which(ww1 == ww2)]]
diffpeak_peak_l_uni = diffpeak_peak_l_uni[names(diffpeak_peak_l_uni) %in% c("diff__wtElba1_q20__elba1Elba1_q20","diff__wtElba2_q20__elba2Elba2_q20","diff__wtElba3_q20__elba3Elba3_q20","diff__wtInsv_q20__insvInsv_q20")]

diffpeak_sets = get_myset(diffpeak_peak_l_uni)

names(diffpeak_sets) = gsub("diff__", "", names(diffpeak_sets), perl=T)
save(diffpeak_sets,  file = paste(file=paste(data_dir, "genesets2.RData",sep="")))













