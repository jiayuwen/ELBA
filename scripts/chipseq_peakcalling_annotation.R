#! /bin/env Rscript

library(Biostrings)
library(ggplot2)
library(getopt)
library(data.table)
library(gdata)
library(GSEABase)
library("annotate") 
library("GOstats") 
library("GO.db") 
library(RBGL)
library(biomaRt)
library(dplyr)
library(tidyr)


current_dir = "/home/jwen/projects/qi/elba/code_NC/"
setwd(current_dir)
data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")
peak_dir = paste(data_dir, "chipseq_peaks/",sep="")
genome_f =  paste(data_dir, "dm3.chrom.sizes",sep="")
genome_fasta = paste(data_dir, "dmel_noUextra.fasta",sep="") 

motif_f = paste(data_dir,"dm_TF_motif_qi.txt",sep="")
mmm = readLines(motif_f )
mmm = mmm[mmm != ""]
mmm = mmm[regexpr("MOTIF", mmm, perl=T) !=-1]
mmm = gsub("MOTIF ", "", mmm, perl=T)

##################################################################
##################################################################

ovlpBeds <- function (inbed, tarfile,do_strand=T, ovlapfrac=1) {
	exprf = tempfile(pattern = paste("tmp_",Sys.getpid(),sep=""), tmpdir = tmp_dir, fileext ="txt") 
	unlink(exprf)
	
	if (ovlapfrac == 1) {
		if (do_strand) {
			system(paste("intersectBed   -s -wo -a ",inbed," -b ",tarfile," > ",exprf,sep=""))
		} else {
			system(paste("intersectBed  -wo -a ",inbed," -b ",tarfile," > ",exprf,sep=""))
		}
	} else {
		if (do_strand) {
			system(paste("intersectBed  -f ", ovlapfrac,"  -s -wo -a ",inbed," -b ",tarfile," > ",exprf,sep=""))
		} else {
			system(paste("intersectBed -f ", ovlapfrac,"  -wo -a ",inbed," -b ",tarfile," > ",exprf,sep=""))
		}
	}
	
	if (file.info(exprf)$size == 0) {
		outi = NULL
	} else {
		outi = read.delim(file=exprf,header=F,stringsAsFactors=F)
	}
	
	unlink(exprf)
	outi
}

annotate_fea <- function (dd, inbed,anno_dir, genomic_regions) {
	
	for (nm in  genomic_regions ){
		
		cat(nm,"\n")
		if (nm =="known_motif_oligoMatch") {
			tarfile = paste("/home/seqdata/people/jwen/project/qi/jaaved_motif/",nm,".bed",sep="")
			outi_gene = ovlpBeds(inbed, tarfile, do_strand = F, ovlapfrac=1)
			outi_gene$V10 = gsub("\\+\\d+|\\-\\d+","", outi_gene$V10, perl=T)
		} else {
			tarfile = paste(anno_dir,nm,".bed",sep="")
			outi_gene = ovlpBeds(inbed, tarfile, do_strand = F, ovlapfrac=0.5)
		}
		
		
		outi_gene  = unique(outi_gene )
		
		if (!is.null(outi_gene)) {
			
			if (nm %in% c("TSS_up2k", "gene_corrd","known_motif_oligoMatch")) {
				outi_gene = subset(outi_gene, select=c(V4, V10))
			} else {
				outi_gene = subset(outi_gene, select=c(V4, V13))
			}
			colnames(outi_gene) = c("peakid","region")
			
			dt<-as.data.table(outi_gene)
			setkey(dt,peakid)
			xx1 = dt[,list(ss=paste(unique(region), collapse=";",sep="")),by=peakid]
			outi_gene  =data.frame(xx1)			
			colnames(outi_gene) = c("peakid", nm)
			
			dd = left_join(dd, outi_gene)
		}
	}
	
	dd
}

dist2TSS <- function (dd, inbed,nm) {
	tarfile = paste(anno_dir,nm,".bed",sep="")
	
	
	exprf = tempfile(pattern = paste("tmp_",Sys.getpid(),sep=""), tmpdir = tmp_dir, fileext ="txt") 
	
#	system(paste('closestBed -D a -t all -g ', genome_f, ' -a ',inbed,' -b ',tarfile,'  > ',exprf,sep=''))
	system(paste('closestBed -D a -t all  -a ',inbed,' -b ',tarfile,'  > ',exprf,sep=''))
	
	if (file.info(exprf)$size == 0) {
		outi = NULL
	} else {
		outi = read.delim(file=exprf,header=F,stringsAsFactors=F)
	}
	
	#outi[outi$V4 == "ham_q20_peak_6954", ]
	
	outi = subset(outi, select=c(V4,V10, V12, V13))
	colnames(outi) = c("peakid", "nearestTSS","nearestTSSStrand", "dist2NearestTSS")
	
	dt<-as.data.table(outi)
	setkey(dt,peakid)
	xx1 = dt[,list(nearestTSS=paste(unique(nearestTSS), collapse=";",sep=""), nearestTSSStrand=paste(unique(nearestTSSStrand), collapse=";",sep=""), dist2NearestTSS= unique(dist2NearestTSS)),by=peakid]
	outi  =data.frame(xx1)			
	
	unlink(exprf)
	outi
	
}

dist2peak <- function (dd, inbed) {
	
	
	exprf = tempfile(pattern = paste("tmp_",Sys.getpid(),sep=""), tmpdir = tmp_dir, fileext ="txt") 
	
#	system(paste('closestBed -D a -t all -g ', genome_f, ' -a ',inbed,' -b ',tarfile,'  > ',exprf,sep=''))
	system(paste('closestBed -D a -io -t all  -a ',inbed,' -b ',inbed,'  > ',exprf,sep=''))
	
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

annotate_genes <- function (dnp_anno) {
	load(file=paste("/seqdata2/people/jwen/projects/12files/RData/gene_anno_merge_fb.RData",sep="") )
#danno,
	danno$desc = gsub("\\\\", "_",  danno$desc,perl=T)
	
#danno_desc = unique(subset(danno, select=c(primary_FBgn, gene_symbol, annotation_ID, Interacting_genes_symbol,Interacting_genes_FBgn, gene_symbol2,gene_FBgn2, desc)))
	danno_desc = unique(subset(danno, select=c(gene_symbol, annotation_ID, Interacting_genes_symbol,gene_symbol2,desc)))
	
	dnp_anno_1 = merge(dnp_anno, danno_desc, by.x="nearestTSS",by.y="gene_symbol")
	dnp_anno_2 = dnp_anno[!dnp_anno$peakid %in% dnp_anno_1$peakid, ]
	dnp_anno_2 = merge(dnp_anno_2, danno_desc, by.x="nearestTSS",by.y="annotation_ID")
	xx = intersect(colnames(dnp_anno_1 ),colnames(dnp_anno_2) )
	dnp_anno_12 = rbind(dnp_anno_1[,xx],dnp_anno_2[,xx])
	dnp_anno_3 = dnp_anno[!dnp_anno$peakid %in% dnp_anno_12$peakid, ]
	dnp_anno_3 = mutate(dnp_anno_3, Interacting_genes_symbol=NA, gene_symbol2=NA, desc=NA)
	
	
#	xx3 = setdiff(colnames(dnp_anno_12),colnames(dnp_anno_3))
#	tt[tt$nearestTSS == "HLHm7", ]
#	danno[regexpr("CG8361", danno$secondary_annotation_IDs, perl=T) !=-1, ]
#	danno[regexpr("HLHm7", danno$gene_symbol, perl=T) !=-1, ]
	
	dnp_anno_12  = rbind(dnp_anno_12 , dnp_anno_3)
	yy = c(colnames(dnp_anno), setdiff(colnames(dnp_anno_12),colnames(dnp_anno)))
	dnp_anno_12 = dnp_anno_12[, yy]
	dnp_anno_12
	
}

create_beds <- function(bn, dnp, nm, dofasta=F) {
	
	dnp = dnp[with(dnp, order(chr,chr_start, chr_end,strand)),]
	dnp[,5] = round(dnp[,5])
	outf = paste(macs2_dir,bn,"/",bn, "_",nm,".bed", sep="")
	unlink(outf)
	write.table(dnp[,1:6],file=outf,row.names = F,col.names=F, quote = F,sep="\t",append = F)
	
#	outf_bb =paste(macs2_dir,bn,"/", bn,"_",nm,".bb", sep="")
#	system(paste("bedToBigBed ",outf," ",genome_f," ",outf_bb,sep=""))
	
	if (dofasta) {
		outfasta = paste(macs2_dir,bn,"/",bn, "_",nm,".fasta", sep="")
		system(paste("fastaFromBed -fi ",genome_fasta," -name -s -bed ",outf," -fo ",outfasta,sep="" ))
	}
	
}

get_parents <- function (ontology,goids){	
	cname = switch(ontology,
			P = "go_biological_process_id",
			F = "go_molecular_function_id",
			C = "go_cellular_component_id"
	)
	
#	goarray = get.go.bp
	parents = NULL
	gids = NULL
	for (gg in unique(goids) ) {
		
		if (ontology == "P") {
			og = GOGraph(gg, GOBPPARENTS)	
		} else if (ontology == "F") {
			og = GOGraph(gg, GOMFPARENTS)	
		} else if (ontology == "C") {
			og = GOGraph(gg, GOCCPARENTS)	
		}
#		checked.go = c("GO:0060111", "GO:0060112")
		op = as.character(rev(tsort(og))) # do topological sort to get more general categories at front of list
#		op = op[3:length(op)] # remove top level go classificatons
#		op = op[op %in% nn]
		op = op[!op %in% "all"]
		
#		parents = rbind(parents,cbind(gg,op))
		parents = c(parents,op)
		gids = c(gids,rep(gg,length(op)))
	}
#	go.parents = cbind(gids,parents)
#	colnames(go.parents) = c(geneIdType,cname)
	parents
}

assign_gene2go<- function (array) {
	geneId = as.character(array[1])
	ontology = as.character(array[2])
	
	get.go = goo[(as.character(goo$geneId) %in% as.character(geneId)) & (goo$catalog %in% ontology) , "GOId"]
	gos = get_parents(ontology,get.go)
	gos  = paste(gos, collapse=";" )
	
	gos
}

getGO_neuro <- function (dnp_anno) {
	genes_all = unique(dnp_anno$nearestTSS)
	goo = goo[goo$geneId %in% genes_all,]
	goo = unique(subset(goo, select=c(geneId,GOId,catalog)))
	goo = merge(goo, subset(gonerv_desc_biomaRt, select=c(go_id,name_1006)),by.x="GOId",by.y="go_id")
	
	dt<-as.data.table(goo)
	setkey(dt,geneId)
	xx1 = dt[,list(GOIds_neuro=paste(unique(GOId), collapse=";",sep=""),GOterms_neuro=paste(unique(name_1006), collapse=";",sep="")),by=geneId]
	xx1  =data.frame(xx1)	
	
	dnp_anno = merge(dnp_anno, xx1, by.x="nearestTSS", by.y="geneId", all.x=T)
	dnp_anno
}


################################################################################################
################################################################################################

#peak calling parameters
# macs2 callpeak -t  wt -c mutant --format=BAM  --gsize=dm  --qvalue=0.05   --cutoff-analysis --call-summits  


################################################################################################
################################################################################################
do_anno=T
if (do_anno) {
	
	genomic_regions = c("CDS", "five_prime_UTR", "three_prime_UTR", "intron","TSS_up2k","gene_corrd","pseudogene","enhancer", "protein_binding_site","regulatory_region") #,"known_motif_oligoMatch"
	
	
	setwd(peak_dir)
	macs_dirs = list.dirs(path = ".", full.names = FALSE, recursive = FALSE)
	macs_dirs = macs_dirs[regexpr("(.*)insv(.*)__insvInsv_q20|(.*)Insv(.*)__insvInsv_q20|(.*)elba1(.*)__elba1Elba1_q20|(.*)Elba1(.*)__elba1Elba1_q20|(.*)elba2(.*)__elba2Elba2_q20|(.*)Elba2(.*)__elba2Elba2_q20|(.*)elba3(.*)__elba3Elba3_q20|(.*)Elba3(.*)__elba3Elba3_q20|wtElba1_q20__elba2Elba1|wtElba1_q20__elba3Elba1|wtElba2_q20__elba1Elba2|wtElba2_q20__elba3Elba2|wtElba3_q20__elba1Elba3|wtElba3_q20__elba2Elba3",macs_dirs, perl=T) != -1]
	
	peak_l = list()
	peak_l_uni = list()
	
	for ( rr in sort(macs_dirs)) { 

		peaktype = "narrow"
		inRfile = paste(peak_dir, rr, "_macs2.RData",sep="")
		
		if (file.exists(inRfile)) {
			load(file=inRfile)
			#dcutoff, dnp, dsumit, 
			dnp = transform(dnp, coord= paste(chr, ":", chr_start, "-", chr_end, sep=""), stringsAsFactors=F)
			bn = gsub("_spikein", "", basename(rr), perl=T)
			
			if (peaktype == "broad") {
				dnp = subset(dnp, select=c(coord, peakid, peakscore, foldchange, pscore, qscore))
				narrowPeak_f = paste(rr,"/", bn, "_broadPeak.bed",sep="")	
			} else if (peaktype == "narrow"){
				dnp = subset(dnp, select=c(coord, peakid, peakscore, foldchange, pscore, qscore, sumitpos2peak))
				narrowPeak_f = paste(rr,"/", bn, "_narrowPeak.bed",sep="")
			}
			
			
			if (file.exists(narrowPeak_f)  && file.info(narrowPeak_f)$size != 0) {
				
				dnp_anno = annotate_fea(dd = dnp, inbed=narrowPeak_f, anno_dir,genomic_regions)

				# motifs from memechip and fimo
				motiff = paste(peak_dir, "fimo_newchip/macs2/", rr, "/fimo.txt", sep="")
				
				
				if (file.exists(motiff)) {
					
					dmotif  = read.delim(motiff,header=T, stringsAsFactors=F)
					
					if (nrow(dmotif) != 0) {
						
						colnames(dmotif) = c("motif.pattern", "peakid", "motif.start", "motif.end", "motif.strand","motif.score","motif.pval","motif.qval","motif.matchedSequence")		
						
						tt = dmotif %>% filter((motif.pattern %in% "insv_CCAATAAG"))
						table(tt$motif.matchedSequence)
						
						tt = dmotif %>% filter((motif.pattern %in% "insv_CCAATTGG"))
						table(tt$motif.matchedSequence)
						
						dmotif = dmotif %>% 
								filter((motif.pattern %in% "insv_CCAATTGG" & motif.matchedSequence %in% c("CCTATTGG","CCAATTAG","CTAATTGG","CCAATTGG")) | (motif.pattern %in% "insv_CCAATAAG" & motif.matchedSequence %in% c("CCGATAAG","CCAATAAG")) | !(motif.pattern %in% c("insv_CCAATTGG","insv_CCAATAAG")))

						dmotif2 = data.frame(peakid = unique(dnp_anno$peakid), stringsAsFactors=F)
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
						
						dnp_anno = left_join(dnp_anno, dmotif2)
					}
				}
							
				dnp_dist2TSS = dist2TSS(dd = dnp, inbed=narrowPeak_f,nm = "TSS")
				
				dnp_peak2peak = dist2peak(dd = dnp, inbed=narrowPeak_f)
				
				dnp_anno = merge(dnp_anno, dnp_dist2TSS, by.x="peakid",  by.y="peakid", all.x=T, sort=F)
				
				dnp_anno = merge(dnp_anno,dnp_peak2peak, by.x="peakid",  by.y="peakid", all.x=T, sort=F)

				
				dnp_anno = getGO_neuro(dnp_anno)
				
				dnp_anno = annotate_genes(dnp_anno)
								
				dt<-as.data.table(dnp_anno)
				setkey(dt,nearestTSS)
				xx1 = dt[,list(Interacting_genes_symbol=paste(unique(Interacting_genes_symbol), collapse=";",sep=""),gene_symbol2=paste(unique(gene_symbol2), collapse=";",sep="")),by=nearestTSS]
				xx1  =data.frame(xx1)			
				
				dnp_anno_uni = merge(unique(subset(dnp_anno,select=-c(Interacting_genes_symbol, gene_symbol2))), xx1, by.x="nearestTSS", by.y="nearestTSS")
				dnp_anno_uni = dnp_anno_uni[, colnames(dnp_anno)]
				
				
				dnp_anno = dnp_anno[with(dnp_anno, order(-peakscore, -qscore,-pscore)), ]
				dnp_anno_uni = dnp_anno_uni[with(dnp_anno_uni, order(-peakscore, -qscore,-pscore)), ]
				
				
				peak_l[[rr]] = dnp_anno
				
				peak_l_uni[[rr]] = dnp_anno_uni
			}
		}
	}
	
	save(peak_l,peak_l_uni, file = paste(data_dir,"diffpeaks_anno_q20TRUE_spikeinFALSE_correctedMotif.RData",sep=""))

}

add_motif=T
if (add_motif) {
	
	load(file= paste(data_dir,"diffpeaks_anno_q20TRUE_spikeinFALSE_correctedMotif.RData",sep=""))
	# peak_l,peak_l_uni, 
	
	
	for (rr in names(peak_l)) {
		
		dnp_anno = peak_l[[rr]] 
		dnp_anno_uni = peak_l_uni[[rr]] 
		
		motiff = paste(peak_dir, "fimo_newchip/macs2/", rr, "/fimo.txt", sep="")
		dmotif  = read.delim(motiff,header=T, stringsAsFactors=F)
		
		if (nrow(dmotif) != 0) {
			
			colnames(dmotif) = c("motif.pattern","motif_alt_id", "peakid", "motif.start", "motif.end", "motif.strand","motif.score","motif.pval","motif.qval","motif.matchedSequence")		
			
			tt = dmotif %>% filter((motif.pattern %in% "insv_CCAATAAG"))
			table(tt$motif.matchedSequence)
			
			tt = dmotif %>% filter((motif.pattern %in% "insv_CCAATTGG"))
			table(tt$motif.matchedSequence)
			
			dmotif = dmotif %>% 
					filter((motif.pattern %in% "insv_CCAATTGG" & motif.matchedSequence %in% c("CCTATTGG","CCAATTAG","CTAATTGG","CCAATTGG")) | (motif.pattern %in% "insv_CCAATAAG" & motif.matchedSequence %in% c("CCGATAAG","CCAATAAG")) | !(motif.pattern %in% c("insv_CCAATTGG","insv_CCAATAAG")))
			
			dmotif2 = data.frame(peakid = unique(dnp_anno$peakid), stringsAsFactors=F)
			for (oo in mmm) {
				
				hh = dmotif  %>% filter(motif.pattern== oo) %>% select(peakid,motif.pattern ) %>% distinct()
				colnames(hh) = c("peakid", oo)
				dmotif2 = dmotif2 %>% left_join(hh) 
#					
			}
			
			dmotif2[is.na(dmotif2$insv_CCAATTGG), "insv_CCAATTGG"] = "" 
			dmotif2[is.na(dmotif2$insv_CCAATAAG), "insv_CCAATAAG"] = "" 
			
			dmotif2 = mutate(dmotif2, CP190="", GAGA="", BEAF32="", SuH="", EBox="")
			dmotif2[(!is.na(dmotif2$CP190_NTGGCMACACTR)), "CP190"] = "CP190" 
			dmotif2[(!is.na(dmotif2$GAGA_AGAGAGMGAGAG)), "GAGA"] = "GAGA" 
			dmotif2[(!is.na(dmotif2$BEAF32_WHATCGATARY)), "BEAF32"] = "BEAF32" 
			dmotif2[(!is.na(dmotif2$SuH_TGTGGGAA)), "SuH"] = "SuH" 
			dmotif2[(!is.na(dmotif2$EBox_CAGCTGWTY)), "EBox"] = "EBox" 
			
			xxx = c("peakid", setdiff(colnames(dmotif2), colnames(dnp_anno) ))
			dmotif2 = dmotif2[, xxx]
			dnp_anno = left_join(dnp_anno, dmotif2)
			dnp_anno_uni = left_join(dnp_anno_uni, dmotif2)
			
			peak_l[[rr]] = dnp_anno
			peak_l_uni[[rr]] = dnp_anno_uni
			
		}
	}
	
	save(peak_l,peak_l_uni, file = paste(data_dir,"diffpeaks_anno_q20TRUE_spikeinFALSE_correctedMotif.RData",sep=""))
	
}




################################################################################################
################################################################################################
