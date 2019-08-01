#! /bin/env Rscript

library(Biostrings)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(venneuler)
library(getopt)
library(magrittr)

current_dir = "/home/jwen/projects/qi/elba/ELBA/"
setwd(current_dir)

data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")

##########################################################################
########################### functions ####################################
##########################################################################

plot_motif_frac <- function (dovlp) {
	
	dovlp$ss= factor(dovlp$ss, levels = unique(dovlp$ss))
	dovlp$motifs= factor(dovlp$motifs, levels = unique(dovlp$motifs))
	
	qqq= ggplot(data =dovlp, aes(x = ss, y = frac,fill=motifs)) +  geom_bar(stat="identity",width=0.9) 
	qqq = qqq + 
			scale_fill_manual(values =c( "gray", "orange", "red", "dodgerblue"),name = "")+
			labs(x="",y="Fraction of peaks") +
			theme_bw()+
			ggtitle("") + 
			theme(legend.position = "right",axis.text.x = element_text(angle =-45,hjust=0,vjust=1))
	
	print(qqq) 		  
}

get_motif_frac_chipnexus <-function (dd, myset) {
	
	dd = dd[,c("mergedPeakId", paste(myset,".motifs",sep=""))]
	colnames(dd) = c("mergedPeakId", "motifs")
	
	dplot = dd %>% 
			dplyr::select(mergedPeakId, motifs) %>% distinct() %>% group_by(motifs) %>% 
			summarise(npeaks=n_distinct(mergedPeakId)) %>% filter(!is.na(motifs)) %>% ungroup %>%
			mutate(frac=npeaks/sum(npeaks)) %>% mutate(ss = myset)  %>% data.frame  %>% arrange(ss, frac)
	dplot[dplot$motifs %in% "", "motifs"]	= "others"
	dplot2 = rbind(dplot[dplot$motifs %in% "others",],  dplot[dplot$motifs %in% "insv_CCAATTGG;insv_CCAATAAG",],  dplot[dplot$motifs %in% "insv_CCAATAAG",], dplot[dplot$motifs %in% "insv_CCAATTGG",] )
	dplot2	
}

get_motif_frac_chipseq <-function (dd, myset) {
	dplot = dd %>%  mutate(motifs=motifsCombine) %>%
			dplyr::select(mergedPeakId, motifs) %>% distinct() %>% group_by(motifs) %>% 
			summarise(npeaks=n_distinct(mergedPeakId)) %>% filter(!is.na(motifs)) %>% ungroup %>%
			mutate(frac=npeaks/sum(npeaks)) %>% mutate(ss = myset)  %>% data.frame  %>% arrange(ss, frac)
	dplot[dplot$motifs %in% "", "motifs"]	= "others"
	dplot2 = rbind(dplot[dplot$motifs %in% "others",],  dplot[dplot$motifs %in% "insv_CCAATTGG;insv_CCAATAAG",],  dplot[dplot$motifs %in% "insv_CCAATAAG",], dplot[dplot$motifs %in% "insv_CCAATTGG",] )
	dplot2	
}

plot2chip_ovlp_frac <- function (dovlp) {
	
	qqq= ggplot(data =dovlp, aes(x = motifs, y = frac,fill=ss)) +  geom_bar(stat="identity",width=0.9) #geom_area(stat="identity") #
	qqq = qqq + 
			scale_fill_manual(values = c("orange","forestgreen", "darkmagenta"),name = "")+
			labs(x="",y="fraction of peaks") +
			facet_wrap(~ myfac,scales = "free_x",ncol=4) +
			theme_bw()+
			theme(legend.position = "top",axis.text.x = element_text(angle =-45,hjust=0,vjust=1))
	
	print(qqq) 		  
}

plot_venneuler <- function (dplot, mytitle="") {
	colnames(dplot) = c("elements", "sets")
	v <- venneuler(dplot)
	plot(v, main=mytitle, col= c("blue","green", "magenta",  "black"))
}


###########################################################################
##### supp figure 4A,C ChIP-seq and ChIP-nexus motif fraction comparison ###
###########################################################################

load(file=paste(data_dir, "chip_ElbaAntibody_supertable_with_RNAseq_correctedMotif.RData",sep=""))	
# Dsuper2,dd,
dfrac_chipseq = rbind(
		get_motif_frac_chipseq(dd=Dsuper2 %>% filter((!is.na(wtElba1__elba1Elba1.peakid))), myset="Elba1_ChIPseq"),
		get_motif_frac_chipseq(dd=Dsuper2 %>% filter((!is.na(wtElba2__elba2Elba2.peakid))), myset="Elba2_ChIPseq"), 
		get_motif_frac_chipseq(dd=Dsuper2 %>% filter((!is.na(wtElba3__elba3Elba3.peakid))), myset="Elba3_ChIPseq"), 
		get_motif_frac_chipseq(dd=Dsuper2 %>% filter((!is.na(wtInsv__insvInsv.peakid))), myset="Insv_ChIPseq")
)

load(file =  paste(data_dir, "chipnexus_q10_flank0_given_anno.RData",sep=""))
dmerge = dd %>% dplyr::rename(mergedPeakId="PeakID")
dfrac_chipnexus = rbind(
		get_motif_frac_chipnexus(dd=dmerge %>% filter((!Elba1_ChIPnexus %in% "") ), myset="Elba1_ChIPnexus"),
		get_motif_frac_chipnexus(dd=dmerge %>% filter((!Elba2_ChIPnexus %in% "") ), myset="Elba2_ChIPnexus"), 
		get_motif_frac_chipnexus(dd=dmerge %>% filter((!Elba3_ChIPnexus %in% "")), myset="Elba3_ChIPnexus"), 
		get_motif_frac_chipnexus(dd=dmerge %>% filter((!Insv_ChIPnexus %in% "")), myset="Insv_ChIPnexus")
)


load(file = paste(data_dir,"merge2chip_q10_given_anno.RData",sep=""))
dmerge = dd %>% mutate(mergedPeakId = PeakID)

dmerge = dmerge %>% mutate(motifs = paste(insv_CCAATTGG,insv_CCAATAAG,CP190,GAGA,sep=";")) %>% 
         mutate(motifs=gsub("\\;{2,}",";", motifs)) %>% mutate(motifs=gsub("^\\;{1,}","", motifs)) %>% 
         mutate(motifs=gsub("\\;{1,}$","", motifs))
dmerge[dmerge$motifs %in% "", "motifs"] = "nomotifs"
dmerge = dmerge %>% mutate(motifsCombine = "")
dmerge[grepl("insv_CCAATTGG;insv_CCAATAAG", dmerge$motifs) &  (dmerge$motifsCombine %in% ""), "motifsCombine"] = "insv_CCAATTGG;insv_CCAATAAG"
dmerge[grepl("insv_CCAATTGG", dmerge$motifs) &  (dmerge$motifsCombine %in% ""), "motifsCombine"] = "insv_CCAATTGG"
dmerge[grepl("insv_CCAATAAG", dmerge$motifs) &  (dmerge$motifsCombine %in% ""), "motifsCombine"] = "insv_CCAATAAG"
table(dmerge$motifsCombine)


dfrac_merge2chip = rbind(
		get_motif_frac_chipseq(dd=dmerge %>% filter((!wtElba1__elba1Elba1 %in% "") & (!Elba1_ChIPnexus %in% "")), myset="Elba1_merge2chip"), 
		get_motif_frac_chipseq(dd=dmerge %>% filter((!wtElba2__elba2Elba2 %in% "") & (!Elba2_ChIPnexus %in% "")), myset="Elba2_merge2chip"), 
		get_motif_frac_chipseq(dd=dmerge %>% filter((!wtElba3__elba3Elba3 %in% "") & (!Elba3_ChIPnexus %in% "")), myset="Elba3_merge2chip"), 
		get_motif_frac_chipseq(dd=dmerge %>% filter((!wtInsv__insvInsv %in% "") & (!Insv_ChIPnexus %in% "")), myset="Insv_merge2chip")
)


outf = paste(res_dir, "supp_fig4ac.chipseq_main_motif_frac_barplot.pdf", sep="")
pdf(file=outf,width=5,height=4)
plot_motif_frac(dfrac_chipseq)
plot_motif_frac(dfrac_chipnexus)
plot_motif_frac(dfrac_merge2chip)
dev.off()
system(paste("pdfjam ",outf,"  --frame false --nup 3x1 --suffix 3up --outfile ",outf,sep=""))



###########################################################################
##### supp figure 4B ChIP-seq and ChIP-nexus peak overlap ###
###########################################################################

load(file = paste(data_dir,"merge2chip_q10_given_anno.RData",sep=""))

dmerge2peak = NULL
flist = list()
for (bn in c("Elba1", "Elba2", "Elba3", "Insv")) {
	
	if (grepl("Elba1",bn)) {
		ddi_chipnexus = dd %>% filter(!Elba1_ChIPnexus %in% "") 
		ddi_chipseq = dd %>% filter(!wtElba1__elba1Elba1 %in% "")
	} else if  (grepl("Elba2",bn)) {
		ddi_chipnexus = dd %>% filter(!Elba2_ChIPnexus %in% "") 
		ddi_chipseq = dd %>% filter(!wtElba2__elba2Elba2 %in% "")
	} else if  (grepl("Elba3",bn)) {
		ddi_chipnexus = dd %>% filter(!Elba3_ChIPnexus %in% "") 
		ddi_chipseq = dd %>% filter(!wtElba3__elba3Elba3 %in% "")
	} else if  (grepl("Insv",bn)) {
		ddi_chipnexus = dd %>% filter(!Insv_ChIPnexus %in% "") 
		ddi_chipseq = dd %>% filter(!wtInsv__insvInsv %in% "")
	}
	
	ddi_chipnexus = ddi_chipnexus %>% dplyr::select(PeakID, insv_CCAATTGG, insv_CCAATAAG ) %>% 
			mutate(ss=paste(bn, "_chipnexus",sep="")) %>% distinct()
	ddi_chipseq = ddi_chipseq %>% dplyr::select(PeakID, insv_CCAATTGG, insv_CCAATAAG) %>% 
			mutate(ss=paste(bn, "_chipseq",sep=""))  %>% distinct()
	ddi_overlap = ddi_chipnexus %>% filter(PeakID %in% ddi_chipseq$PeakID) %>% mutate(ss=paste(bn, "_overlap",sep=""))
	
	dmerge2peak = rbind(dmerge2peak,ddi_chipnexus,ddi_chipseq,ddi_overlap )
	
	flist[[paste(bn, "_chipnexus",sep="")]] = unique(ddi_chipnexus$PeakID)
	flist[[paste(bn, "_chipseq",sep="")]] = unique(ddi_chipseq$PeakID)
	flist[[paste(bn, "_overlap",sep="")]] = unique(ddi_overlap$PeakID)
}
colnames(dmerge2peak)[1] = "peakId"

dmerge2peak_sum = dmerge2peak %>%  
		group_by(ss, insv_CCAATTGG, insv_CCAATAAG) %>% summarise(n=n_distinct(peakId)) %>% ungroup %>%
		mutate(myfac=gsub("_(.*)","", ss))  %>%  mutate(ss=gsub("(.*)_(.*)","\\2", ss)) 
dmerge2peak_tot  = dmerge2peak  %>% 
		mutate(myfac=gsub("_(.*)","", ss))  %>%  mutate(ss=gsub("(.*)_(.*)","\\2", ss)) %>% 
		group_by(myfac,insv_CCAATTGG, insv_CCAATAAG)  %>% 
		summarise(ntot=n_distinct(peakId))
dmerge2peak_sum2 = dmerge2peak_sum %>% left_join(dmerge2peak_tot) %>%
		mutate(nfrac = n/ntot)  %>% mutate(motifs = paste( insv_CCAATTGG,"_", insv_CCAATAAG, sep="")) %>% 
		mutate(motifs=gsub("^\\_$","nomotifs", motifs)) %>%  mutate(motifs=gsub("^\\_|\\_$","", motifs)) %>% ungroup %>%
		dplyr::select(myfac, motifs,ss, nfrac, ntot ) %>% spread(ss,  nfrac)  %>%
		mutate(chipnexus=chipnexus-overlap, chipseq=chipseq-overlap) %>%
		gather(ss, frac,-matches("myfac|motifs|ntot"))  %>%
		arrange(myfac, motifs) %>% mutate(motifs = paste(motifs, "(",ntot,")", sep=""))  

outf = paste(res_dir, "supp_fig4b.chipseq_chipnexus_overlap_frac.pdf", sep="")
pdf(file=outf,width=6,height=5)
plot2chip_ovlp_frac(dmerge2peak_sum2 )
dev.off()


###########################################################################
##### supp figure 4D ChIP-seq and ChIP-nexus 4 factor overlap #############
###########################################################################



setwd(data_dir)
lnfiles_elba = list.files(path=".",pattern = "diff(.*)sumit.bed", full=F)

pdff = paste(res_dir, "supp_fig4d.Elba4factor_venn_chipseq.pdf", sep="")
pdf(file=pdff ,width=7,height=7)
for (dist in c(10, 25,50)) {
	lnfiles_str = paste(c(lnfiles_elba), sep="", collapse=" ")
	outf = paste("chipseq_",dist,".txt", sep="")
	vennf = paste("chipseq_",dist,"_venn.txt", sep="")
	system(paste("mergePeaks -d ", dist, " ", lnfiles_str ," -venn ",vennf, " > ", outf,sep=""))
	
	dout = read.delim(outf, stringsAsFactors=F, header=T)
	colnames(dout)[1] ="peakid"
	colnames(dout) = gsub("diff__|_q20_sumit_flank100.bed|_q20|_q20_sumit.bed","", colnames(dout) )
	
	flist = list()
	elba1 = dout %>% dplyr::filter(!wtElba1__elba1Elba1 %in% "") %>% dplyr::select(peakid) %>% distinct() %>% mutate(ss = paste("Elba1.allpeaks(",n(), ")",sep=""))
	elba2 = dout %>% dplyr::filter(!wtElba2__elba2Elba2 %in% "") %>% dplyr::select(peakid) %>% distinct() %>% mutate(ss = paste("Elba2.allpeaks(",n(), ")",sep=""))
	elba3 = dout %>% dplyr::filter(!wtElba3__elba3Elba3 %in% "") %>% dplyr::select(peakid) %>% distinct() %>% mutate(ss = paste("Elba3.allpeaks(",n(), ")",sep=""))
	insv = dout %>% dplyr::filter(!wtInsv__insvInsv %in% "") %>% dplyr::select(peakid) %>% distinct() %>% mutate(ss = paste("Insv.allpeaks(",n(), ")",sep=""))
	flist[["elba1"]] = elba1 %$%  peakid
	flist[["elba2"]] = elba2  %$%  peakid
	flist[["elba3"]] = elba3  %$%  peakid
	flist[["insv"]] = insv %$%  peakid
	
	dm_peak = elba1 %>% bind_rows(elba2)  %>% bind_rows(elba3)  %>% bind_rows(insv)
	dm_peak %>% group_by(ss) %>%  summarise(n_distinct(peakid))
	
	plot_venneuler(dplot=dm_peak, mytitle= paste("chipseq_peaks_dist", dist, sep=""))	

}
dev.off()
system(paste("pdfjam ",pdff,"  --frame false --nup 3x1 --suffix 3up --outfile ",pdff,sep=""))



lnfiles_elba = list.files(path=".",pattern = "(.*)_ChIPnexus_sumit_q(.*).bed", full=F)
pdff = paste(res_dir, "supp_fig4d.Elba4factor_venn_chipnexus.pdf", sep="")
pdf(file=pdff ,width=7,height=7)
for (dist in c(10, 25,50)) {
	lnfiles_str = paste(c(lnfiles_elba), sep="", collapse=" ")
	outf = paste("chipseq_",dist,".txt", sep="")
	vennf = paste("chipseq_",dist,"_venn.txt", sep="")
	system(paste("mergePeaks -d ", dist, " ", lnfiles_str ," -venn ",vennf, " > ", outf,sep=""))
	
	dout = read.delim(outf, stringsAsFactors=F, header=T)
	colnames(dout)[1] ="peakid"
	colnames(dout) = gsub("_sumit_q10.bed|_sumit_q5.bed","", colnames(dout) )
	
	flist = list()
	elba1 = dout %>% dplyr::filter(!Elba1_ChIPnexus %in% "") %>% dplyr::select(peakid) %>% distinct() %>% mutate(ss = paste("Elba1.allpeaks(",n(), ")",sep=""))
	elba2 = dout %>% dplyr::filter(!Elba2_ChIPnexus %in% "") %>% dplyr::select(peakid) %>% distinct() %>% mutate(ss = paste("Elba2.allpeaks(",n(), ")",sep=""))
	elba3 = dout %>% dplyr::filter(!Elba3_ChIPnexus %in% "") %>% dplyr::select(peakid) %>% distinct() %>% mutate(ss = paste("Elba3.allpeaks(",n(), ")",sep=""))
	insv = dout %>% dplyr::filter(!Insv_ChIPnexus %in% "") %>% dplyr::select(peakid) %>% distinct() %>% mutate(ss = paste("Insv.allpeaks(",n(), ")",sep=""))
	flist[["elba1"]] = elba1 %$%  peakid
	flist[["elba2"]] = elba2  %$%  peakid
	flist[["elba3"]] = elba3  %$%  peakid
	flist[["insv"]] = insv %$%  peakid
	
	dm_peak = elba1 %>% bind_rows(elba2)  %>% bind_rows(elba3)  %>% bind_rows(insv)
	dm_peak %>% group_by(ss) %>%  summarise(n_distinct(peakid))
	
	plot_venneuler(dplot=dm_peak, mytitle= paste("chipnexus_peaks_dist", dist, sep=""))	
	

}
dev.off()
system(paste("pdfjam ",pdff,"  --frame false --nup 3x1 --suffix 3up --outfile ",pdff,sep=""))

