#! /bin/env Rscript


library(Biostrings)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(venneuler)
library(getopt)
library(corrplot)
library(RColorBrewer)

current_dir = "/home/jwen/projects/qi/elba/ELBA/"
setwd(current_dir)

data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")


##########################################################################
######################################## functions #######################
##########################################################################


##### figure 5.A memechip for de novo motif search ######################


##### figure 5.B Elba/Insv overlap with other insulators #################

setwd(data_dir)
lnfiles_elba = list.files(path=".",pattern = "(.*)_merge2chip_given.bed", full=F)
lnfiles_others =  list.files(path=".",pattern = "(.*)flank100.bed", full=F)
lnfiles_str = paste(c(lnfiles_elba, lnfiles_others ), sep="", collapse=" ")

dist = "given"
outf = paste("chipchip_insulator_",dist,".txt", sep="")
vennf = paste("chipchip_insulator_",dist,"_venn.txt", sep="")
system(paste("mergePeaks -d ", dist, " ", lnfiles_str ," -venn ",vennf, " > ", outf,sep=""))


dout = read.delim(outf, stringsAsFactors=F, header=T)
colnames(dout)[1] ="mergeid"
colnames(dout) = gsub("X.home.jwen.projects.qi.elba.modencode_insulator..|_merge2chip_given.bed|_table854.bed|_NO_1_fdr_consecutive_peak_filtered_flank100.bed|_peaks.narrowPeak|_peaks_flank100.bed|__wt_mutant_cycle12.14", "", colnames(dout))
colnames(dout) = gsub("GSM\\d+\\_", "",colnames(dout))
dout = dout %>% dplyr::select(mergeid, "Elba3","Elba2","Elba1" ,"Insv","BEAF","CP190","CTCF_C","CTCF_N","GAF","MDJ4","suHw","suHw_pam") %>% gather(ss, peaks, -mergeid) %>% filter(!peaks %in% "") %>% dplyr::select(-peaks)
dnn = dout %>% group_by(ss)  %>% summarise(npeaks = n_distinct(mergeid)) %>% data.frame
dout = dout %>% left_join(dnn)  %>%  mutate(ss=paste(ss,"(", npeaks,")",sep=""))  %>%  dplyr::select(-npeaks)

xx = unique(dout$ss)
dmat = matrix(0, nrow=length(xx), ncol=length(xx))
rownames(dmat) = xx; colnames(dmat) = xx
dmat_frac = dmat


for (ii in 1:(length(xx)-1)) {
	for (jj in 2:length(xx)) {
		
		sii = xx[ii]
		sjj = xx[jj]
		
		dij = dout %>% filter(ss %in% c(sii, sjj)) 
		dij_indiv = dij %>% group_by(ss)  %>% summarise(npeaks = n_distinct(mergeid)) %>% data.frame
		dij_ovlp = dij %>% group_by(mergeid)  %>% summarise(npeaks = n_distinct(ss))  %>% dplyr::filter(npeaks > 1)  %>% n_distinct()
		
		for (gg in dij_indiv$ss) {
			dmat[gg, gg] = dij_indiv[dij_indiv$ss == gg, "npeaks"]
			dmat_frac[gg, gg] =1
		}
		dmat[sii, sjj] = dij_ovlp
		dmat[sjj, sii] = dij_ovlp
		
		dmat_frac[sii, sjj] = dij_ovlp/min(dij_indiv$npeaks)
		dmat_frac[sjj, sii] = dij_ovlp/min(dij_indiv$npeaks)
		
	}
}


pdff = paste(res_dir, "figure5b.chipchip_insulator.", dist, "_corrplot.pdf", sep="" )
pdf(file=pdff,width=10, height=10)
corrplot(dmat_frac*100, method = "pie",col= brewer.pal(11,"RdYlBu"), is.corr = FALSE, type="lower",cl.lim =c(0,100), title="", diag=TRUE,  tl.srt=45, tl.col="black") 
dev.off()


