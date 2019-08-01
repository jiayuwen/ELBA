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

plot_peak_motif <- function (dplot) {
	
	dplot = dplot %>% gather(xx, pct,-contains("ss")) %>% 
			filter(xx %in% "insv_motifs.frac") %>% mutate(genotype=gsub("__(.*)","", ss)) %>% 
			mutate(against=gsub("(.*)__","", ss)) %>% mutate(against=ifelse(!grepl("Input|IgG", against), "mutant",against))
	dplot = rbind(dplot[!grepl("IgG|Input", dplot$ss),], dplot[grepl("Input", dplot$ss),],dplot[grepl("IgG", dplot$ss),])
	dplot$against = factor(dplot$against, levels = unique(dplot$against ))
	
	qqq= ggplot(data =dplot, aes(x = genotype, y = pct, fill=factor(against))) + geom_bar(position="dodge",stat="identity",width=0.7)
	qqq = qqq + 
			labs(x="",y="proportion of insv motifs") +
			scale_fill_manual(values = c("blue","dimgray","gray"),name = "") +
			theme_classic()+
			theme(legend.position = "right",axis.text.x = element_text(angle =-45,hjust=0,vjust=1))
	
	print(qqq) 
	
}

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

get_motif_frac <-function (dd, myset) {
	dplot = dd %>%  mutate(motifs=motifsCombine) %>%
			dplyr::select(mergedPeakId, motifs) %>% distinct() %>% group_by(motifs) %>% 
			summarise(npeaks=n_distinct(mergedPeakId)) %>% filter(!is.na(motifs)) %>% ungroup %>%
			mutate(frac=npeaks/sum(npeaks)) %>% mutate(ss = myset)  %>% data.frame  %>% arrange(ss, frac)
	dplot[dplot$motifs %in% "", "motifs"]	= "others"
	dplot2 = rbind(dplot[dplot$motifs %in% "others",],  dplot[dplot$motifs %in% "insv_CCAATTGG;insv_CCAATAAG",],  dplot[dplot$motifs %in% "insv_CCAATAAG",], dplot[dplot$motifs %in% "insv_CCAATTGG",] )
	dplot2	
}


###########################################################################
##### supp figure 1C proportion of motifs mutant, input, igg ##############
###########################################################################

load(file = paste(data_dir,"Allpeaks_anno_q20TRUE_spikeinFALSE_correctedMotif.RData",sep=""))
# peak_l,peak_l_uni, 

Dsum = NULL
for (ii in names(peak_l_uni)) {
	
	ddi = peak_l_uni[[ii]]
	if (!"dm3-blacklist" %in% colnames(ddi)) {
		ddi = cbind(ddi, `dm3-blacklist` = NA)
	}
	if (!"insv_CCAATTGG" %in% colnames(ddi)) {
		ddi = cbind(ddi, insv_CCAATTGG = "")
	}
	if (!"insv_CCAATAAG" %in% colnames(ddi)) {
		ddi = cbind(ddi, insv_CCAATAAG = "")
	}

	ddi  = ddi %>% 
		   dplyr::select(peakid, `dm3-blacklist`,insv_CCAATTGG, insv_CCAATAAG,dist2NearestTSS) %>% 
           mutate(peakid=gsub("[a-z]$","", peakid))  %>% distinct()
   
	ddi_sum = bind_cols(
			ss = ii, 
			peaks_tot = n_distinct(ddi$peakid), 
			blacklist = n_distinct(ddi[!is.na(ddi$`dm3-blacklist`), "peakid"]), 
			insv_motifs = n_distinct(ddi[(!ddi$insv_CCAATTGG %in% "") | (!ddi$insv_CCAATAAG %in% "" ), "peakid"]), 
			npromoter2K = n_distinct(ddi[abs(ddi$dist2NearestTSS) <=2000, "peakid"]),
	        npromoter2K_insv_motif = n_distinct(ddi[abs(ddi$dist2NearestTSS) <=2000 & ((!ddi$insv_CCAATTGG %in% "") | (!ddi$insv_CCAATAAG %in% "" )), "peakid"])
	)
	Dsum = Dsum %>% bind_rows(ddi_sum)
}

Dsum = Dsum %>% mutate(ss=gsub("diff__|_q20","", ss)) %>% filter(!grepl("old", ss)) %>%  data.frame %>% 
	   mutate(blacklist.frac=blacklist/peaks_tot, insv_motifs.frac= insv_motifs/peaks_tot, npromoter.frac=npromoter2K/peaks_tot, npromoter2K_insv_motif.frac=npromoter2K_insv_motif/npromoter2K) %>% 
	   arrange(ss) %>% data.frame 

dplot = Dsum %>% dplyr::filter(grepl("(.*)wtInsv(.*)(wt|insvInsv)|(.*)wtElba1(.*)(wt|elba1Elba1)|(.*)wtElba2(.*)(wt|elba2Elba2)|(.*)wtElba3(.*)(wt|elba3Elba3)", ss))%>%
        dplyr::select(ss, insv_motifs.frac) %>%  mutate(nomotif = 1-insv_motifs.frac)
pdff = paste(res_dir, "supp_fig1c.npeak_wt_mutant_input_igg.pdf", sep="" )
pdf(file=pdff,width=4, height=6)
plot_peak_motif(dplot)
dev.off()


###########################################################################
##### supp figure 1D subset density center at peak summits #####
###########################################################################


###########################################################################
##### supp figure 1E subset motif fraction  #####
###########################################################################

load(file =paste(data_dir,"chip_elba_supertable_subsets.RData",sep=""))
# chip_subsets, 
chip_subsets = chip_subsets[c("elba123_ovlp.ovlp_insv", "elba123_ovlp.notovlp_insv", "elba3insv_noelba12", "insv_only", "elba3_only")]
names(chip_subsets) = c("Elba123_Insv","Elba123_noInsv","Elba3Insv_noElba12", "Insv-unique", "Elba3-unique")

dfrac = NULL
for (ii in names(chip_subsets)) {
	dfrac = dfrac %>% bind_rows(get_motif_frac(dd=chip_subsets[[ii]], myset=ii))
}

outf = paste(res_dir, "supp_fig1e.subsets_motif_frac_barplot.pdf", sep="")
pdf(file=outf,width=5,height=4)
plot_motif_frac(dfrac)
dev.off()


