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

plot_insulator_motif_frac <- function (dovlp) {

	dovlp= dovlp %>% mutate(ccc = paste(insulator, "_", xx, sep=""))
	
	dovlp$ss= factor(dovlp$ss, levels = unique(dovlp$ss))
	dovlp$insulator = factor(dovlp$insulator, levels = unique(dovlp$insulator))
	dovlp$ccc= factor(dovlp$ccc, levels = unique(dovlp$ccc))
	
	
	qqq= ggplot(data =dovlp, aes(x = insulator, y = frac,fill=ccc)) +  geom_bar(stat="identity",width=0.9) #geom_area(stat="identity") #
	qqq = qqq + 
			scale_fill_manual(values = c( "lightgray","darkorchid3", "lightgray", "thistle3", "lightgray","coral2"),name = "")+
			labs(x="",y="Fraction of peaks") +
			facet_wrap(~ ss,scales = "free_x",ncol=5) +
			theme_bw()+
			ggtitle("") + 
			theme(legend.position = "none",axis.text.x = element_text(angle =-90,hjust=0,vjust=1))

	
	print(qqq) 		  
}

get_insulator_motif_frac <- function (dd, myset) {

	dd = dd %>% dplyr::select(mergedPeakId, CP190, GAGA, BEAF32) %>% gather(insulator,xx, -mergedPeakId) %>% mutate(xx=ifelse(xx %in% "", 0, 1)) 
	dnn  = dd %>% group_by(insulator, xx) %>% summarise(nn=n_distinct(mergedPeakId)) %>% data.frame()
	dtot = dd %>%  group_by(insulator) %>% summarise(ntot=n_distinct(mergedPeakId))
	dplot = dnn %>% inner_join(dtot) %>% mutate(frac=nn/ntot) %>% arrange(insulator, xx) %>% mutate(ss = myset) 
	dplot2 = rbind(dplot[dplot$insulator %in% "CP190",],  dplot[dplot$insulator %in% "GAGA",],  dplot[dplot$insulator %in% "BEAF32",] )
	
	dplot2 
}


###########################################################################
##### supp figure 6A subset overlap insulator motifs #####################
###########################################################################

load(file =paste(data_dir,"chip_elba_supertable_subsets.RData",sep=""))
# chip_subsets, 
chip_subsets = chip_subsets[c("elba123_ovlp.ovlp_insv", "elba123_ovlp.notovlp_insv", "elba3insv_noelba12", "insv_only", "elba3_only")]
names(chip_subsets) = c("Elba123_Insv","Elba123_noInsv","Elba3Insv_noElba12", "Insv-unique", "Elba3-unique")


dfrac = NULL
for (ii in names(chip_subsets)) {
	dfrac = dfrac %>% bind_rows(get_insulator_motif_frac(dd=chip_subsets[[ii]], myset=ii))
}

outf = paste(res_dir, "supp_fig6a.subsets_insulator_motif_frac.pdf", sep="")
pdf(file=outf,width=7,height=3.7)
plot_insulator_motif_frac(dfrac)
dev.off()
