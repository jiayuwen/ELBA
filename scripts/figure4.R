#! /bin/env Rscript


library(Biostrings)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(venneuler)
library(getopt)

current_dir = "/home/jwen/projects/qi/elba/ELBA/"
setwd(current_dir)

data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")

##########################################################################
########################### functions ####################################
##########################################################################

plot_mean_errorbar <- function (dplot) {

	dplot = dplot   %>% 
			group_by(ee, ss)  %>% 		
			summarise(FC.mean = mean(logFC), n=n(), FC.sd = sd(logFC))   %>%
			mutate(FC.se = FC.sd/sqrt(n))   %>%
			mutate(se_plus = FC.mean-FC.se, se_minus=FC.mean+FC.se)   
	
	dplot$ee= factor(dplot$ee, levels = unique(dplot$ee))
	dplot$ss= factor(dplot$ss, levels = unique(dplot$ss))	
	
	qqq= ggplot(data = dplot, aes(x = ss, y = FC.mean,fill=ss)) + geom_bar(stat="identity",width=0.6)
	
	qqq = qqq + scale_fill_manual(values = c("black", "green", "red"),name = "")+
			labs(x="",y="log2FC (mutant/wt)") +
			geom_errorbar(aes(ymin=se_plus, ymax=se_minus), width=0.4) + 
			facet_wrap(~ ee,scales = "free_x",ncol=4) +
			theme_bw()+
#			theme(legend.position = "right",axis.text.x = element_text(angle =-45,hjust=0,vjust=1))
			theme(legend.position = "right",axis.text.x =element_blank(),  axis.ticks.x=element_blank())
	
	print(qqq) 
	
}

###########################################################################
##### figure 4.B RNA-seq fold change for genes with peaksets ##############
###########################################################################
load( file = paste(file=paste(data_dir, "genesets2.RData",sep="")))
#diffpeak_sets
diffpeak_sets = diffpeak_sets[!grepl("tile|top100",names(diffpeak_sets))]
names(diffpeak_sets)

# combine DF fold change from each factors
load(file=paste(data_dir, "elba__gene.rep.sig_final.RData",sep=""))
dfc = NULL
for (ss in names(diffpeak_sets)) {
	if (grepl("Elba1", ss)) {
		ee = "elba1_vs_yw"	
	} else if (grepl("Elba2", ss)) {
		ee = "elba2_vs_yw"
	} else if (grepl("Elba3", ss)) {
		ee = "elba3_vs_yw"
	} else if (grepl("Insv", ss)) {
		ee = "insv_vs_yw"
	}
	
	dfc_ss = sigGenes_list[[ee]] %>% dplyr::select(gene, logFC) %>% 
			   dplyr::filter(gene %in% diffpeak_sets[[ss]]) %>% mutate(ee=ee, ss=gsub("(.*)__(.*)","\\2", ss)) 
	dfc = dfc %>% bind_rows(dfc_ss)
}

pdff = paste(res_dir, "fig4b.elba_peaksets_RNAseq_FC.pdf",sep="")
pdf(file=pdff,width=6,height=4)
plot_mean_errorbar(dplot=dfc)
dev.off()










