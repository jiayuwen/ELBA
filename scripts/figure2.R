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


########################################################################
#######################  functions  ####################################
########################################################################

plot_venneuler <- function (dplot, mytitle="", ccolours="") {
	colnames(dplot) = c("elements", "sets")
	v <- venneuler(dplot)
	plot(v, main=mytitle, col= ccolours)
}



########################################################################
##### figure 2.A & suppl figure 3.A #peaks for wt and mutant ###########
########################################################################

load(file=paste(data_dir,"chip_ElbaAntibody_supertable_with_RNAseq_correctedMotif.RData",sep="") )
colnames(dm_peak) = c("peakId", "ss")
dm_peak = dm_peak %>% group_by(ss) %>% summarise(npeaks = n_distinct(peakId)) %>% ungroup %>% 
		  separate(ss, into=c("genotype_antibody", "mutant")) %>% 
		  mutate(genotype=gsub("(.*)(Elba|Insv)(.*)", "\\1", genotype_antibody))  %>% 
		  mutate(antibody=gsub("(.*)(Elba|Insv)(.*)", "\\2\\3", genotype_antibody)) %>% 
		  dplyr::select(genotype, antibody,npeaks) %>% 
		  bind_rows(tibble(genotype="elba3",antibody="Elba1",npeaks=0))  
  

## figure 2.A
pdf(file=paste(res_dir, "fig2a.Elba_npeaks.pdf",sep=""), width=5, height=4)
dplot = dm_peak %>%  arrange(antibody, desc(genotype)) 
dplot$genotype = factor(dplot$genotype, levels = unique(dplot$genotype))
dplot$antibody = factor(dplot$antibody, levels = unique(dplot$antibody))

qqq= ggplot(data = dplot, aes(x = genotype, y = npeaks, fill=antibody)) + geom_bar(stat="identity",width=0.5)
qqq = qqq + 
		  labs(x="",y="Number of peaks") +
		  scale_fill_manual(values =c("cornflowerblue","chartreuse4","mediumorchid1","gray50"), name="" ) +
		  facet_wrap(~ antibody,scales = "free",ncol=4) +
		  theme_classic()+
		  theme(legend.position = "bottom",axis.text.x = element_text(angle =90,hjust=1,vjust=0))
  
print(qqq) 
dev.off()

## suppl figure 3.A
pdf(file=paste(res_dir, "supp_fig3a.Elba_npeaks.pdf",sep=""), width=6, height=4)
dplot = dm_peak  %>%  arrange(desc(antibody), desc(genotype)) 
dplot$genotype = factor(dplot$genotype, levels = unique(dplot$genotype))
dplot$antibody = factor(dplot$antibody, levels = unique(dplot$antibody))

qqq= ggplot(data = dplot, aes(x = antibody, y = npeaks, fill=antibody)) + geom_bar(stat="identity",width=0.5)
qqq = qqq + 
		labs(x="",y="Number of peaks") +
		scale_fill_manual(values =c("gray50","mediumorchid1","cornflowerblue","chartreuse4"), name="" ) +
		facet_wrap(~ genotype,scales = "free",ncol=5) +
		theme_classic()+
		theme(legend.position = "bottom",axis.text.x = element_text(angle =90,hjust=1,vjust=0))

print(qqq) 
dev.off()


########################################################################
##### figure 2.B peak overlaps of Elba3 binding ########################
########################################################################

wtElba3 = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid))  %>% mutate(ss=paste("wtElba3(",n_distinct(mergedPeakId), ")",sep=""))
elba2Elba3 = Dsuper2 %>% filter(!is.na(elba2Elba3__elba3Elba3.peakid))  %>% mutate(ss=paste("elba2Elba3(",n_distinct(mergedPeakId), ")",sep=""))
elba1Elba3 = Dsuper2 %>% filter(!is.na(elba1Elba3__elba3Elba3.peakid)) %>% mutate(ss=paste("elba1Elba3(",n_distinct(mergedPeakId), ")",sep=""))
dpeak = wtElba3 %>% bind_rows(elba2Elba3) %>% bind_rows(elba1Elba3) %>% dplyr::select(mergedPeakId, ss)

# plot overlpping venn diagram for 4 factors
outf = paste(res_dir, "fig2b.Elba3_venn.pdf", sep="")
pdf(file=outf,width=7,height=7)
plot_venneuler(dpeak, mytitle= paste("chip-seq Elba3 wt/mutant peaks", sep=""), ccolours=c("blue","red", "yellow"))
dev.off()


