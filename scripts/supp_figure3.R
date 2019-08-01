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

plot_venneuler <- function (dplot, mytitle="",selcolors) {
	colnames(dplot) = c("elements", "sets")
	v <- venneuler(dplot)
	plot(v, main=mytitle, col= selcolors)
}

print_pie_gr <- function (dd, myset) {
	
	dpie = dd %>% dplyr::select(mergedPeakId,genomic_region, Distance.to.TSS) %>% distinct()
	table(dpie$genomic_region)
	dpie[is.na(dpie$Distance.to.TSS), "Distance.to.TSS"] = 100000
	dpie[(dpie$Distance.to.TSS >= -2000 & dpie$Distance.to.TSS <=500), "genomic_region" ] = "promoter-TSS"
	dpie[(dpie$Distance.to.TSS > -10000 & dpie$Distance.to.TSS < -2000), "genomic_region" ] = "upstream-10K"
	dpie[(dpie$genomic_region %in% c("exon","3'UTR","5'UTR","non-coding","TTS")), "genomic_region" ] = "gene_body"
	dpie[(dpie$genomic_region %in% c("intron")), "genomic_region" ] = "intron"
	dpie[is.na(dpie$genomic_region) | dpie$genomic_region %in% "Intergenic","genomic_region"] = "integenic"
	table(dpie$genomic_region)
	
	set.seed(8)
	piecolors = c("deepskyblue","black",  "lemonchiffon1", "deeppink2","darkolivegreen3" )
		
	dplot = dpie %>% dplyr::select(mergedPeakId, genomic_region) %>% distinct() %>% group_by(genomic_region) %>% 
			summarise(npeaks=n_distinct(mergedPeakId)) %>% filter(!is.na(genomic_region)) %>%
			mutate(frac=npeaks/sum(npeaks)*100) %>% mutate(mylab=paste(genomic_region, ":", npeaks, "(",signif(frac, digits = 0), "%)",sep=""))
	
	pie(dplot$frac,main=paste(myset, "(", nrow(dd), " peaks)", sep=""),labels=dplot$mylab,col = piecolors,  border = NA) 
	
}

print_pie_motifs <- function (dd, myset) {
	
	dpie = dd %>% dplyr::select(mergedPeakId,motifs) %>% distinct()
	table(dpie$motifs)
	
	set.seed(8)
	piecolors = c("lightgray","red","dodgerblue", "orange")
	dplot = dpie %>% dplyr::select(mergedPeakId, motifs) %>% distinct() %>% group_by(motifs) %>% 
			summarise(npeaks=n_distinct(mergedPeakId)) %>% filter(!is.na(motifs))  %>%
			mutate(frac=npeaks/sum(npeaks)*100) %>% dplyr::mutate(mylab=paste(motifs, ":", npeaks, "(",signif(frac, digits = 0), "%)",sep=""))
	
	pie(dplot$frac,main=paste(myset, "(", nrow(dd), " peaks)", sep=""),labels=dplot$mylab,col = piecolors,  border = NA) 
	
}

ovlp_insv <- function (dd, myset) {
	
	dovlp = dd %>% dplyr::filter(!is.na(wtInsv__insvInsv.peakid))  %>% nrow()
	dovlp = data.frame(myset=myset, tot = nrow(dd), novlp = dovlp, ovlpfrac=dovlp/nrow(dd), notovlpfrac=(nrow(dd)-dovlp)/nrow(dd))
	dovlp
}

plot_ovlp <- function (dovlp) {
	
	dovlp = dovlp %>% dplyr::select(-tot) %>% gather(xx, pct,-matches("myset|novlp"))
	
	qqq= ggplot(data =dovlp, aes(x = myset, y = pct,fill=xx)) + geom_bar(stat="identity",width=0.5)
	qqq = qqq + 
			scale_fill_manual(values = c("slategray2","navy"),name = "")+
			labs(x="",y="Fraction of peaks") +
			theme_classic()+
			theme(legend.position = "right",axis.text.x = element_text(angle =-45,hjust=0,vjust=1))
	
	print(qqq) 
	
}

###########################################################################


load(file=paste(data_dir,"chip_ElbaAntibody_supertable_with_RNAseq_correctedMotif.RData",sep="") )
# flist,dm_peak,
colnames(dm_peak) = c("peakId", "ss")


###########################################################################
##### supp figure 3B insv binding #########################################
###########################################################################

outf = paste(res_dir, "supp_fig3b.Insv_binding_venn.pdf", sep="")
pdf(file=outf,width=7,height=7)
psel = dm_peak %>% filter(grepl("insvInsv", ss))
plot_venneuler(dplot=psel, mytitle= paste("Insv", sep=""), selcolors=c("blue","green","red","orange"))	
dev.off()

###########################################################################
##### supp figure 3C-E Elba3 binding ######################################
###########################################################################

elba3_elba12dep = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid)  & is.na(elba1Elba3__elba3Elba3.peakid) & is.na(elba2Elba3__elba3Elba3.peakid)) 
elba3_elba12indep = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid)  & (!is.na(elba1Elba3__elba3Elba3.peakid) | !is.na(elba2Elba3__elba3Elba3.peakid))) 


pief = paste(res_dir, "supp_fig3cde.Elba3_binding.pdf",sep="")
pdf(file=pief,width=6,height=6)

#motif
print_pie_motifs(dd=elba3_elba12dep %>% mutate(motifs=motifsCombine), myset="elba3_elba12dep_motifs")
print_pie_motifs(dd=elba3_elba12indep %>% mutate(motifs=motifsCombine), myset="elba3_elba12indep_motifs")

#gloc
print_pie_gr(dd=elba3_elba12dep, myset="elba3_elba12dep_gr")
print_pie_gr(dd=elba3_elba12indep, myset="elba3_elba12indep_gr")


#ovlp insv
dovlp = rbind(ovlp_insv(dd=elba3_elba12dep, myset="elba3_elba12dep"), 
		      ovlp_insv(dd=elba3_elba12indep, myset="elba3_elba12indep"))
plot_ovlp(dovlp)

dev.off()

system(paste("pdfjam ",pief,"  --frame false --nup 2x3 --suffix 2up --outfile ",pief,sep=""))


###########################################################################
##### supp figure 3 F-I Elba1 binding #####################################
###########################################################################

elba1_wt = Dsuper2 %>% filter((!is.na(wtElba1__elba1Elba1.peakid)))
elba1_elba2mut = Dsuper2 %>% filter((!is.na(elba2Elba1__elba1Elba1.peakid)) )
elba3_elba2mut = Dsuper2 %>% filter((!is.na(elba2Elba3__elba3Elba3.peakid)))

pief = paste(res_dir, "supp_fig3fgh.Elba1_binding.pdf",sep="")
pdf(file=pief,width=6,height=6)

# overlap venn
psel = dm_peak %>% filter(grepl("wtElba1__elba1Elba1|elba2Elba1__elba1Elba1|elba2Elba3__elba3Elba3", ss))
plot_venneuler(dplot=psel, mytitle= paste("Elba1", sep=""), selcolors = c("black","yellow", "red"))

#motif
print_pie_motifs(dd=elba1_wt %>% mutate(motifs=motifsCombine), myset="elba1_wt")
print_pie_motifs(dd=elba1_elba2mut %>% mutate(motifs=motifsCombine), myset="elba1_elba2mut")

#gloc
print_pie_gr(dd=elba1_wt, myset="elba1_wt")
print_pie_gr(dd=elba1_elba2mut, myset="elba1_elba2mut")

#ovlp insv
dovlp = rbind(ovlp_insv(dd=elba1_wt, myset="elba1_wt"), 
		      ovlp_insv(dd=elba1_elba2mut, myset="elba1_elba2mut"), 
		      ovlp_insv(dd=elba3_elba2mut, myset="elba3_elba2mut"))
plot_ovlp(dovlp)

dev.off()
system(paste("pdfjam ",pief,"  --frame false --nup 2x3 --suffix 2up --outfile ",pief,sep=""))


