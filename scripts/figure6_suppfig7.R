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


########################################################################################
######################################## functions #####################################
########################################################################################

plot_boxplot <- function (dplot, mytitle, ylims) {
	
	p = ggplot(dplot, aes(x=paircate, y=pro.promoter_absFC, fill=ss)) + geom_boxplot() 
	
	if (!is.null(ylims)) {
		p = p + coord_cartesian(ylim = ylims)
	}
	
	p = p + 	labs(x="", y = "promoter FC") +
		    scale_fill_manual(values =c("red","green","orange","gray", "white"),name = "")+
			ggtitle(paste(mytitle," promoter",sep="")) +
			theme_classic()
	print(p)
	
	p = ggplot(dplot, aes(x=paircate, y=pro.body_absFC, fill=ss)) + geom_boxplot() 
	
	if (!is.null(ylims)) {
		p = p + coord_cartesian(ylim = ylims)
	}
	
	p = p + 
			labs(x="", y = "body FC") +
			scale_fill_manual(values =c("red","green","orange","gray", "white"),name = "")+
			ggtitle(paste(mytitle," body",sep="")) +
			theme_classic()
	print(p)
}

plot_mean_errorbar <- function (dplot, ylab="",xlab="", mytitle="",ylims=NULL){
	
	dplot$ss= factor(dplot$ss, levels = unique(dplot$ss))	
	qqq= ggplot(data = dplot, aes(x = ss, y = log2FC.mean, fill=ss))  + geom_bar(stat="identity",width=0.6)
	
	if (!is.null(ylims)) {
		qqq = qqq + coord_cartesian(ylim = ylims)
	}
	
	qqq = qqq + 
			scale_fill_manual(values =c("red","green","orange","gray", "white"),name = "") +
			labs(x="",y=ylab) +
			geom_errorbar(aes(ymin=log2FC_se_plus, ymax=log2FC_se_minus), width=0.3) +
			facet_wrap(~ paircate,scales = "fixed",ncol=3) +
			geom_hline(yintercept = 0,linetype=1,colour= I("gray50"),size=0.2)  +
			theme_classic()+
			ggtitle(paste(mytitle," promoter meanFC",sep="")) + 
			theme(legend.position = "right",axis.text.x = element_text(angle =-90,hjust=0,vjust=1))
	
	print(qqq) 
	
}

plot_peak_ratio <- function (dplot, mytitle) {
	
	dplot$genotype = factor(dplot$genotype,level=  unique(dplot$genotype))
	mycolor=c("red","green", "orange", "darkgray","black")
	qqq= ggplot(data =dplot, aes(x = genotype, y = FC,fill=genotype ,colour=I("black"))) + geom_bar(stat="identity",width=0.6)
	qqq = qqq + 
			scale_fill_manual(values = mycolor,name = "")+
			labs(x="",y="ratio of peaks") +
			theme_classic()+
			ggtitle(mytitle) + 
			theme(legend.position = "right",axis.text.x = element_text(angle =-45,hjust=0,vjust=1))
	
	print(qqq) 
}

get_propeak_ratio <- function (peaks_sel) {

	dpeak_sel = dPI %>% filter(PeakID %in% peaks_sel) %>% 
			dplyr::select(PeakID, genotype,  promoter_log2TPM) %>%
			spread(PeakID,promoter_log2TPM) 
	colnames(dpeak_sel) = c("genotype","peak1","peak2")
	dpeak_sel = dpeak_sel %>% mutate(FC =2^abs(peak1-peak2))
}

########################################################################################
########################### figure 6A,B and supp_fig7A-E ###############################
########################################################################################

load(file=paste(data_dir, "merge_chipseq_chipnexus_q10__pro_round2.RData",sep=""))
dpro = dist_peak %>% dplyr::select(mergedPeakId,upstreamGene.genotype, paircate,upstreamGene.Distance, downstreamGene.Distance,  matches("absFC")) %>% 
		mutate(distance=upstreamGene.Distance+downstreamGene.Distance)%>% distinct()

for (whichside in c("greater", "less")) {
	
	pdff = paste(res_dir, "fig6ab.supp_fig7.proseq_divergent_neighboring_absFC_boxplot.",whichside,".pdf",sep="")
	pdf(file=pdff,width=8,height=7,useDingbats=FALSE)
	sink(file =paste(res_dir, "fig6ab.supp_fig7.proseq_divergent_neighboring_absFC_boxplot.",whichside,".txt",sep="") )
	
	for (FC_cutoff in c(log2(4), log2(1))) {
		
		dist = 10000
		if (whichside == "greater") {
			dpro_sel = dpro %>% dplyr::filter(upstreamGene.genotype %in% "yw_pro") %>% filter(pro.promoter_absFC >= FC_cutoff)
			mytitle=paste("FC>",2^FC_cutoff,"-fold", sep="")
		}
		
		if (whichside == "less") {
			dpro_sel = dpro %>% dplyr::filter(upstreamGene.genotype %in% "yw_pro") %>% filter(pro.promoter_absFC <= FC_cutoff)	
			mytitle=paste("FC<",2^FC_cutoff,"-fold", sep="")
		}
		
		dpro_high = dpro %>% filter(mergedPeakId %in% dpro_sel$mergedPeakId)
		
		Dplot_box = dpro_high %>% dplyr::filter(distance <= dist) %>%  rename(ss=upstreamGene.genotype) %>% arrange(mergedPeakId, ss)  %>%
				group_by(ss, paircate)  %>% mutate(nloci = n_distinct(mergedPeakId)) %>% ungroup  %>% 
				mutate(paircate = paste(paircate,"(",nloci,")",sep="")) %>% dplyr::select(-nloci) %>% data.frame
		
		plot_boxplot(Dplot_box, mytitle, ylims=c(-0.01, 12.5))		
		
		
		cat("\n\nFC_cutoff>",2^FC_cutoff, " distance<", dist,"\n")
		for (jj in unique(Dplot_box$paircate)) {
			Dplot_ij = Dplot_box %>% filter(paircate %in% jj)
			
			for (ff in c("elba1_pro","elba2_pro", "elba3_pro", "insv_pro")) {
				cat(jj,": ",ff, " -- yw promoter : ", t.test(Dplot_ij[Dplot_ij$ss %in% ff,"pro.promoter_absFC"], Dplot_ij[Dplot_ij$ss %in% "yw_pro","pro.promoter_absFC"],alternative ="less")$p.value*12,"\n")			
				cat(jj,": ",ff,"  -- yw body : ", t.test(Dplot_ij[Dplot_ij$ss %in% ff,"pro.body_absFC"], Dplot_ij[Dplot_ij$ss %in% "yw_pro","pro.body_absFC"],alternative ="less")$p.value*12,"\n")	
			}
		}
		
		if (whichside == "greater") {
			Dplot_bar = dpro_high %>%  rename(ss=upstreamGene.genotype)%>% distinct() %>% dplyr::filter(distance <= dist) %>% 
					dplyr::select(-pro.body_absFC, -pro.PI_absFC, -upstreamGene.Distance, -downstreamGene.Distance)  %>% 
					spread(ss,pro.promoter_absFC) %>% 
					mutate(elba1_yw=elba1_pro-yw_pro,elba2_yw=elba2_pro-yw_pro,elba3_yw=elba3_pro-yw_pro,insv_yw=insv_pro-yw_pro) %>% 
					dplyr::select(-matches("_pro")) %>% gather(ss, log2FC, -c(mergedPeakId,paircate, distance))    %>%  
					group_by(ss,paircate) %>%
					summarise(nn=n_distinct(mergedPeakId), log2FC.mean=mean(log2FC,na.rm = TRUE), log2FC.sd=sd(log2FC,na.rm = TRUE)) %>% 
					ungroup() %>%   mutate(log2FC.se =  log2FC.sd/sqrt(nn)) %>% 
					mutate(paircate = paste(paircate, "(", nn,")",sep=""))  %>% 
					mutate(log2FC_se_plus = log2FC.mean+log2FC.se, log2FC_se_minus = log2FC.mean-log2FC.se)  %>%   data.frame()		
			plot_mean_errorbar(dplot = Dplot_bar, ylab="log2FC (mean+-se)", xlab="", mytitle) 
		}
	
	}
	sink()
	dev.off()
}

########################################################################################
################################# figure 6C-E, 7D ######################################
########################################################################################


load(file=paste(data_dir,"/proseq_peaks/combined_pro__ProseqE_tssFold4_bodyFold3_annot_round2.RData",sep=""))
load(file=paste(data_dir,"combined_pro__ProseqE_tssFold4_bodyFold3_annot_round2.RData",sep=""))

# dpeak,dPI, dPI_promoter ,dPI_coding, dPI_noncoding, 
dPI = dPI %>% dplyr::filter(grepl("AlignToDmOnly",genotype)) %>% mutate(genotype=gsub("_sorted.AlignToDmOnly", "", genotype))


pdff = paste(res_dir, "fig6cde_fig7d.PROseq_peakratio_selected.pdf",sep="")
pdf(file=pdff,width=3, height=6)

dpeak_sel = get_propeak_ratio( peaks_sel= c("chr2L-2072-1", "chr2L-916-0"))
plot_peak_ratio(dplot=dpeak_sel, mytitle = "Wnt4:Wg")

dpeak_sel = get_propeak_ratio( peaks_sel= c("chr2L-846-1","chr2L-845-1"))
plot_peak_ratio(dplot=dpeak_sel, mytitle = "CG31809:CG6012")

dpeak_sel = get_propeak_ratio( peaks_sel= c("chr3R-997-1","chr3R-2550-0"))
plot_peak_ratio(dplot=dpeak_sel, mytitle = "CG31145:GILT3")

dpeak_sel = get_propeak_ratio( peaks_sel= c("chr3R-997-1","chr3R-2550-0"))
plot_peak_ratio(dplot=dpeak_sel, mytitle = "CG31145:GILT3")

dpeak_sel = get_propeak_ratio( peaks_sel= c("chr2R-2414-0","chr2R-387-1"))
plot_peak_ratio(dplot=dpeak_sel, mytitle = "CG11275:CG11170")

dpeak_sel = get_propeak_ratio( peaks_sel= c("chr2R-2414-0","chr2R-2415-0"))
plot_peak_ratio(dplot=dpeak_sel, mytitle = "CG11275:MED16")

dpeak_sel = get_propeak_ratio( peaks_sel= c("chr2R-2414-0","chr2R-2416-0"))
plot_peak_ratio(dplot=dpeak_sel, mytitle = "CG11275:Vps35")

dev.off()

system(paste("pdfjam ",pdff," --frame false --nup 3x3 --suffix 3up --outfile ",pdff,sep=""))


