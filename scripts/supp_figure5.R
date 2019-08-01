#! /bin/env Rscript


library(Biostrings)
library(ggplot2)
library(tidyverse)
library(magrittr)
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

plot_venneuler <- function (dplot, mytitle="") {
	colnames(dplot) = c("elements", "sets")
	v <- venneuler(dplot)
	plot(v, main=mytitle, col= c("orange","green", "skyblue",  "purple"))
}


plot_mean_errorbar <- function (dplot) {
	
	#  dplot = DDD
	dplot =  dplot %>% group_by(num)  %>% mutate(nloci = n_distinct(gene)) %>% ungroup  %>% 
			mutate(ss = paste(ss, "(",nloci,")",sep="")) %>% dplyr::select(-nloci) %>% data.frame
	
	dplot = dplot   %>% 
			group_by(num, ss, ee)  %>% 		
			summarise(FC.mean = mean(logFC), n=n(), FC.sd = sd(logFC))   %>%
			mutate(FC.se = FC.sd/sqrt(n))   %>%
			mutate(se_plus = FC.mean-FC.se, se_minus=FC.mean+FC.se) 

	dplot$ee= factor(dplot$ee, levels = unique(dplot$ee))
	dplot$ss= factor(dplot$ss, levels = unique(dplot$ss))
	
	
	qqq= ggplot(data = dplot, aes(x = ss, y = FC.mean,fill=ss)) + geom_bar(stat="identity",width=0.6)
	
	qqq = qqq + scale_fill_manual(values = c("magenta","blue","burlywood3", "green","red"),name = "")+
			labs(x="",y="log2FC (mutant/yw)") +
			geom_errorbar(aes(ymin=se_plus, ymax=se_minus), width=0.4) +
#				geom_bar(width=0.5) +  
			facet_wrap(~ee,scales = "fixed",ncol=5) +
			geom_hline(yintercept = 0,linetype=1,colour= I("gray50"),size=0.2)  +
			theme_bw()+
			theme(legend.position = "right",axis.text.x =element_blank(),  axis.ticks.x=element_blank())
	
	print(qqq) 
	
	
}


do_gesa_peakcall <- function (v_region, GG,  pp, cont.matrix, mycolor) {
	
	C2t <- ids2indices(GG, rownames(v_region)) 
	
	rr_camera <- camera(as.matrix(v_region),C2t,design=design,contrast=cont.matrix[,pp,drop=F],inter.gene.cor=0.01)
	rr_camera = rr_camera %>% mutate( geneset=rownames(rr_camera), pp=pp) #%>%filter(PValue <= 0.05)
	
	dxx =  res[[pp]]%>% rownames_to_column(var="GeneID")
	stat = dxx[order(match(dxx$GeneID, rownames(v_region))), "logFC"]
	for (i in names(C2t)) {
		barcodeplot(stat, C2t[[i]], col.bars=mycolor, main=paste("DF:",pp, " -- ",i," ",rr_camera[rr_camera$geneset == i, "Direction"], " Pvalue = ",signif(rr_camera[rr_camera$geneset == i, "PValue"], 2), sep=""))		
	}
	
	rr_camera
}



###########################################################################
##### supp figure 5.AB RNA-seq DF gene overlap ############################
###########################################################################


load(file=paste(data_dir, "chip_supertable_with_RNAseq_correctedMotif.RData",sep=""))
# Dsuper2
Dsuper2 = Dsuper2 %>%
		  dplyr::filter(abs(Distance.to.TSS) <= 2000)  %>% 
		  dplyr::filter((insv_CCAATTGG %in% "") | (insv_CCAATAAG %in% "") ) %>% 
		  dplyr::select(mergedPeakId, Gene.Name, matches("peakid") ) %>% 
		  distinct()

load(file=paste(data_dir, "elba__gene.rep.sig_final.RData", sep=""))
# sigGenes_list
for (gg in names(sigGenes_list)) {

	rnaseq = sigGenes_list[[gg]]
	rnaseq = rnaseq %>% mutate(elba1=1/3*(elba1_1+elba1_2+elba1_3),elba2=1/3*(elba2_1+elba2_2+elba2_3), 
			 elba3=1/3*(elba3_1+elba3_2+elba3_3),insv=1/3*(insv_1+insv_2+insv_3), yw=1/3*(yw_1+yw_2+yw_3) )
	if (gg == "elba1_vs_yw") {
		rnaseq = rnaseq %>% rename(elba1_yw.logFC = logFC, elba1_yw.adj.P.Val=adj.P.Val )  %>% 
				select(gene, yw, elba1, elba1_yw.logFC, elba1_yw.adj.P.Val)
	} else if (gg == "elba2_vs_yw") {
		rnaseq = rnaseq %>% rename(elba2_yw.logFC = logFC, elba2_yw.adj.P.Val=adj.P.Val )  %>% 
				select(gene, elba2, elba2_yw.logFC, elba2_yw.adj.P.Val)
	} else if (gg == "elba3_vs_yw") {
		rnaseq = rnaseq %>% rename(elba3_yw.logFC = logFC, elba3_yw.adj.P.Val=adj.P.Val )  %>% 
				select(gene, elba3, elba3_yw.logFC, elba3_yw.adj.P.Val)
	} else if (gg == "insv_vs_yw") {
		rnaseq = rnaseq %>% rename(insv_yw.logFC = logFC, insv_yw.adj.P.Val=adj.P.Val )  %>% 
				select(gene, insv, insv_yw.logFC, insv_yw.adj.P.Val)
	}
	Dsuper2 = left_join(Dsuper2, rnaseq, by=c("Gene.Name"="gene"))	
	
}

dfc_peaks = NULL
for (bn in c("Elba1", "Elba2", "Elba3", "Insv")) {
	
	if (grepl("Elba1",bn)) {
		ddi = Dsuper2 %>% filter(!is.na(wtElba1__elba1Elba1.peakid) )  %>%
				rename(logFC=elba1_yw.logFC, adjP=elba1_yw.adj.P.Val)  %>%
				dplyr::select(Gene.Name,logFC, adjP) %>% mutate(ss="elba1_yw")
		
	} else if  (grepl("Elba2",bn)) {
		ddi = Dsuper2 %>% filter(!is.na(wtElba2__elba2Elba2.peakid) )  %>%
				rename(logFC=elba2_yw.logFC, adjP=elba2_yw.adj.P.Val)  %>%
				dplyr::select(Gene.Name,logFC, adjP)  %>% mutate(ss="elba2_yw")
	} else if  (grepl("Elba3",bn)) {
		ddi = Dsuper2 %>% filter(!is.na(wtElba3__elba3Elba3.peakid) )  %>%
				rename(logFC=elba3_yw.logFC, adjP=elba3_yw.adj.P.Val)  %>%
				dplyr::select(Gene.Name,logFC, adjP)  %>% mutate(ss="elba3_yw")
	} else if  (grepl("Insv",bn)) {
		ddi = Dsuper2 %>% filter(!is.na(wtInsv__insvInsv.peakid) )  %>%
				rename(logFC=insv_yw.logFC, adjP=insv_yw.adj.P.Val)  %>%
				dplyr::select(Gene.Name,logFC, adjP) %>% mutate(ss="insv_yw")
	}
	dfc_peaks = rbind(dfc_peaks, ddi)
	
}
dfc_peaks = dfc_peaks[complete.cases(dfc_peaks), ]

Dsig_fdr = NULL
for (fdr_cutoff in c( 0.1, 0.2)) {
	for (fc_cutoff in c(1,1.3, 1.5)) {
		
		dsig_up = dfc_peaks %>%  filter(adjP <= fdr_cutoff & logFC >= log2(fc_cutoff)) %>%
				mutate(fdr_cutoff=fdr_cutoff, fc_cutoff=fc_cutoff, direction="up")  
		dsig_up = dsig_up %>% left_join(dsig_up %>% group_by(ss) %>% summarise(nn=n_distinct(Gene.Name))) %>%
				mutate(genotype=paste(ss,"(",nn,")",sep="") )   %>%
				dplyr::select(Gene.Name, genotype,fdr_cutoff,fc_cutoff,direction)
		
		dsig_down = dfc_peaks %>%  filter(adjP <= fdr_cutoff & logFC < log2(1/fc_cutoff)) %>%
				mutate(fdr_cutoff=fdr_cutoff, fc_cutoff=fc_cutoff, direction="down")  
		dsig_down = dsig_down %>% left_join(dsig_down %>% group_by(ss) %>% summarise(nn=n_distinct(Gene.Name))) %>%
				mutate(genotype=paste(ss,"(",nn,")",sep="") )   %>%
				dplyr::select(Gene.Name, genotype,fdr_cutoff,fc_cutoff,direction)
		
		Dsig_fdr = Dsig_fdr %>% bind_rows(dsig_up) %>% bind_rows(dsig_down)
	}
}	

pdff = paste(res_dir, "fig5a.RNAseq_chipseq_venn.pdf",sep="")
pdf(file=pdff,width=8, height=8)
for (ii in unique(Dsig_fdr$fdr_cutoff)) {
	for (jj in unique(Dsig_fdr$fc_cutoff)) {
		dsig = Dsig_fdr %>% filter(fdr_cutoff == ii & fc_cutoff == jj)
		dsig_up = dsig %>% filter( direction %in% "up") %>% dplyr::select(Gene.Name, genotype)
		dsig_down = dsig %>% filter( direction %in% "down") %>% dplyr::select(Gene.Name, genotype)
		
		plot_venneuler(dplot= dsig_up %>% dplyr::select(Gene.Name, genotype), mytitle=paste("Up-regulation: FDR <=",ii," & FC >", jj,"-fold",sep=""))
		plot_venneuler(dplot= dsig_down %>% dplyr::select(Gene.Name, genotype), mytitle=paste("Down-regulation: FDR <=",ii," & FC >", jj,"-fold",sep=""))
	}
}
dev.off()
system(paste("pdfjam ",pdff," --frame false --nup 2x2 --suffix 2up --outfile ",pdff,sep=""))


###########################################################################
##### supp figure 5.C RNA-seq fold change for genes with peak subsets #####
###########################################################################

load(file =paste(data_dir,"chip_elba_supertable_subsets.RData",sep=""))
# chip_subsets, 
chip_subsets = chip_subsets[c("elba123_ovlp.ovlp_insv", "elba123_ovlp.notovlp_insv", "elba3insv_noelba12", "insv_only", "elba3_only")]
names(chip_subsets) = c("Elba123_Insv","Elba123_noInsv","Elba3Insv_noElba12", "Insv-unique", "Elba3-unique")


# combine DF fold change from each factors
load(file=paste(data_dir, "elba__gene.rep.sig_final.RData",sep=""))
Dexpr = NULL
for (mm in names(sigGenes_list)) {
	sigGenes = sigGenes_list[[mm]]
	sigGenes = sigGenes %>%  dplyr::select(gene, logFC) %>% mutate(ee = mm)
	Dexpr = Dexpr %>% bind_rows(sigGenes)
}

DDD = NULL
cc = 0
for (ii in  names(chip_subsets)) {
	cc =  cc+1
	dcc = chip_subsets[[ii]] %>% dplyr::filter(!motifsCombine %in% "") %>%  dplyr::filter(abs(Distance.to.TSS) <= 2000)
	dii = Dexpr %>%  dplyr::filter(gene %in% dcc$Gene.Name ) %>% mutate(ss = ii, num=cc)
	DDD = DDD %>% bind_rows(dii)
}

# plot log2FC for each peak subset
pdf(file= paste(res_dir, "supp_fig5c.elba_subsets_RNAseq_FC.pdf",sep=""),width=8,height=3)
plot_mean_errorbar(dplot=DDD)
dev.off()


# Geneset enrichment analysis for each peak subset
diffpeak_sets = list()
for (ss in names(chip_subsets)) {
	diffpeak_sets[[paste(ss, ".motif", sep="")]] =  chip_subsets[[ss]] %>% filter(!motifsCombine %in% "") %>% filter(abs(Distance.to.TSS) <= 2000) %$% Gene.Name
}
names(diffpeak_sets)

pdff = paste(res_dir, "supp_fig5c.GESA_elba_subsets_ranks.pdf",sep="")
pdf(file = pdff, height=5, width=7)

res_list = list()
for (bn in c("elba1","elba2", "elba3", "insv" )) {
	
	pp = paste(bn, "_vs_yw", sep="")
	cat(pp, "\n")

	GG = diffpeak_sets
	lapply(GG, length)
	
	v_region = v$E
	res_list[[bn]] = do_gesa_peakcall(v_region, GG, pp, cont.matrix,mycolor="black")
	
}
write.xlsx(res_list, file =paste(res_dir, "supp_fig5c.GESA_elba_subsets_ranks.xlsx",sep=""), asTable = rep(TRUE,length((res_list)), rowNames=F, firstRow=TRUE, firstCol=TRUE, colWidths="auto"))

dev.off()

