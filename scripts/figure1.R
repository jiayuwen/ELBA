#! /bin/env Rscript

library(Biostrings)
library(ggplot2)
library(tidyverse)
library(openxlsx)
library(venneuler)
library(getopt)

current_dir = "/home/jwen/projects/qi/elba/ELBA/"
data_dir  = paste(current_dir, "data/",sep="")
res_dir  = paste(current_dir, "res/",sep="")
coverage_dir  = paste(data_dir, "coverage/",sep="")

##########################################################################
########################### functions ####################################
##########################################################################

plot_venneuler <- function (dplot, mytitle="") {
	colnames(dplot) = c("elements", "sets")
	v <- venneuler(dplot)
	plot(v, main=mytitle, col= c("blue","green", "red",  "gray"))
}

get_chipseq_peak_set <- function (dd,scorecutoff=0) {
	
	elba1 = dd %>% filter(!is.na(wtElba1__elba1Elba1.peakid) & wtElba1__elba1Elba1.maxPeakscore>=scorecutoff) %>%  mutate(ss=paste("Elba1(",n_distinct(mergedPeakId), ")",sep=""))
	elba2 = dd %>% filter(!is.na(wtElba2__elba2Elba2.peakid) & wtElba2__elba2Elba2.maxPeakscore>=scorecutoff ) %>%  mutate(ss=paste("Elba2(",n_distinct(mergedPeakId), ")",sep=""))
	elba3 = dd %>% filter(!is.na(wtElba3__elba3Elba3.peakid) & wtElba3__elba3Elba3.maxPeakscore>=scorecutoff) %>%  mutate(ss=paste("Elba3(",n_distinct(mergedPeakId), ")",sep=""))
	insv = dd %>% filter(!is.na(wtInsv__insvInsv.peakid) & wtInsv__insvInsv.maxPeakscore>=scorecutoff) %>%  mutate(ss=paste("Insv(",n_distinct(mergedPeakId), ")",sep=""))

	dpeak = elba1 %>% bind_rows(elba2) %>% bind_rows(elba3) %>% bind_rows(insv) %>% dplyr::select(mergedPeakId, ss)
	dpeak
}

get_motif_ntiles <-function (ddx, nm, nbins=5) {	
	ddx = ddx %>% rename(xx=nm) %>% 
			dplyr::select(PeakID, xx,insv_motif) %>% arrange(-xx) %>% 
			filter(xx>0) %>% mutate(sntile= ntile(xx,nbins)) %>% 
			group_by(sntile) %>% summarise(nmotif = sum(insv_motif), nn=n()) %>% arrange(-sntile) %>% mutate(motif_frac=nmotif/nn)  %>% 
			mutate(ss = gsub("_peakscore", "", nm))  %>% mutate(sntile=1:nbins)
	ddx
}

call_spatial_heatmap  <- function (bwfiles,bedf="", averageTypeBins="max", sortRegions="descend", sortUsingSamples="None", pngf, pngheight=20,pngwidth=5, zmin="None", zmax="None", WW_up=2000, WW_down=2000, do_colorMap=T, mycolor="YlGnBu", rlabel="genes",ylabel="chip density", samplesLabel="", startLabel="TSS", heatmaponly=F, plotFileFormat="png") {
	
	matrixf = gsub(".png|.pdf", ".matrix.mat.gz", pngf, perl=T)
	
	if (WW_up>2000) {
		mystep= 50
	} else if (WW_up>1000) {
		mystep= 25
	} else {
		mystep= 10
	}
	
	if (!file.exists(matrixf)) {
		if (bedf == ""){
			system(paste("computeMatrix reference-point -S ", bwfiles," -R ",dm_gene_bed," -a ",WW_down," -b ",WW_up," -o ",matrixf," -bs ",mystep," -p 6 --skipZeros --averageTypeBins ",averageTypeBins, sep="")) #--minThreshold 1  
			
		} else {
			system(paste("computeMatrix reference-point -S ", bwfiles," -R ",bedf," -a ",WW_down," -b ",WW_up," -o ",matrixf," -bs ",mystep," -p 6 --skipZeros --averageTypeBins ",averageTypeBins, sep="")) #--minThreshold 1  
		}
	}
	
		if (heatmaponly) {
			system(paste("plotHeatmap -m ",matrixf ," -out ",pngf," --colorList ",mycolor," --sortRegions ",sortRegions, " --refPointLabel ", startLabel," --regionsLabel ",rlabel," --yAxisLabel ",  ylabel," --samplesLabel ",samplesLabel, " --heatmapHeight ",pngheight," --heatmapWidth ", pngwidth,"  --sortUsing max --xAxisLabel ''  --whatToShow 'heatmap and colorbar' --plotFileFormat ",plotFileFormat,sep=""))
		} else {
			system(paste("plotHeatmap -m ",matrixf ," -out ",pngf," --colorList ",mycolor," --sortRegions ",sortRegions, " --sortUsingSamples ", sortUsingSamples, " --refPointLabel ", startLabel," --regionsLabel ",rlabel," --yAxisLabel ",  ylabel," --samplesLabel ",samplesLabel, " --heatmapHeight ",pngheight," --heatmapWidth ", pngwidth,"  --zMin ",zmin," --zMax ",zmax," --sortUsing max --xAxisLabel '' --plotFileFormat ",plotFileFormat,sep=""))	
		}

}

########################################################################
##### figure 1.B coverage heatmap of 4 factors around peak summits #####
########################################################################
dm_gene_bed = paste("gene2.bed",sep="")
peak_bed = paste("peakRankByElba3_gene.bed", sep="")

setwd(coverage_dir)

bedGraph_files = list.files(path=coverage_dir,pattern = "wt(.*)log2ratio.bw$", full=F)
bedGraph_files = bedGraph_files[grepl("^wt", basename(bedGraph_files))]
for (bedGraphf in bedGraph_files) {
	if (grepl("Elba2",bedGraphf )) {
		mycolors = "white,white,aliceblue,blue,navy"
	} else {
		mycolors = "white,aliceblue,blue,navy"
	}
	pdff =  paste("diffdensity_",gsub("__log2ratio.bw", "", bedGraphf),"_rankByElba3_log2ratio.up3000.dn3000.pdf", sep="")
	call_spatial_heatmap(bwfiles=bedGraphf, bedf=peak_bed, averageTypeBins="max", sortRegions="no", pngheight=20,pngwidth=3, pngf =pdff,  WW_up=3000, WW_down=3000, do_colorMap=F, mycolor=mycolors, samplesLabel=gsub(".bw","",bedGraphf,perl=T),rlabel="gene",ylabel="chip_enrichment",heatmaponly=T,plotFileFormat="pdf")
}

########################################################################
##### figure 1.C 4 factor overlaps #####################################
########################################################################
load(file=paste(data_dir, "chip_supertable_with_RNAseq_correctedMotif.RData",sep=""))

# obtain the peak sets for 4 factors
dpeak_all= get_chipseq_peak_set(Dsuper,scorecutoff=0)

# plot overlpping venn diagram for 4 factors
outf = paste(res_dir, "fig1c.Elba4factor_intersection_venn.pdf", sep="")
pdf(file=outf,width=7,height=7)
plot_venneuler(dpeak_all, mytitle= paste("chip-seq 4 factor peaks", sep=""))
dev.off()


########################################################################
##### figure 1.D 4 peak score and motif fraction correlation ###########
########################################################################

Dmotif = Dsuper %>% dplyr::select(mergedPeakId, insv_CCAATTGG,insv_CCAATAAG, matches("maxPeakscore")) %>% 
		mutate(insv_motif = ifelse(!insv_CCAATTGG %in% "" | !insv_CCAATAAG %in% "", 1, 0)) %>% distinct()
colnames(Dmotif) = c("PeakID", "insv_CCAATTGG", "insv_CCAATAAG",  "elba1_peakscore", "elba2_peakscore", "elba3_peakscore", "insv_peakscore", "insv_motif")
Dmotif[is.na(Dmotif)] = 0

nbins = 10
d1 = get_motif_ntiles(Dmotif, nm = "elba1_peakscore", nbins=nbins)
d2 = get_motif_ntiles(Dmotif, nm = "elba2_peakscore", nbins=nbins)
d3 = get_motif_ntiles(Dmotif, nm = "elba3_peakscore", nbins=nbins)
d4 = get_motif_ntiles(Dmotif, nm = "insv_peakscore", nbins=nbins)
dd = d1 %>% bind_rows(d2) %>% bind_rows(d3) %>% bind_rows(d4)


outf = paste(res_dir, "fig1d.Elba4factor_motifs_peakscore_correlation_nbins",nbins,".pdf", sep="")
pdf(file=outf,width=6,height=5, useDingbats=FALSE)
mycolor = c("blue","darkgreen","deeppink", "darkgray", "brown", "pink","orange","hotpink","seagreen3","yellow", "lightblue","black","gray")
qqq= ggplot(data = dd, aes(x = sntile, y = motif_frac,colour=ss)) + geom_line()+geom_point()

qqq = qqq + scale_colour_manual(values = mycolor,name = "motif")+
		labs(x="peak score high to low bins",y="Insv/Elba motif fraction") +
		scale_x_continuous(breaks=1:nbins) +
		theme_classic()+
		ggtitle("") + 
		theme(legend.position = "right")
print(qqq) 
dev.off()


########################################################################
##### figure 1.E anchor at TSS, read coverage for genes with 
##### peaks containing Elba/Insv motif  
########################################################################
