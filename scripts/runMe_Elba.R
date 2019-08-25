#! /bin/env Rscript


cd /Users/jean_macpro_2017/projects/Qi_chip/code/
rsync -avzpP -e 'ssh -o TCPKeepAlive=yes -o ServerAliveInterval=100' Elba_NatureComm kite2:/home/jwen/projects/qi/code/
		
########################################################################
########################################################################

code_dir = "/home/jwen/projects/qi/code/Elba_NatureComm/"

system(paste(code_dir,"chipseq_peakcalling_annotation.R",sep=""))

system(paste(code_dir,"chipnexus_peakcalling_annotation.R",sep=""))

system(paste(code_dir,"chipseq_merge_peaksets.R",sep=""))

system(paste(code_dir,"chipseq_chipnexus_merge_peaksets.R",sep=""))

system(paste(code_dir,"limma_elba_RNAseq.R",sep=""))

system(paste(code_dir,"PROseq_insulator_annotation.R",sep=""))
