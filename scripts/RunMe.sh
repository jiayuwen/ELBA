#!/bin/bash

## required tools: MACS2, HOMER2, memechip, fimo, bedtools, deeptools
## required R packages: tidyverse, ggplot2, Biostrings

cd ELBA/

## RNA-seq DF
./limma_elba_RNAseq.R

## ChIP-seq peak calling and annotation
./chipseq_peakcalling_annotation.R

## ChIP-nexus peak calling and annotation
./chipnexus_peakcalling_annotation.R

## combine chip-seq peaks from wt, mutant for 4 factors
./chipseq_merge_peaksets.R

## combine chip-nexus peaks for 4 factors
./chipseq_chipnexus_merge_peaksets.R

## PRO-seq peak annotation
./PROseq_peakcalling_annotation.R

## combine chip-nexus peaks for 4 factors
./PROseq_insulator_annotation.R


## Figures for the paper
./figure1.R

./figure2.R

./figure3.R

./figure4.R

./figure5.R

./figure6_suppfig7.R

./supp_figure1.R

./supp_figure2.R

./supp_figure3.R

./supp_figure4.R

./supp_figure5.R

./supp_figure6.R

