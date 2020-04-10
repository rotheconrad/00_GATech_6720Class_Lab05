#!/usr/bin/env Rscript

# @author  Roth Conrad
# @license Artistic-2.0
# Apologies for the quick and dirty code.

#= Load stuff
suppressPackageStartupMessages(library(enveomics.R))
args <- commandArgs(trailingOnly=TRUE)

infile=args[1]
load(infile)

rp <- enve.recplot2.changeCutoff(rp, 95)

mn=mean(enve.recplot2.seqdepth(rp)) # <- Average
md=median(enve.recplot2.seqdepth(rp)) # <- Median
tad80=enve.truncate(enve.recplot2.seqdepth(rp), 0.8) # <- 80% Central Truncated Mean
anir95=enve.recplot2.ANIr(rp, c(95,100)) # <- Reads above 95%

message(
	sprintf("\n\n"),
	sprintf("Mean Sequencing Depth (Coverage): %f\n", mn),
	sprintf("Median Sequencing Depth: %f\n", md),
	sprintf("80%% Truncated Average Sequcing Depth (TAD80): %f\n", tad80),
	sprintf("ANI of reads >= 95%% Sequence Identity (ANIr): %f", anir95),
	sprintf("\n\n")
	)

