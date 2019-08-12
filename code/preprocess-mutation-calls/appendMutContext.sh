#!/bin/sh

. ./code/utils/variables.sh

## WGS
$RSCRIPT code/preprocess-mutation-calls/appendContextToMutations.R Data/mut-files/rdas/mut-wgs-fdr5.bcsqFilter.rda Data/mut-files/triplet-seq-context/seqContext-wgs.fa Data/mut-files/triplet-seq-context/mut-wgs-fdr5.context.rda
$RSCRIPT code/preprocess-mutation-calls/appendContextToMutations.R Data/mut-files/rdas/mut-wgs-fdr5.bcsqFilter.wsharinfo.rda Data/mut-files/triplet-seq-context/seqContext-wgs.fa Data/mut-files/triplet-seq-context/mut-wgs-fdr5.context.wsharinfo.rda

## WES
$RSCRIPT code/preprocess-mutation-calls/appendContextToMutations.R Data/mut-files/rdas/mut-wes-fdr5.bcsqFilter.rda Data/mut-files/triplet-seq-context/seqContext-wes.fa Data/mut-files/triplet-seq-context/mut-wes-fdr5.context.rda
$RSCRIPT code/preprocess-mutation-calls/appendContextToMutations.R Data/mut-files/rdas/mut-wes-fdr5.bcsqFilter.wsharinfo.rda Data/mut-files/triplet-seq-context/seqContext-wes.fa Data/mut-files/triplet-seq-context/mut-wes-fdr5.context.wsharinfo.rda

## HWES
$RSCRIPT code/preprocess-mutation-calls/appendContextToMutations.R Data/mut-files/rdas/mut-hwes-fdr5.bcsqFilter.rda Data/mut-files/triplet-seq-context/seqContext-hwes.fa Data/mut-files/triplet-seq-context/mut-hwes-fdr5.context.rda
$RSCRIPT code/preprocess-mutation-calls/appendContextToMutations.R Data/mut-files/rdas/mut-hwes-fdr5.bcsqFilter.wsharinfo.rda Data/mut-files/triplet-seq-context/seqContext-hwes.fa Data/mut-files/triplet-seq-context/mut-hwes-fdr5.context.wsharinfo.rda
