#!/usr/bin/env bash

. ./code/utils/variables.sh

message "$0: Creating rds files with mutation sharing info in Data/mutShar"
$RSCRIPT code/mutShar/makeMutSharDataStruct.R Data/mut-files/rdas/mut-wes-fdr5.collapseDinucs.rda Data/mutShar/mut-wes-fdr5.collapseDinucs.mutShar.rds
$RSCRIPT code/mutShar/makeMutSharDataStruct.R Data/mut-files/rdas/mut-hwes-fdr5.collapseDinucs.rda Data/mutShar/mut-hwes-fdr5.collapseDinucs.mutShar.rds
$RSCRIPT code/mutShar/makeMutSharDataStruct.R Data/mut-files/rdas/mut-wgs-fdr5.collapseDinucs.rda Data/mutShar/mut-wgs-fdr5.collapseDinucs.mutShar.rds
message "$0: Creating text files with lists of donors that have two lines in Data/mutShar"
$RSCRIPT code/mutShar/makeListOfDonorsWithTwoLines.R Data/mutShar/mut-wes-fdr5.collapseDinucs.mutShar.rds Data/mutShar/wes-2donors.txt
$RSCRIPT code/mutShar/makeListOfDonorsWithTwoLines.R Data/mutShar/mut-hwes-fdr5.collapseDinucs.mutShar.rds Data/mutShar/hwes-2donors.txt
$RSCRIPT code/mutShar/makeListOfDonorsWithTwoLines.R Data/mutShar/mut-wgs-fdr5.collapseDinucs.mutShar.rds Data/mutShar/wgs-2donors.txt
