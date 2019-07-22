## Attempt to collect all the relevant analyses for the first draft of the paper into a single pipeline
## This script should unpack everything and requires:
## 1. A raw mutation call directory (with results of Fisher Tests, and bcfs) located in Data/mut-files/raw
## 2. FDR threshold files (one per data set, but each using the same thresholds) located in Data/mut-files/fdr-thresholds
## 3. A CNV location file: Data/cnv/atab_20161018.tsv

## 1. Create dataset specific threshold files - these are now the same, but there's a separate file for historical reasons
ln -s fdrThresholds.txt Data/mut-files/fdr-thresholds/fdr-wes.txt
ln -s fdrThresholds.txt Data/mut-files/fdr-thresholds/fdr-hwes.txt
ln -s fdrThresholds.txt Data/mut-files/fdr-thresholds/fdr-wgs.txt

## 2. Preprocess call files
bsub -M6000 -q normal -J process.WES -o farmOut/process.WES.%J.stdout -e farmOut/process.WES.%J.stderr -R"select[mem>6000] rusage[mem=6000]" ./code/preprocess-mutation-calls/preprocess-call-files.sh Data/mut-files/raw/exomes.ft.txt.gz Data/mut-files wes 0
bsub -M6000 -q normal -J process.HWES -o farmOut/process.HWES.%J.stdout -e farmOut/process.HWES.%J.stderr -R"select[mem>6000] rusage[mem=6000]" ./code/preprocess-mutation-calls/preprocess-call-files.sh Data/mut-files/raw/high-vs-low-exomes.ft.txt.gz Data/mut-files hwes 0
bsub -M10000 -q normal -J process.WGS -o farmOut/process.WGS.%J.stdout -e farmOut/process.WGS.%J.stderr -R"select[mem>10000] rusage[mem=10000] span[hosts=1]" -n4 ./code/preprocess-mutation-calls/preprocess-call-files.sh Data/mut-files/raw/wgs.396.ft.txt.gz Data/mut-files wgs 0

## 3. Make the source R data structs
./code/preprocess-mutation-calls/make-all-data-structs.sh

## 4. Make figure 1 panels

## 3.1 Add the shared / line specific field R data structures
./code/masterAnalysis/masterAnalysis-v5/setUpStructuresForSharingAnalysis-v2.sh
./code/masterAnalysis/masterAnalysis-v5/submitFindSharingJobs-v2.sh
./code/masterAnalysis/masterAnalysis-v5/concatenateAndAddSharInfoToRds.sh

## 3.2 Add sequence context to mutation rda files (v4 now removes X chromosome mutations)
Rscript code/nmf/getCoordsOfPointMuts-v4.R Data/newCalls2/rdas/mut-wgs-fdr5.bcsqFilter.rda Data/nmf/wgs/coords.bed
Rscript code/nmf/getCoordsOfPointMuts-v4.R Data/newCalls2/rdas/mut-wes-fdr5.bcsqFilter.rda Data/nmf/wes/coords.bed
Rscript code/nmf/getCoordsOfPointMuts-v4.R Data/newCalls2/rdas/mut-hwes-fdr5.bcsqFilter.rda Data/nmf/hwes/coords.bed
twoBitToFa -noMask ~/local/ucsc/hg19/hg19.2bit -bed=Data/nmf/wgs/coords.bed Data/nmf/wgs/seqContext.fa && gzip -f Data/nmf/wgs/seqContext.fa
twoBitToFa -noMask ~/local/ucsc/hg19/hg19.2bit -bed=Data/nmf/wes/coords.bed Data/nmf/wes/seqContext.fa && gzip -f Data/nmf/wes/seqContext.fa
twoBitToFa -noMask ~/local/ucsc/hg19/hg19.2bit -bed=Data/nmf/hwes/coords.bed Data/nmf/hwes/seqContext.fa && gzip -f Data/nmf/hwes/seqContext.fa
./code/masterAnalysis/masterAnalysis-v5/appendMutContext-v2.sh

## 4. Do polyclonal analysis
Rscript code/polyclonal/fitMixModel.R Data/newCalls2/rdas/mut-wgs-fdr5.collapseDinucs.rda Data/polyclonal/mixModelRun.rds 

## 5. Discover mutation signatures
./code/masterAnalysis/masterAnalysis-v5/doMutSigs.sh

## 5.1 Get exposures to discovered signatures, and get differences between lines
./code/masterAnalysis/masterAnalysis-v5/getExposures.sh

## 5.2 Get exposures in shared and line specific (uniq) mutations
./code/masterAnalysis/masterAnalysis-v5/getLineSpSharedExposures-v2.sh

## 5. Make fitness impact files
BCFTOOLS=$HOME/local/src/bcftools-git/bcftools/bcftools
SRCDIR=/lustre/scratch116/vr/user/pd3/hipsci/exome-point-mutations/unified-set/final-set-2018-02-05.all
BCF=$SRCDIR/exomes.bcf
$BCFTOOLS query -f '%CHROM\t%POS\t%AF_1KG\t%AF_EXAC\t%CSQ\t%BCSQ\n' $BCF | ./code/fitnessImpact/parseCsq.pl stdin | gzip > Data/fitnessImpact/wes.all.sift-pphen.txt.gz
BCF=$SRCDIR/wgs.396.bcf
$BCFTOOLS query -f '%CHROM\t%POS\t%AF_1KG\t%AF_EXAC\t%CSQ\t%BCSQ\n' $BCF | ./code/fitnessImpact/parseCsq.pl stdin | awk '$5!="NA" || $7!="NA"' | gzip > Data/fitnessImpact/wgs.all.sift-pphen.txt.gz
BCF=$SRCDIR/high-vs-low-exomes.bcf
$BCFTOOLS query -f '%CHROM\t%POS\t%AF_1KG\t%AF_EXAC\t%CSQ\t%BCSQ\n' $BCF | ./code/fitnessImpact/parseCsq.pl stdin | gzip > Data/fitnessImpact/hwes.all.sift-pphen.txt.gz
./code/fitnessImpact/makeCosmicCancerMutList.sh Data/cancer/CosmicGenomeScreensMutantExport.tsv.gz | gzip > Data/fitnessImpact/cosmicMutList.bed.gz
Rscript code/fitnessImpact/makeCosmicMutDataStruct.R
Rscript code/fitnessImpact/makeVepCsqImpactDataStructs.R
Rscript code/fitnessImpact/computeSiftPphenLofSumStat-v2.R Data/fitnessImpact/wes.all.sift-pphen.rds Data/newCalls2/rdas/mut-wes-fdr5.bcsqFilter.rda Data/fitnessImpact/siftPhenLofSumStat-wes.rds
Rscript code/fitnessImpact/computeSiftPphenLofSumStat-v2.R Data/fitnessImpact/wgs.all.sift-pphen.rds Data/newCalls2/rdas/mut-wgs-fdr5.bcsqFilter.rda Data/fitnessImpact/siftPhenLofSumStat-wgs.rds
Rscript code/fitnessImpact/makeCosmicOverlapDataStruct-v2.R Data/fitnessImpact/cosmicOverlap.rds
Rscript code/fitnessImpact/doCosmicTissueFisherTest.R Data/fitnessImpact/cosmicOverlap.rds Data/fitnessImpact/cosmicTissueFisherTest.rds
Rscript code/fitnessImpact/doCosmicHistologyFisherTest.R Data/fitnessImpact/cosmicOverlap.rds Data/fitnessImpact/cosmicHistologyFisherTest.rds

## 6. Do DnDs C/V analysis
cp /nfs/team78pc20/im3/Projects/Selection_paper/Mutation_calls_to_share/Public_TCGA_calls_with_tissue.txt Data/dndscv/ ## Now missing, chase with Inigo
Rscript code/dndscv/makeTcgaStruct.R Data/dndscv/Public_TCGA_calls_with_tissue.txt Data/dndscv/tcga.rds
Rscript code/dndscv/getTcgaSampleSizes.R
./code/dndscv/doAllDnDsCancerGlobal.sh Data/dndscv/tcga Data/dndscv/tcga-dndsRes/ debug noclobber
Rscript code/dndscv/makeIpsCallTables-v4.R
./code/dndscv/getGermline.sh
Rscript code/dndscv/getGermline.R
Rscript code/dndscv/makeGermlineCallTables.R
./code/masterAnalysis/masterAnalysis-v5/runDnDsCv.sh

## 7. Mutation sharing
Rscript code/mutShar/makeMutSharDataStruct.R Data/newCalls2/rdas/mut-wes-fdr5.collapseDinucs.rda Data/mutShar/mut-wes-fdr5.collapseDinucs.mutShar.rds
Rscript code/mutShar/makeMutSharDataStruct.R Data/newCalls2/rdas/mut-wgs-fdr5.collapseDinucs.rda Data/mutShar/mut-wgs-fdr5.collapseDinucs.mutShar.rds
Rscript code/mutShar/makeSharedUniqueTab-v2.R
Rscript code/mutShar/countLineSpecificMaf0CCTT-v2.R Data/mutShar/mut-wgs-fdr5.lineSpCounts.rds Data/mutShar/lineSpecificMaf0CCTT.rds
Rscript code/mutShar/satelliteTest.R Data/mutShar/mut-wgs-fdr5.lineSpCounts.rds Data/mutShar/mut-wgs-fdr5.collapseDinucs.mutShar.rds Data/mutShar/satelliteTest.rds 
Rscript code/mutShar/makeLineSpSharedFreqGrid.R Data/mutShar/mut-wgs-fdr5.lineSpCounts.rds Data/mutShar/mut-wgs-fdr5.collapseDinucs.mutShar.rds Data/mutShar/afGrid.rds

## Rscript code/mutShar/countLineSpecificMaf0CCTT.R Data/newCalls2/rdas/mut-wgs-fdr5.bcsqFilter.wsharinfo.rda Data/mutShar/lineSpecificMaf0CCTT.rds

## Old

# ./code/newCalls/doAllMakeNewCallsDataStruct.sh Data/newCalls2 Data/newCalls2/rdas wes farm clobber
# ./code/newCalls/doAllMakeNewCallsDataStruct.sh Data/newCalls2 Data/newCalls2/rdas hwes farm clobber
# ./code/newCalls/doAllMakeNewCallsDataStruct.sh Data/newCalls2 Data/newCalls2/rdas wgs farm clobber

# Rscript code/nmf/makeSharedUniq-v2.R Data/nmf/wgs/mut-wgs-fdr5.context.wsharinfo.rda Data/nmf/wgs/mut-wgs-fdr5.context
# Rscript code/nmf/makeSharedUniq-v2.R Data/nmf/wes/mut-wes-fdr5.context.wsharinfo.rda Data/nmf/wes/mut-wes-fdr5.context
# Rscript code/nmf/makeMatrix-v2.R Data/nmf/wgs/mut-wgs-fdr5.context-shar.rds  Data/nmf/wgs/mut-wgs-fdr5-shar.mutMat.rds collapseOnDonor
# Rscript code/nmf/makeMatrix-v2.R Data/nmf/wgs/mut-wgs-fdr5.context-uniq.rds  Data/nmf/wgs/mut-wgs-fdr5-uniq.mutMat.rds
# Rscript code/sigfit/getExposure.R Data/nmf/wgs/mut-wgs-fdr5-uniq.mutMat.rds Data/sigfit/extractedSigs-wgs-expo+sigs.rds Data/sigfit/wgs-exposures-uniq.rds
# Rscript code/sigfit/getExposure.R Data/nmf/wgs/mut-wgs-fdr5-shar.mutMat.rds Data/sigfit/extractedSigs-wgs-expo+sigs.rds Data/sigfit/wgs-exposures-shar.rds


# Rscript code/sigfit/makeMutationMatrices-sharedUniq.R Data/nmf/wgs/mut-wgs-fdr5.context.rda Data/sigfit/wgs-matrix
# Rscript code/sigfit/doSigSearch.R Data/sigfit/sharedUnique/wgs-matrix-shared.rds Data/sigfit/sharedUnique/sigSearchOutput-shared.rds
# Rscript code/sigfit/doSigSearch.R Data/sigfit/sharedUnique/wgs-matrix-unique.rds Data/sigfit/sharedUnique/sigSearchOutput-unique.rds
# Rscript code/sigfit/doGofPlot.R Data/sigfit/sharedUnique/sigSearchOutput-shared.rds Data/sigfit/sharedUnique/wgs-matrix-shared.rds plots/sigfit/gof-shared.png
# Rscript code/sigfit/doGofPlot.R Data/sigfit/sharedUnique/sigSearchOutput-unique.rds Data/sigfit/sharedUnique/wgs-matrix-unique.rds plots/sigfit/gof-unique.png
# Rscript code/sigfit/fitSigN.R Data/sigfit/sharedUnique/wgs-matrix-shared.rds 5 Data/sigfit/sharedUnique/extractedSigs-wgs-shared.rds
# Rscript code/sigfit/fitSigN.R Data/sigfit/sharedUnique/wgs-matrix-unique.rds 5 Data/sigfit/sharedUnique/extractedSigs-wgs-unique.rds
# Rscript code/sigfit/getSigsExpsFromMcmcFit.R Data/sigfit/sharedUnique/extractedSigs-wgs-shared.rds Data/sigfit/sharedUnique/extractedSigs-wgs-shared-expo+sigs.rds
# Rscript code/sigfit/getSigsExpsFromMcmcFit.R Data/sigfit/sharedUnique/extractedSigs-wgs-unique.rds Data/sigfit/sharedUnique/extractedSigs-wgs-unique-expo+sigs.rds
