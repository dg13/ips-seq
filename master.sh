#!/bin/sh

## 1. Copy over raw call files
SRCDIR=/lustre/scratch116/vr/user/pd3/hipsci/exome-point-mutations/unified-set/final-set-2018-02-05.all
rsync -av --progress $SRCDIR/exomes.ft.txt.gz Data/newCalls2/raw/
rsync -av --progress $SRCDIR/high-vs-low-exomes.ft.txt.gz Data/newCalls2/raw/
rsync -av --progress $SRCDIR/wgs.396.ft.txt.gz Data/newCalls2/raw/
