#!/bin/sh

./code/preprocess-mutation-calls/doAllMakeNewCallsDataStruct.sh Data/mut-files Data/mut-files/rdas wes farm clobber
./code/preprocess-mutation-calls/doAllMakeNewCallsDataStruct.sh Data/mut-files Data/mut-files/rdas hwes farm clobber
./code/preprocess-mutation-calls/doAllMakeNewCallsDataStruct.sh Data/mut-files Data/mut-files/rdas wgs farm clobber
