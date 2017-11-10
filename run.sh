#!/bin/bash

mapfile=refseq.map.txt

makeg=scripts/makegenemaps.sh
statg=scripts/genestats.py

chmod +x $makeg
chmod +x $statg

if ! [ -e $mapfile ]; then
    echo "Pulling data from UCSC"
    $makeg hg38 $mapfile --introns --frames --direction --refseq-gene-names
    echo ""
fi

echo "Processing..."
$statg $mapfile && echo "" &&
    echo "Stats written to:" &&
    ls gene*
