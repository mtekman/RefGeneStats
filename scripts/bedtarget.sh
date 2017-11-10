#!/bin/bash

if [ $# -lt 1 ]; then
	echo "This is a wrapper script for the bedtarget2 binary,\
 to remove illegaly specified loci (i.e. start==stop position) that was too much trouble to\
 check for in the binary (believe it or not)." >&2
	exit
fi

scripts/bedtarget2 $* | awk '{if ($2!=$3) print $0}' 
