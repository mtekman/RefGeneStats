#!/bin/bash

usage() {
	echo "Create a single map file containing the start and stop positions of exons genes"
	echo ""
	echo "`basename $0` (hgVer) <filename> [bedtarget args]"
	exit
}

folder="$1"  # dbver
name="$2"    # absolute filename


[ "$1" = "" ] && usage

extrargs="${@:3}"

mkdir -p `dirname $name`

rm $name 2> /dev/null

scripts/bedtarget.sh chr1-22,chrX --exons --database $folder $extrargs > $name
#bedtarget.sh chr16 --exons --database $folder $extrargs > $name ##  testing only


dos2unix $name 2>/dev/null
mac2unix $name 2>/dev/null

echo "FIN" >&2
