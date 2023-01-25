#! /bin/bash

typ=$1
file=$2

# Command line argument checking
if [ $# -eq 0 ]
then
	echo "Error: 0 argument given."
	echo "Usage: ./find_gw_ft.sh <dm/dark> <datafile name>"
	exit 1
elif [ $# -eq 1 ]
then
	echo "Error: only 1 argument given."
	echo "Usage: ./find_gw_ft.sh <dm/dark> <datafile name>"
	exit 1
fi

target=GW_$file

echo finding if $target exist...
cd data/gw/$typ/
if [ -e $target ]
then
	echo GW file found
else
	echo GW file not found, now producing the GW and FT file...
	cd ../../../
	./run_gw_ft.sh $typ $file
fi
