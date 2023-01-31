#! /bin/bash

# Check if "PSD/" folder exist, if not, then creat the folder
if [ ! -d "signal_psd" ]
then
	mkdir signal_psd
fi

typ=$1
file=$2

# Command line argument checking
if [ $# -eq 0 ]
then
    echo "Error: 0 argument given."
    echo "Usage: ./single_snr.sh <dm/dark> <fourier datafile name>"
    exit 1
elif [ $# -eq 1 ]
then
    echo "Error: only 1 argument given."
    echo "Usage: ./single_snr.sh <dm/dark> <fourier datafile name>"
    exit 1
fi

echo finding if $file exist...
cd ../data/fourier/$typ/
if [ -e $file ]
then
	echo F_GW file found, calculating SNR..
	cp $file ../../../SNR/
	cd ../../../SNR/
	./snr $file
fi

# Cleaning the working director
rm F_GW*


