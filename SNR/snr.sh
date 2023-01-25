#!/bin/bash

output_file="snr_graph.txt" 

# Checking if output_file exist
if [ -e $output_file ]
then
	echo "$output_file exist, removing the file....."
	rm -v $output_file 
fi

# Compiling with gcc
gcc -g -w main.c tool.c simp.c -lgsl -lgslcblas -lm -o snr

# Preparing list of fourier transformed gw file
# Edit the file_path and file_list for required fourier transformed data file
file_path=../data/fourier/dark/

file_list=$(ls $file_path | grep pA | grep MA1.40 | grep MB1.40)
if [ "$(ls $file_path | wc -l)" -eq 0 ]
then
	echo "The folder in $file_path is empty!"
	exit 1
fi

# Obtain the x_axis value
# Edit the variable x to read the point for x_axis in graphing SNR from the data file name
echo Writing new snr_graph.txt file...
for file in $file_list
do
	cp $file_path$file ./
	x=$(echo $file | cut -f 14 -d '_' | cut -c 3-)
	if [ $x == 0.80 ]; then
		break
	fi
	y=$(./snr $file)
	echo $x $y >> $output_file
	rm $file
done
echo Done SNR calculation
