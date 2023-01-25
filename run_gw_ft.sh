#! /bin/bash

cd gw_ft/
echo "Compiling cal_gw.c with gcc..."
gcc -o exec_gw cal_gw.c -lgsl -lgslcblas -lm

typ=$1
file=$2
cp ../data/$typ/$file ./
if [ $? -ne 0 ]
then
	echo "Probably something wrong with datafile name or directory <dm/dark>, exiting..."
	exit 1
fi

# $file will be processed to get GW...

echo "Running gw..."
./exec_gw $file

# "Renaming Output.txt to GW_file..."
mv Output.txt GW_$file

echo "Compiling cal_ft.c with gcc..."
gcc -o exec_ft cal_ft.c -lgsl -lgslcblas -lm

echo "Running ft..."
./exec_ft GW_$file

# "Renaming FT_Output.txt to FT_file..."
mv FT_Output.txt F_GW_$file

# "Moving the GW file to database..."
mv GW_$file ../data/gw/$typ/

# "Moving the FT file to database..."
mv F_GW_$file ../data/fourier/$typ/

echo "Cleaning the working directory..."
rm -vf $file

