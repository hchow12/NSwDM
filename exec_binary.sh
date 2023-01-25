#! /bin/bash

cp binary_conf.c binary/
cd binary/
echo "Compiling with gcc..."
gcc -o run binary_conf.c evolve.c -lgsl -lgslcblas -lm 

echo "Evolving the binary..."
./run

# Move the binary configurations file to folder "head/" in ../
echo "Saving the binary configurations file..."
mv HEAD* ../head/

# Copy data in binary/ to the folder "data/" in ../
cp -r data ../

# Cleaning the folder "binary/data/"
rm -vf data/bh/no_dark/*
rm -vf data/bh/dark/*
rm -vf data/dark/*
rm -vf data/dm/*
