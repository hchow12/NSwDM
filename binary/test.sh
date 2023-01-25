#! /bin/bash

echo "Compiling test.c file with gcc..."
gcc -o test test.c -lgsl -lgslcblas -lm

echo "Executing the test executable..."
./test

