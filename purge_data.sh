#! /bin/bash

echo "Deleting all data file including: fourier, gw_ft, head..."

rm -f head/*
rm -f data/bh/dark/*
rm -f data/bh/no_dark/*
rm -f data/dark/*
rm -f data/dm/*
rm -f data/fourier/dark/*
rm -f data/fourier/dm/*
rm -f data/gw/dark/*
rm -f data/gw/dm/*
