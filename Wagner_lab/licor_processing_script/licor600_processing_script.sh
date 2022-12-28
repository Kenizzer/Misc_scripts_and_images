#!/bin/bash

read -e -p "Enter filename, use tab for completion: " file


# Reduced output
# only obs, gsw, e_apparent, Phips2, tleaf, ETR, fm'

cat $file | cut -f 1,2,10,13,28,35,31,27 -d "," | tail -n +2 | grep -v -w HHMMSS > ${file%.csv}.thin_output.csv


# Fat output
# Basically everything but the blank columns and junk

cat $file | cut -f 1,2,10-21,26-40 -d "," | tail -n +2 | grep -v -w HHMMSS >  ${file%.csv}.fat_output.csv

# Clean up all the individual files
# I don't currently have a use for these.
mkdir -p indivdual_files
find PSF* -type f -exec mv -t indivdual_files {} +

