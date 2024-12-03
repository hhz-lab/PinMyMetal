#!/bin/sh
#
#  File:  binary.sh
#
#  Purpose: Prepares binary data from the compressed ASCII data.
#
#  Assumptions: Script is started from within the etc directory.

# Change mode of data directory recursively to allow full access.
chmod -R 777 ../data
 
# Go to ASCII data directory
cd ../data/ascii

# Unbundle ASCII data files
echo Extracting ASCII data files...
for file in `ls *.gz`
do
    gunzip $file
done
echo Done.

echo Creating binary data files from ASCII data files. Please wait...

# Go to root module directory
cd ../..

# Get component.cif from data/ascii directory

cp data/ascii/component.cif .

# Create component.odb

./bin/connect -all

# Move the result to the binary directory
mv component_new.odb data/binary/component.odb

# Cleanup all files
rm -f component_all.cif
rm -f component.cif
rm -f connect_admin.err
rm -f CifParser.log


# Start main binary production
./bin/binary_main -dic -path .


# Cleanup all files
rm -f CifParser.log
rm -f DDLParser.log
rm -f DICParser.log

# Change mode of data directory recursively to allow full access.
# This is so that the newly created files are included.
chmod -R 777 data
 
echo Done.
