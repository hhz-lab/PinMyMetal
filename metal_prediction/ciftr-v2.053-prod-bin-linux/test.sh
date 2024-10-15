#!/bin/sh

ls -1 testdata/*.Z > tmplist
$RCSBROOT/bin/CIFTr -f tmplist -extension pdb -uncompress uncompress -compress gzip -output_path testresult
rm tmplist
rm CifParser.log
rm *.err
