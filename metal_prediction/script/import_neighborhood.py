#!/usr/bin/python
import getopt
import os, sys

script_dir = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['intput='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid

data_dir = os.path.join(script_dir, f"{pdbid}_nb_result")

def psql_db(sqlfile):
    dbname = "metal"+str(pdbid)
    os.system("psql " + dbname + " < "+sqlfile)

psql_db(os.path.join(script_dir, "createNeighborhoodPartialTables.sql"))

sqlfile = os.path.join(data_dir, "copyNeighborhoodData.sql")

psql_db(sqlfile)
