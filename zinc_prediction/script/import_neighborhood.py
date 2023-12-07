#!/usr/bin/python
import getopt
import os, sys


options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['intput='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg

pdbid=inputid

script_dir="/zinc_prediction/script"


data_dir = script_dir + '/' +str(pdbid)+ "_nb_result"

def psql_virusMED(sqlfile):
    dbname = "zinc"+str(pdbid)
    os.system("psql " + dbname + " < "+sqlfile)

psql_virusMED(os.path.join(script_dir, "createNeighborhoodTables.sql"))

sqlfile = os.path.join(data_dir, "copyNeighborhoodData.sql")

psql_virusMED(sqlfile)
