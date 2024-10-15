import psycopg2 as pg
from collections import Counter
import os, sys
import getopt

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'get_premetal_chainid.csv')
print_log = open(file_path,'w')
conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()


sql = """
    select distinct 'edh' as sitetype, id,chainid_a,chainid_b,chainid_c,chainid_d,chainid_e,chainid_f from metal_predict_sites_2 where mvc_result=1 union
    select distinct 'ch' as sitetype, id,chainid_a,chainid_b,chainid_c,chainid_d,'NaN' as chainid_e, 'NaN' as chainid_f from zncu_predict_sites_2 where mvc_result=1
"""

cur.execute(sql)
data = cur.fetchall()

# Function to select the chain based on given criteria
def select_chain(*chains):
    # Count occurrences of each chain
    counter = Counter(chains)
    # Remove 'NaN' from the count
    if 'NaN' in counter:
        del counter['NaN']
    # If no valid chains are left, return None
    if not counter:
        return None
    # Get the most common chains
    most_common = counter.most_common()
    # If the top two counts are equal, select the smallest value
    if len(most_common) > 1 and most_common[0][1] == most_common[1][1]:
        return sorted([chain for chain, count in most_common if count == most_common[0][1]])[0]
    # Otherwise, return the most common chain
    return most_common[0][0]

# Process each row of data
for row in data:
    sitetype= row[0]
    pid = row[1]
    chains = row[2:8]  # Get chainid_a to chainid_f
    selected_chain = select_chain(*chains)
    print("%s\t%s\t%s"%(sitetype,pid,selected_chain), file=print_log)
