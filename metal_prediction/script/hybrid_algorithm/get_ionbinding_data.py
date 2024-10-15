import os, sys
import psycopg2 as pg
import getopt

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")
cur = conn.cursor()

data_dir = os.path.join(script_directory, '../', f"{pdbid}_nb_result")
file_a = os.path.join(data_dir, 'presite_feature_ion.txt')
file_b = os.path.join(data_dir, 'presite_feature_info.txt')
output_a = open(file_a, "w")
output_b = open(file_b, "w")

presites_neighbor_dir = os.path.join(data_dir,'premetal')

query = """
SELECT DISTINCT id, sitetype, pdbid, chainid_ion, resseq_ion, residueid_ion from presites_info_neighbor order by pdbid 
"""
cur.execute(query)
results = cur.fetchall()


for result in results:
    id, sitetype, pdbid, chainid_ion, resseq_ion, residueid_ion = result

    try:
        # Read ION_BINDINGSITE_LIGATOM.data
        ligatom_file_path = os.path.join(presites_neighbor_dir, "RESIDUE.data")
        with open(ligatom_file_path, "r") as ligatom_file:
            for line in ligatom_file:
                parts = line.strip().split('|')
                if parts[4] == chainid_ion.strip() and int(parts[5]) == resseq_ion and parts[3] =="_ZN":
                    pdbfileid = parts[0]
                    residueid = parts[1]

                    try:
                        # Read ION_BINDINGSITE.data
                        bindingsite_file_path = os.path.join(presites_neighbor_dir, "ION_BINDINGSITE.data")
                        with open(bindingsite_file_path, "r") as bindingsite_file:
                            for bs_line in bindingsite_file:
                                bs_parts = bs_line.strip().split('|')
                                if bs_parts[3] == residueid:
                                    output_a.write(bs_line)
                                    output_b.write(f"{id}|{sitetype}|{pdbid}|{chainid_ion}|{resseq_ion}|{pdbfileid}|{residueid}\n")
                    except FileNotFoundError:
                        print(f"File not found: {bindingsite_file_path}, skipping to next result.")
                        continue

    except FileNotFoundError:
        print(f"File not found: {ligatom_file_path}, skipping to next result.")
        continue


output_a.close()
output_b.close()

