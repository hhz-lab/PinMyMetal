import pandas as pd
import numpy as np
from scipy.spatial.distance import cdist
import psycopg2 as pg
from itertools import combinations
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

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'metalsites_cluster.csv')
conn = pg.connect("dbname="+DBname+" password='' port='5432' host='/var/run/postgresql'")

# Execute SQL queries
sql1 = """
SELECT DISTINCT 'edh' AS sitetype, pdbid, id, metal_x, metal_y, metal_z, residueid_a, residueid_b, residueid_c, residueid_d, residueid_e, residueid_f,
    proba_ismetal
    FROM metal_predict_sites_2 where mvc_result=1
"""
sql2 = """
SELECT DISTINCT 'ch' AS sitetype, pdbid, id, metal_x, metal_y, metal_z, residueid_a, residueid_b, residueid_c, residueid_d,
    proba_ismetal
    FROM zncu_predict_sites_2 where mvc_result=1
"""

sql3 = """
select distinct 'exp' as sitetype,pdbid,id, x_ion as metal_x, y_ion as metal_y, z_ion as metal_z,
    string_to_array(conc_comma(residueid_lig::text),', ') as residueids_lig,atomname_ion from
    (select distinct pdbid,id, atomname_ion, x_ion,y_ion,z_ion, residueid_lig from transition_metal_coord) a
    group by pdbid,id,atomname_ion, x_ion,y_ion,z_ion
"""
# where resname_lig in ('CYS','HIS','GLU','ASP')

# Read SQL query results into pandas DataFrames
df1 = pd.read_sql(sql1, conn)
df2 = pd.read_sql(sql2, conn)
df3 = pd.read_sql(sql3, conn)

# Combine the data from both queries
df = pd.concat([df1, df2, df3],ignore_index=True)


def generate_cluster_id_within_pdbid(data):
    pdbid_list = data['pdbid'].unique()
    cluster_id = np.zeros(len(data), dtype=int)
    current_cluster = 1

    for pdbid in pdbid_list:
        pdbid_data = data[data['pdbid'] == pdbid]
        coordinates = pdbid_data[['metal_x', 'metal_y', 'metal_z']].values

        for i, j in combinations(range(len(pdbid_data)), 2):
            if np.linalg.norm(coordinates[i] - coordinates[j]) <= 2.5:
                if cluster_id[pdbid_data.index[i]] == 0 and cluster_id[pdbid_data.index[j]] == 0:
                    cluster_id[pdbid_data.index[i]] = current_cluster
                    cluster_id[pdbid_data.index[j]] = current_cluster
                    current_cluster += 1
                elif cluster_id[pdbid_data.index[i]] == 0:
                    cluster_id[pdbid_data.index[i]] = cluster_id[pdbid_data.index[j]]
                elif cluster_id[pdbid_data.index[j]] == 0:
                    cluster_id[pdbid_data.index[j]] = cluster_id[pdbid_data.index[i]]
                else:
                    cluster_id[cluster_id == cluster_id[pdbid_data.index[j]]] = cluster_id[pdbid_data.index[i]]

        # Assign a new cluster_id for points not assigned to any cluster within the pdbid
        for k in range(len(pdbid_data)):
            if cluster_id[pdbid_data.index[k]] == 0:
                cluster_id[pdbid_data.index[k]] = current_cluster
                current_cluster += 1

    return cluster_id


df['cluster_id'] = generate_cluster_id_within_pdbid(df)

# Filter out -9999 ligands and generate a list of ligands for each row
def filter_ligands(row):
    if row['sitetype'] in ['edh', 'ch']:
        ligands = [row['residueid_a'], row['residueid_b'], row['residueid_c'], row['residueid_d']]
        if 'residueid_e' in row:
            ligands.append(row['residueid_e'])
        if 'residueid_f' in row:
            ligands.append(row['residueid_f'])
    elif row['sitetype'] == 'exp':
        ligands = row['residueids_lig']
    return [int(ligand) for ligand in ligands if pd.notna(ligand) and ligand != -9999]

df['ligands'] = df.apply(filter_ligands, axis=1)

df['is_representative'] = False


for cluster_id, group in df.groupby('cluster_id'):
    exp_sites = group[group['sitetype'] == 'exp']
    if not exp_sites.empty:
        non_exp_sites = group[group['sitetype'] != 'exp']
        if not non_exp_sites.empty:
            coords_exp = exp_sites[['metal_x', 'metal_y', 'metal_z']].values
            coords_non_exp = non_exp_sites[['metal_x', 'metal_y', 'metal_z']].values
            dist_matrix = cdist(coords_non_exp, coords_exp, metric='euclidean')
            min_dist_idx = dist_matrix.min(axis=1).argmin()
            rep_idx = non_exp_sites.iloc[min_dist_idx].name
            min_dist = dist_matrix.min(axis=1)[min_dist_idx]
            df.loc[rep_idx, 'is_representative'] = True
            df.loc[rep_idx, 'min_distance_to_exp'] = min_dist
    else:
        df['atomname_ion'] = 'NaN'
        df['exp_site'] = False
        df['min_distance_to_exp'] = -1
        max_proba_idx = group['proba_ismetal'].idxmax()
        df.loc[max_proba_idx, 'is_representative'] = True

df = df.dropna(subset=['sitetype'])
df['exp_site'] = df.groupby('cluster_id')['sitetype'].transform(lambda x: 'exp' in x.values)
#df['exp_site'] = df.groupby('cluster_id')['sitetype'].apply(lambda x: 'exp' in x.values).reset_index(drop=True)

output_df = df[['pdbid', 'cluster_id', 'id', 'sitetype', 'ligands', 'is_representative', 'atomname_ion', 'exp_site','min_distance_to_exp']]

output_df.to_csv(file_path, sep='\t', index=False)

