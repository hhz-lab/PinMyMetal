import pandas as pd
import numpy as np
import os, sys
import getopt
sys.path.append(os.path.join(os.path.dirname(__file__), "../utils"))
import db_utils

# Get the directory where the script file resides
script_directory = os.path.dirname(os.path.abspath(__file__))

options, remainder = getopt.getopt(sys.argv[1:], 'i:o', ['input='])
for opt, arg in options:
    if opt in ('-i', '--input'):
        inputid = arg
pdbid=inputid
DBname="metal"+str(pdbid)

file_path = os.path.join(script_directory, '../', f"{pdbid}_nb_result", 'match_exp2pre.csv')

# Get a database connection
conn = db_utils.create_connection(DBname)

sql_exp = """
SELECT DISTINCT 'exp' AS sitetype, pdbid, id, x_ion AS metal_x, y_ion AS metal_y, z_ion AS metal_z,resname_ion FROM transition_metal_coord
"""

sql_edh = """
SELECT DISTINCT 'edh' AS sitetype, pdbid, id, metal_x, metal_y, metal_z,
    proba_ismetal 
FROM metal_predict_sites_2
WHERE mvc_result=1
"""

sql_ch = """
SELECT DISTINCT 'ch' AS sitetype, pdbid, id, metal_x, metal_y, metal_z, 
    proba_ismetal
FROM zncu_predict_sites_2
WHERE mvc_result=1
"""

exp_sites = pd.read_sql(sql_exp, conn)
edh_sites = pd.read_sql(sql_edh, conn)
ch_sites = pd.read_sql(sql_ch, conn)

pred_sites = pd.concat([edh_sites, ch_sites])

results = []
for pdbid in exp_sites['pdbid'].unique():
    exp_subset = exp_sites[exp_sites['pdbid'] == pdbid]
    pred_subset = pred_sites[pred_sites['pdbid'] == pdbid]

    for _, exp_row in exp_subset.iterrows():
        exp_coords = np.array([exp_row['metal_x'], exp_row['metal_y'], exp_row['metal_z']])

        min_dist = float('inf')
        closest_pred = None

        for _, pred_row in pred_subset.iterrows():
            pred_coords = np.array([pred_row['metal_x'], pred_row['metal_y'], pred_row['metal_z']])
            distance = np.linalg.norm(exp_coords - pred_coords)

            if distance < min_dist:
                min_dist = distance
                closest_pred = pred_row

        if closest_pred is not None:
            result = {
                'exp_sitetype': exp_row['sitetype'],
                'exp_pdbid': exp_row['pdbid'],
                'exp_id': exp_row['id'],
                'exp_metal_x': exp_row['metal_x'],
                'exp_metal_y': exp_row['metal_y'],
                'exp_metal_z': exp_row['metal_z'],
                'exp_metalname':exp_row['resname_ion'],
                'pred_sitetype': closest_pred['sitetype'],
                'pred_pdbid': closest_pred['pdbid'],
                'pred_id': closest_pred['id'],
                'pred_metal_x': closest_pred['metal_x'],
                'pred_metal_y': closest_pred['metal_y'],
                'pred_metal_z': closest_pred['metal_z'],
                'pred_proba_ismetal': closest_pred['proba_ismetal'],
                'distance': min_dist
            }
            results.append(result)

results_df = pd.DataFrame(results)
results_df.to_csv(file_path, sep='\t', index=False)



