#!/bin/bash

filename="$1"
pdb="$2"

# Get the directory where the script is located
SCRIPT_DIR=$(dirname "$(readlink -f "$0")")

input_path="${SCRIPT_DIR}/execute/input_data"
output_path="${SCRIPT_DIR}/execute/output_data"

pdbfiledir="${SCRIPT_DIR}/${pdb}_nb_result"
mkdir -p ${pdbfiledir}
touch ${pdbfiledir}/XXXX.cluster
mv ${input_path}/${filename}.ent  ${pdbfiledir}/pdb${pdb}.ent

# create database $DB_NAME 
DB_NAME=metal${pdb}
DB_USER="myuser"
DB_PASSWORD="123456"
DB_HOST="postgresql_container"
DB_PORT=5432
export PGPASSWORD=$DB_PASSWORD
psql -h $DB_HOST -U $DB_USER -d postgres -p $DB_PORT -c "CREATE DATABASE $DB_NAME;"

# Connecting to database $DB_NAME
# Define a function to execute an SQL script
run_sql_script() {
  local sql_file=$1  # The SQL file to execute is passed as an argument
  psql -h $DB_HOST -U $DB_USER -d ${DB_NAME} -p $DB_PORT < "$sql_file"  # Execute the SQL file using psql
}

# Call the function to execute multiple SQL files
run_sql_script "${SCRIPT_DIR}/utils/util.sql"
run_sql_script "${SCRIPT_DIR}/utils/c_value_exp_chedh.sql"

#SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
#echo "$SCRIPT_DIR"
# Set the NEIGHBORHOOD environment variable
export NEIGHBORHOOD="${SCRIPT_DIR}/neighborhood_0.8.0_src/"

# Change directory to the src folder within the neighborhood_0.8.0_src directory
cd ${SCRIPT_DIR}/neighborhood_0.8.0_src/src/

# Run the neighborhood program
./neighborhood ${pdbfiledir} pdb${pdb}.ent XXXX.cluster NULL NULL NULL
python3 ${SCRIPT_DIR}/import_neighborhood.py -i ${pdb}


cd ${pdbfiledir}

# prediction of metal ligands
run_sql_script "${SCRIPT_DIR}/script1_predict.sql"

# get metal coord
python3 ${SCRIPT_DIR}/metal_coord_pre/CH_metal_coord23.py -i ${pdb}
python3 ${SCRIPT_DIR}/metal_coord_pre/CH_metal_coord4.py -i ${pdb}
python3 ${SCRIPT_DIR}/metal_coord_pre/EDH_metal_coord2345.py -i ${pdb}
python3 ${SCRIPT_DIR}/metal_coord_pre/EDH_metal_coord6.py -i ${pdb}

run_sql_script "${SCRIPT_DIR}/script2_angle_dist.sql"


# hydrophobicity
python3 ${SCRIPT_DIR}/hybrid_algorithm/feature_pre_chedh.py -i ${pdb}
run_sql_script "${SCRIPT_DIR}/hybrid_algorithm/script3_hydrophobic.sql"

# hydro_proba
python3 ${SCRIPT_DIR}/hybrid_algorithm/hydro_proba.py -i ${pdb}
run_sql_script "${SCRIPT_DIR}/hybrid_algorithm/script4_hydro_proba.sql"

# geometric feature
python3 ${SCRIPT_DIR}/hybrid_algorithm/calculate_angle_dist_site234.py -i ${pdb}
run_sql_script "${SCRIPT_DIR}/hybrid_algorithm/script5_geome_feature.sql"

# machine learning
python3 ${SCRIPT_DIR}/hybrid_algorithm/predict_ch_metal.py -i ${pdb}
python3 ${SCRIPT_DIR}/hybrid_algorithm/predict_edh_metal.py -i ${pdb}
run_sql_script "${SCRIPT_DIR}/hybrid_algorithm/script6_ml_result.sql"

# Determine which chain the predicted metal belongs to
python3 ${SCRIPT_DIR}/hybrid_algorithm/get_premetal_chainid.py -i ${pdb}
# metals cluster
python3 ${SCRIPT_DIR}/hybrid_algorithm/match_exp_to_pre.py -i ${pdb}
python3 ${SCRIPT_DIR}/hybrid_algorithm/exp_presites_clusterer.py -i ${pdb}
run_sql_script "${SCRIPT_DIR}/hybrid_algorithm/script7_metal_chain_cluster.sql"

# get premetal pdbfile
python3 ${SCRIPT_DIR}/hybrid_algorithm/add_metal_site_to_pdb.py -i ${pdb}
run_sql_script "${SCRIPT_DIR}/hybrid_algorithm/script8_premetal_pdbfile.sql"

# execute neighborhood
cd ${SCRIPT_DIR}/neighborhood_0.8.0_src/src/
touch ${pdbfiledir}/premetal/XXXX.cluster
# Run the neighborhood program
./neighborhood ${pdbfiledir}/premetal pdb${pdb}.ent XXXX.cluster NULL NULL NULL

cd ${pdbfiledir}

# chemical feature
python3 ${SCRIPT_DIR}/hybrid_algorithm/classmodel_features.py -i ${pdb}
# get ionbinding features
python3 ${SCRIPT_DIR}/hybrid_algorithm/get_ionbinding_data.py -i ${pdb}
run_sql_script "${SCRIPT_DIR}/hybrid_algorithm/script9_get_ionbinding_data.sql"

# chedh clasmodel
python3 ${SCRIPT_DIR}/hybrid_algorithm/predict_ched_metalclass.py -i ${pdb}
run_sql_script "${SCRIPT_DIR}/hybrid_algorithm/script10_clasmodel_result.sql"

# show metal sites
python3 ${SCRIPT_DIR}/add_metal_site_to_pdb2.py -i ${pdb}

# The script did not filter out the metal ions already present in the structure when removing duplicate predicted sites.
#python3 ${SCRIPT_DIR}/add_metal_site_to_pdb_contain_expsites.py -i ${pdb}


mv ${output_path}/${pdb}_metal.pdb ${output_path}/${filename}_metal.pdb
mv ${output_path}/${pdb}_output.csv ${output_path}/${filename}_output.csv

psql -h $DB_HOST -U $DB_USER -d postgres -p $DB_PORT -c "DROP DATABASE IF EXISTS $DB_NAME;"

echo 'end'
