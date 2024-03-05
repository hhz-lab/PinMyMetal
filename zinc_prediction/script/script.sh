#!/bin/bash

filename="$1"
pdb="$2"

input_path=./excute/input_data
output_path=./excute/output_data

cd ../
data_dir=./

pdbfiledir=${pdb}_nb_result
mkdir -p ${pdbfiledir}
touch ${pdbfiledir}/XXXX.cluster
mv $input_path/${filename}.ent  $pdbfiledir/pdb${pdb}.ent

dbname=zinc${pdb}
createdb ${dbname}


psql ${dbname} < ${data_dir}util.sql
psql ${dbname} < ${data_dir}c_value_exp.sql 

SCRIPT_DIR=$(dirname "$(readlink -f "$0")")
echo "$SCRIPT_DIR"
export NEIGHBORHOOD="$SCRIPT_DIR/script/neighborhood_0.8.0_src/"
cd neighborhood_0.8.0_src/src/

./neighborhood ./../../${pdbfiledir} pdb${pdb}.ent XXXX.cluster NULL NULL NULL
python3 ../../import_neighborhood.py -i ${pdb}


cd ./../../${pdbfiledir}

# prediction of zinc ligands
psql ${dbname} < ../script1_predict.sql

# get zinc coord
python3 ../zinc_coord_pre/His_zinc_coord_23.py -i ${pdb}
python3 ../zinc_coord_pre/cys_his_zinc_coord23.py -i ${pdb}
python3 ../zinc_coord_pre/cys_zinc_coord23.py  -i ${pdb}
python3 ../zinc_coord_pre/zinc_coord4.py  -i ${pdb}
psql ${dbname} < ../script2_angle_dist.sql


# hydrophobicity
python3 ../zinc_hydro_pre.py -i ${pdb}
psql ${dbname} < ../script3_hydrophobic.sql
python3 ../hydro_proba.py -i ${pdb} 

#machine learning
python3 ../ml_predict.py -i ${pdb}

psql ${dbname} < ../script4_expzinc.sql


# remove redundant site
python3 ../zinc_coord_pre/exp_pre_dist.py  -i ${pdb}
psql ${dbname} < ../zinc_coord_pre/de_redundant_site_v1.sql
python3 ../zinc_coord_pre/pre_pre_dist.py  -i ${pdb}
psql ${dbname} < ../zinc_coord_pre/de_redundant_site_v2.sql
python3 ../zinc_coord_pre/contact_ab.py  -i ${pdb}
psql ${dbname} < ../zinc_coord_pre/de_redundant_site_v3.sql


python3 ../predict_zinc_script.py -i ${pdb}

mv ../${output_path}/${pdb}_zinc.pdb ../${output_path}/${filename}_zinc.pdb
mv ../${output_path}/${pdb}_output.csv ../${output_path}/${filename}_output.csv

echo 'end'

