#!/bin/bash
rm /var/www/django/zinc/log/*.txt
filename="$1"
pdb="$2"

input_path=/var/www/django/zinc/static/split_chain 

#pdb=$(echo ${filename: 3: 4})

log_dir=/var/www/django/zinc/log
data_dir=/zinc_prediction/script
pdbfiledir=${data_dir}/${pdb}_nb_result
mkdir -p ${pdbfiledir}
touch ${pdbfiledir}/XXXX.cluster
mv $input_path/$filename $pdbfiledir/pdb${pdb}.ent
mv /var/www/django/zinc/static/input_data/$filename $pdbfiledir/pdb2${pdb}.ent
cd ${data_dir}

dbname=zinc${pdb}
createdb ${dbname}
touch $log_dir/01.txt

psql ${dbname} < ${data_dir}/util.sql 
psql ${dbname} < ${data_dir}/c_value_exp.sql # new add

export NEIGHBORHOOD=${data_dir}/neighborhood_0.8.0_src/
cd ${data_dir}/neighborhood_0.8.0_src/src/

./neighborhood ${pdbfiledir} pdb${pdb}.ent XXXX.cluster NULL NULL NULL 
python3 ../../import_neighborhood.py -i ${pdb} 
touch $log_dir/02.txt

cd ${pdbfiledir}

# prediction of zinc ligands
psql ${dbname} < ../script1_predict.sql

# get zinc coord
python3 ../zinc_coord_pre/His_zinc_coord_23.py -i ${pdb} 
python3 ../zinc_coord_pre/cys_his_zinc_coord23.py -i ${pdb}
python3 ../zinc_coord_pre/cys_zinc_coord23.py  -i ${pdb}
python3 ../zinc_coord_pre/zinc_coord4.py  -i ${pdb}
psql ${dbname} < ../script2_angle_dist.sql


# hydrophobicity
python3 ../zinc_hydro_pre2.py -i ${pdb} # different
psql ${dbname} < ../script3_hydrophobic.sql
python3 ../hydro_proba.py -i ${pdb} # new add
touch $log_dir/03.txt
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
touch $log_dir/04.txt
output_file=/var/www/django/zinc/static/output_data
OLD_IFS="$IFS"
IFS="."
array=(${filename})
IFS="$OLD_IFS"

mv ${output_file}/${pdb}_zinc.pdb ${output_file}/${array[0]}_zinc.pdb
mv ${output_file}/${pdb}_output.csv ${output_file}/${array[0]}_output.csv

dropdb zinc${pdb}
rm -r ${pdbfiledir}
