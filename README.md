# PinMyMetal: A hybrid learning system to accurately model transition metal binding sites in macromolecules


If using this work please cite:
>PinMyMetal: A hybrid learning system to accurately model transition metal binding sites in macromolecules

# How to use PinMyMetal 
For easy and quick use without installation, use [the PinMyMetal web server](https://PMM.biocloud.top)

You can also run it locally. Command line usage is described below.

# How to install and run locally


## Create a PostgreSQL database

### If you already have a PostgreSQL database with a version greater than or equal to 10, you can skip this step.

**Linux Platform installation**
```
sudo apt-get install postgresql-15
```
For other installation methods, please visit [PostgreSQL official website](https://www.postgresql.org/)


## To configure the environment for running PinMyMetal, you can create a Conda environment using the provided environment.yml file. Follow these steps:

**0. Ensure Conda is installed on your system.**

**1. Run the following commands in your terminal:** 
```
cd PinMyMetal/metal_prediction
conda env create -f environment.yml
```
**2. Activate the new environment using the following command:**
```
conda activate PinMyMetal
```
**3. optimizing Atomium**

#### Atomium is a molecular modeller and file parser, capable of reading from and writing to .pdb, .cif and .mmtf files.
#### Please perform the following actions to update the structures.py file in the PinMyMetal environment, as we have made modifications to the source code of Atomium to enhance its execution speed.

```
conda info --envs  # Find the location of the PinMyMetal environment
cd /path/to/PinMyMetal/environment  # Navigate to the PinMyMetal environment directory
find . -type d -name "atomium"  # Find the installation directory of the atomium package
cp /your_path/PinMyMetal/metal_prediction/structures.py /path/to/atomium/structures.py  # Replace the structures.py file
```
**4. Run the shell script to complete the prediction**
```
cd script/excute
python3 excute.py -p PDB_id 
python3 excute.py -u uniprot_id
python3 excute.py -f PDB_file
```

### For example, using the PDB structure 3mnd as input, the output results are saved in the 'output_data' directory.

`python3 excute.py -p 2zp9`

### When using the -f option to run your own uploaded file, please ensure the file is uploaded to the input_data directory. For example, to use an existing file in that directory, run:
`python3 excute.py -f 7pw5.pdb`

### One PDB file containing the predicted metal information is named X_metal.pdb, and another file storing detailed information about the binding sites is named X_output.file. For example, you may find 7pw5_metal.pdb and 7pw5_output.csv in the output_data directory. Generally, you can use only the X_metal.pdb file to analyze the prediction results. If any metal binding sites are predicted, they will appear at the end of the file, with the insertion code in column 27 marked with "@" to indicate the metals predicted by PMM. Corresponding ligands can be found in the LINK section.

No non-standard hardware is required.
Installation and prediction process can typically be completed within 1 hour on a "normal" desktop computer.

# Data
The PDB codes and relevant characteristic data used for training and testing are available in `data_model`.

# Note
#### Please replace the original structures.py file in the Atomium software with the provided structures.py file to ensure the use of customized structure processing features and the proper functioning of the program.

If you want to predict CIF (Crystallographic Information File) formatted files, you need to convert them to PDB (Protein Data Bank) format. Follow the steps below:
```
cd metal_prediction
zcat ciftr-v2.053-prod-bin-linux.tar.gz | tar -xf -
RCSBROOT=metal_prediction/ciftr-v2.053-prod-bin-linux/
export RCSBROOT
PATH="$RCSBROOT/bin:"$PATH
export PATH
CIFTr [-i filename] -extension pdb
```
# License
All code is licensed under MIT license.

