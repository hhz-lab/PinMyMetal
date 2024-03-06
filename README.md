# PinMyMetal: Accurately model metal binding sites through a hybrid learning system


If using this work please cite:
>PinMyMetal: A hybrid learning system to accurately model metal binding sites in macromolecules

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
cd PinMyMetal/zinc_prediction
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
cp /your_path/PinMyMetal/zinc_prediction/structures.py /path/to/PinMyMetal/environment/atomium/structures.py  # Replace the structures.py file
```
**4. Run the shell script to complete the prediction**
```
cd zinc_prediction/script/excute
python3 excute.py -p PDB_id 
python3 excute.py -u uniprot_id
python3 excute.py -f PDB_file
```

### For example, using the PDB structure 3mnd as input, the output results are saved in the 'output_data' directory.

`python3 excute.py -p 3mnd`

No non-standard hardware is required.
Installation and prediction process can typically be completed within 1 hour on a "normal" desktop computer.

# Data
The PDB codes and relevant characteristic data used for training and testing are available in `data`.

# Note
If you want to predict CIF (Crystallographic Information File) formatted files, you need to convert them to PDB (Protein Data Bank) format. Follow the steps below：
```
cd PinMyMetal/zinc_prediction
zcat ciftr-v2.053-prod-bin-linux.tar.gz | tar -xf -
cd ciftr-v2.053-prod-bin-linux
RCSBROOT=your_path/PinMyMetal/zinc_prediction/ciftr-v2.053-prod-bin-linux/
export RCSBROOT
PATH="$RCSBROOT/bin:"$PATH
export PATH
CIFTr [-i filename] -extension pdb
```
For more detailed information on using CIFTr to translate files from mmCIF format to PDB format, please refer to [RCSB PDB CIFTr](https://sw-tools.rcsb.org/apps/CIFTr/)
# License
All code is licensed under MIT license.

