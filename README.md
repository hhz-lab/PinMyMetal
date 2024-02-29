# PinMyMetal: Accurately model metal binding sites through a hybrid learning system


If using this work please cite:
>PinMyMetal: A hybrid learning system to accurately model metal binding sites in macromolecules

# How to use PinMyMetal 
For easy and quick use without installation, use [the PinMyMetal web server](https://PMM.biocloud.top)

You can also run it locally. Command line usage is described below.

# How to install and run locally


0. Create a PostgreSQL database

**Linux Platform installation**
```
sudo apt install postgresql-client-common
sudo apt install postgresql
cd /var/lib/postgresql/
sudo apt install postgresql
```
For other installation methods, please visit [PostgreSQL official website](https://www.postgresql.org/)


1. Preparing for predictions involves two steps to ensure a smooth and accurate process. Here are the steps you should take:

**First, move the zinc_prediction folder to your root directory, and then modify the directory permissions to make it executable and writable.**
```
sudo mv zinc_prediction /zinc_prediction
chmod -R 777 /zinc_prediction
```

Then optimizing Atomium
atomium is a molecular modeller and file parser, capable of reading from and writing to .pdb, .cif and .mmtf files.
As we have made modifications to the source code of Atomium to enhance its execution speed, please perform the following actions:

**Install the Atomium package using pip3.**

`pip3 install atomium`

**Use the command to find the installation location of the Atomium package. Look for the "Location" field in the output, which indicates the installation path of the Atomium package.**

`pip3 show atomium`

**This command copies the modified structures.py file to the Atomium package directory, replacing the existing file.**

`cp /zinc_prediction/zinc_prediction/structures.py /path/to/atomium/`

2. Run the shell script to complete the prediction
```
cd zinc_prediction/script/excute
python3 excute.py -p PDB_id 
python3 excute.py -u uniprot_id
python3 excute.py -f PDB_file
```

**For example, using the PDB structure 3mnd as input, the output results are saved in the 'output_data' directory.**

`python3 excute.py -p 3mnd`

No non-standard hardware is required.
Installation and prediction process can typically be completed within 1 hour on a "normal" desktop computer.

# Data
The PDB codes and relevant characteristic data used for training and testing are available in `data`.

# Note
If you want to predict CIF (Crystallographic Information File) formatted files, you need to convert them to PDB (Protein Data Bank) format. Follow the steps below：
```
cd /zinc_predict
zcat ciftr-v2.053-prod-bin-linux.tar.gz | tar -xf -
RCSBROOT=/zinc_predict/ciftr-v2.053-prod-bin-linux/
export RCSBROOT
PATH="$RCSBROOT/bin:"$PATH
export PATH
CIFTr [-i filename] -extension pdb
```
# License
All code is licensed under MIT license.

