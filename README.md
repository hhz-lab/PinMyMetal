# PinMyMetal: A hybrid learning system to accurately model transition metal binding sites in macromolecules


If using this work please cite:
>PinMyMetal: A hybrid learning system to accurately model transition metal binding sites in macromolecules

# How to use PinMyMetal 
For easy and quick use without installation, use [the PinMyMetal web server](https://PMM.biocloud.top)

You can also run it locally. Command line usage is described below.

# How to install and run locally
This repository provides a fully containerized environment for the PinMyMetal project. Users can either:

- **Pull pre-built Docker images from GitHub Container Registry** for quick setup, or  
- **Build the images locally** using the provided `Dockerfile` and configuration files. 

Both options ensure seamless replication of the project environment for reliable local execution.

## Installation Options: Pre-built Images or Local Build

### Step 1: Install Docker
Ensure that Docker is installed on your system. You can follow the official [Docker installation guide](https://docs.docker.com/get-docker/) for your platform.

### Step 2: Build the Docker Image for Conda
You can either use the pre-pulled image from GitHub Container Registry or build the image locally.
- **Option 1: Using Pre-built Image**
  If you prefer not to build the image yourself, you can pull the pre-built image from GitHub Container Registry:
  ```bash
  docker pull ghcr.io/hhz-lab/pinmymetal-conda-env:latest
  ```
- **Option 2: Building the Image Locally**
  Ensure you are in the `PinMyMetal/metal_prediction` directory and run:
  ```bash
  docker build -t pinmymetal-conda-env -f conda_dockerfile .
  ```
Both options will create functionally equivalent images, but the local build allows for customization of the image.

### Step 3: Modify Docker Configuration (if needed)
- **If using the pre-built image, no changes are needed to**
  `docker-compose.yml`.
- **If building the image locally, comment out the image:**
  `ghcr.io/hhz-lab/pinmymetal-conda-env:latest` line and replace it with **image**: `pinmymetal-conda-env`.
  
### Step 4: Start Containers
Navigate to the `metal_prediction` directory:
```bash
cd metal_prediction
```
Start the containers using `docker-compose`:
```bash
docker-compose up -d
```

### Step 5: Access the Conda Environment
Enter the Conda container:
```bash
docker exec -it conda_container bash
```
Activate the Conda environment:
```bash
conda activate PinMyMetal
```

### Step 6: Execute Scripts
Navigate to the `execute` directory and run the required scripts:
```bash
cd execute
python3 excute.py -p PDB_id 
python3 excute.py -u uniprot_id
python3 excute.py -f PDB_file
```
#### Example 1: Using PDB Structure as Input
To predict metal binding sites for a specific PDB structure (e.g., `2zp9`), use the `-p` option:
```bash
python3 execute.py -p 2zp9
```
The results will be saved in the `output_data` directory, with the predicted metal information stored in `2zp9_metal.pdb` and detailed binding site data in `2zp9_output.csv`.

#### Example 2: Using Your Own Uploaded File
To use your own PDB file (e.g., `7PW5.pdb`), place it in the `input_data` directory and run:
```bash
python3 execute.py -f 7PW5.pdb
```
The output files `7PW5_metal.pdb` and `7PW5_output.csv` will be saved in the `output_data` directory.

#### Output Files Description
- **`X_metal.pdb`**: This file contains the predicted metal binding information. If any metal binding sites are predicted, they will appear at the end of the file. The insertion code in column 27 will be marked with "@" to indicate the metals predicted by PMM. Corresponding ligands can be found in the LINK section.
- **`X_output.csv`**: This file stores detailed information about the predicted binding sites.

You can generally use only the `X_metal.pdb` file to analyze the prediction results. The `X_output.csv` file provides additional details for further analysis.


## File Descriptions

- **`conda_dockerfile`**: Dockerfile for building the Conda container with all dependencies for PinMyMetal.
- **`docker-compose.yml`**: Configuration file to set up and manage the Conda and PostgreSQL containers.
- **`execute/`**: Directory containing the input files, execution scripts, and output files for the PinMyMetal project.


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

