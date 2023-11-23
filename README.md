# cgeisolate

# Installation

The following Conda channels are required:

`conda config --add channels defaults`

`conda config --add channels bioconda`

`conda config --add channels conda-forge`

Mamba can be installed with:

`conda install -c conda-forge mamba`

For a fast install of cgeisolate, use mamba:
`mamba install -c genomicepidemiology cgeisolate`


# Datebase

Download the cge_db database:

`wget https://cge.food.dtu.dk/services/MINTyper/cge_db.tar.gz`

`tar -xvzf cge_db.tar.gz`

`sudo mkdir -m 777 /opt/cge`

`mv cge_db /opt/cge/cge_db`


# Usage

Standard Usage:
`cgeisolate -i <input_file> -o <output_file>`

If you have a folder of many fastq.gz files:
`cgeisolate -f <input_folder> -name <sample_name> -o <output_folder>`

If your cge_db is not stored in /opt/cge/cge_db:
`cgeisolate -i <input_file> -o <output_file> -db_dir <path_to_cge_db>`

Help Message:
`cgeisolate -h`
