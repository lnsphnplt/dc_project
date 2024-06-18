# Digital Chemistry Project
## Getting started

Prerequisits:
* Python
* Docker [[https://www.docker.com/products/docker-desktop/]]

1. Start docker engine by starting docker desktop 
2. navigate to the repo root directory and execute "python startVirtualEnv.py"
3. A browser with jupyter lab opens. Enjoy!

PS: The first time you start it might take a while because docker needs to build the necessairy images.

## Structure

* 00-Utils
    * Contains tools used for setting up the virtual environements
    * all docker files

* 01-Data_Extraction
    * Scraping of html databases
    * Uniform formating of all accessible datasets
    * Matching of pos/neg samples in sequence-length distributions

* 02-Data_Curation
    * Checks data integrity
    * converts csv files into unfied sqlite database
    * checks for contradictory data

* 03-Data_exploration
    * analyize data

* 04-Data_processing
    * provides the scripts necessairy to calculate the descriptors localy or on euler cluster

## Snakemake
For optimal reproducablity this entire project can be run via snakemake. However there are some caviats and things to adress:

* Make sure you are using the newest snakemake version:

```bash
(snakemake) residual@XPS:/mnt/z/projects/dc_project$ snakemake --version
8.11.6
```
Older versions like the one you get wenn using:
```bash
(base) residual@XPS:/mnt/z/projects/dc_project$ sudo apt install snakemake
(base) residual@XPS:/mnt/z/projects/dc_project$ snakemake --version
6.15.1
```
do not work. There are recent changes in the snakemake/docker syntax. Make sure you are using at least verion 8. I recommend installing snakemake as recommended in their docks using mamba: [Installtion guide](https://snakemake.readthedocs.io/en/stable/getting_started/installation.html) 

* Make sure apptainer is installed:
Sadly apptainer does not support windows. We run this on a windows machine using WSL. Make sure you have the correct version

```bash
(snakemake) residual@XPS:/mnt/z/projects/dc_project$ apptainer --version
apptainer version 1.3.1
```

* Lastly make sure docker is installed and the docker engine is running.
```bash
(snakemake) residual@XPS:/mnt/z/projects/dc_project$ docker --version
Docker version 26.0.0, build 2ae903e
```

* Adjust the absult path of the root directory in 

```path
dc_project/workflow/config.yaml

```
### Execution

For the execution of the project we recommend the following command:

```bash
snakemake --sdm apptainer -c all
```

This executes the entire project using all avalible cores.  Depending on your internet speed it takes a moment for the docker image to be downloaded.