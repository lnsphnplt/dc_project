# Digital Chemistry Project
## numpy version 2.0.0
On the 16.6 numpy released version 2 this breaks rdkit. We recomend using an older numpy version. <2

```bash
pip install --force-reinstall -v "numpy==1.26.4"
```
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

### Setup
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

To reproduce our project we recommend the following steps:

1. Data extraction


It is theoretically possible to run the webscraping in a docker environment, however everytime we tried to implement this we had a lot of problems with a) the virtualisation of a headless browser in docker and b) Spam/bot detection. This is why we recommend running this locally so your browerser is going to have a more unique fingerprint and avoids the bot detection. To do this:

```bash
pip install 00-Utils/Docker-ipnyb-execution/requirements.txt
jupyter nbconvert --execute 01-Data_Extraction/scraping.ipynb --to notebook
jupyter nbconvert --execute 01-Data_Extraction/data_extraction.ipynb --to notebook
jupyter nbconvert --execute 01-Data_Extraction/unified_data.ipynb --to notebook
```

To skip this step:
```bash
mv 01-Data_Extraction/new_data_pos_neg_backup.csv 01-Data_Extraction/new_data_pos_neg.csv
```

2. Create database and validate data

```bash
snakemake --until EDA --sdm apptainer -c all
```

This only executes the first part of the project. Depending on your internet speed it takes a moment for the docker image to be downloaded. The image only needs to downloaded once.

To skip this step:
```bash
mv 02-Data_Curation/unified-curated_backup.db 02-Data_Curation/unified-curated.db
```

3. Calculate descriptors

You can let this run automatically in snakemake however we dont recomend this. This is a very performance intensive task. We have optimized the code to make use of multithreading. Running this in a docker container would waste computing power by the virtualization used by docker. This calculation will max out your CPU, making your PC unusable for the next 15-60 min (Depending on the performance of your PC). We recommend letting this run on the euler cluster. 


To skip this step:
```bash
mv unified-final_backup.db unified-final.db
```

### Calculate on Euler
SSH into euler (you need to be in a eth network or use the VPN)
```bash
ssh username@euler.ethz.ch
```

We recommend working the scratch directory. This will be automatically delted after some time. This ensures you dont clutter up euler. Replace username with your eth creds.

```bash
cd $SCRATCH
```
__On your local machine:__

Copy the relevant files to euler(replace username with your username):
scripts:
```bash
scp 04-Data_calculate_Descriptors/copyToEuler/* username@euler.ethz.ch:/cluster/scratch/username/
```
Data file:
```bash
scp 02-Data_Curation/unified-curated.db username@euler.ethz.ch:/cluster/scratch/username/
```

__On euler:__

Now this should look like this:
```bash
[emathier@eu-login-48 emathier]$ pwd
/cluster/scratch/emathier
[emathier@eu-login-48 emathier]$ ls
CD2.py  runOnEuler.sh  unified-curated.db
```

Submit into the batch system:
```bash
sbatch runOnEuler.sh
```

Monitor resource usage
```bash
mjobs

Job information
 Job ID                          : 62516304
 Job name                        : calculate_descriptors
 Status                          : RUNNING
 Running on node                 : eu-g9-024-2
 User                            : emathier
 Shareholder group               : es_chab
 Slurm partition (queue)         : normal.4h
 Command                         : sbatch runOnEuler.sh
 Working directory               : /cluster/scratch/emathier
Requested resources
 Requested runtime               : 00:25:00
 Requested cores (total)         : 48
 Requested nodes                 : 1
 Requested memory (total)        : 14400 MiB
Job history
 Submitted at                    : 2024-06-18T16:35:45
 Started at                      : 2024-06-18T16:35:55
 Queue waiting time              : 10 s
Resource usage
 Wall-clock                      : 00:02:02
 Total CPU time                  : 01:23:52
 CPU utilization                 : 85.92%
 Total resident memory           : 5452.34 MiB
 Resident memory utilization     : 37.86%
```

Monitor progress:
```bash
tail -f output.log
Processing:  62%|██████▏   | 7969/12828 [04:33<04:32, 17.84it/s]
```

#### Result:

```bash
Job information
 Job ID                          : 62516304
 Job name                        : calculate_descriptors
 Status                          : COMPLETED
 Running on node                 : eu-g9-024-2
 User                            : emathier
 Shareholder group               : es_chab
 Slurm partition (queue)         : normal.4h
 Command                         : sbatch runOnEuler.sh
 Working directory               : /cluster/scratch/emathier
Requested resources
 Requested runtime               : 00:25:00
 Requested cores (total)         : 48
 Requested nodes                 : 1
 Requested memory (total)        : 14400 MiB
Job history
 Submitted at                    : 2024-06-18T16:35:45
 Started at                      : 2024-06-18T16:35:55
 Queue waiting time              : 10 s
Resource usage
 Wall-clock                      : 00:11:28
 Total CPU time                  : 06:52:08
 CPU utilization                 : 74.87%
 Total resident memory           : 5709.97 MiB
 Resident memory utilization     : 39.65%
```
Looking at the execution stats it took around 11 min to complete. If we look at Total CPU time we see that it would take around 7 hours on a single core.

Finally copy back the database file:

On your local machine
```bash
scp username@euler.ethz.ch:/cluster/scratch/username/unified-curated.db .
```

### Local execution

This is not recommend however if you want to run it locally: Execute the following commands from the root directory. (Maybe do this in a venv)
Install dependencies
 ```bash
sudo pip install -r 00-Utils/Docker-ipnyb-execution/requirements.txt
python3 04-Data_calculate_Descriptors/CD2.py
 ```
