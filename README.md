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


