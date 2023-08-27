## Workflow to programatically retrieve patent information by SMILES/SDF, allowing for exact and inexact matches.

### Instructions presume user has Anaconda/Miniconda package manager installed on their machine. Requires Python == 3.11.

*DOCKER (using files inside ./Docker directory)*
1. install Docker Desktop and WSL 1
2. build Docker image: "docker build -t my-postgres ."
3. start Docker container and initialize database: "docker compose up"
    * run Docker container in the background: add "-d" flag
    * turn off Docker container: "docker compose down"
        * destroy containers and volumes: add "-v" flag

*PYTHON (using files inside ./src directory)*
1. conda create -n "dawnenv" python=3.11
2. conda activate dawnenv
3. "pip install -r requirements.txt"
4. "python _1_get_sure_chembl.py"
5. "python _2_split_backfile.py"
6. "python _3_populate_db.py"
7. DAWNSTAR: "python _4_run_dawnstar.py"

*DAWNSTAR*
* example usage: "python _4_run_dawnstar.py -if ../compound_data/metformin.txt -of result -st exact"
    * get argument help: "python _4_run_dawnstar.py --help"

*MISC*
* default Postgres database connection parameters: inside ./Docker/.env and ./src/config.py
* Docker database files location on local machine: "\\\wsl.localhost\docker-desktop-data\data\docker\volumes\\{volume_name}\\_data"
    * replace {volume_name} with name of your volume
