Workflow to programatically retrieve patent information by SMILES/SDF, allowing for exact and inexact matches.

EXAMPLE USAGE:
1. docker build -t my-postgres .
2. docker compose up

3. create secret.py with dbinfo dictionary (containing db connection params)
4. run all cells in populate_db.ipynb

5. execute "python run_dawnstar.py -if ./data/metformin.txt -is y -a 20210101 -b 20230101 -of result"
6. run all cells in "analyze_dawnstar.ipynb" to inspect ./data/result.txt for information to analyze

Requires Python == 3.11
