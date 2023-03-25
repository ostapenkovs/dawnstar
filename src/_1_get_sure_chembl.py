import os
import time
from tqdm import tqdm
import ftplib

def main() -> None:
    '''Download all SureChembl data files from Chembl FTP'''
    
    start = time.time()

    data_folder = '../sure_chembl_data'
    if not os.path.exists(data_folder):
        os.makedirs(data_folder)
    
    server = 'ftp.ebi.ac.uk'
    datadir = 'pub/databases/chembl/SureChEMBL/data/map'

    session = ftplib.FTP(server)

    session.login()
    session.cwd(datadir)

    files = [f for f in session.nlst() if '.txt.gz' in f]
    files.reverse()

    for file in tqdm(files):
        with open(f'{data_folder}/{file}', 'wb') as f:
            session.retrbinary('RETR ' + file, f.write)

    session.quit()

    print(f'File download took {time.time() - start} seconds.')

if __name__ == '__main__':
    main()
