import psycopg2
from psycopg2.extensions import AsIs

import glob
import pandas as pd
import io

import time

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, SaltRemover
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from config import db_info

def main() -> None:
    '''Populate database within Docker container with our SureChembl files.'''
    
    schema_name = 'sc'

    conn = psycopg2.connect(**db_info)
    conn.autocommit = True

    curs = conn.cursor()
    curs.execute(f'set search_path to {schema_name}')

    data_folder = '../sure_chembl_data/'
    files = glob.glob(data_folder + '*.txt.gz')

    colnames = ['id', 'smiles', 'inchi_key', 'corpus_freq', 'patent_id', \
                'pub_date', 'field_type', 'field_freq']

    compound = 'compound'
    cmp_cols = ['id', 'smiles']

    patent = 'patent'
    pat_cols = ['patent_id', 'pub_date']

    field_freq = 'field_freq'
    fld_cols = ['id', 'patent_id', 'field_type', 'field_freq']

    remover = SaltRemover.SaltRemover()

    ### POPULATING THE DATABASE ###
    for f_idx, f in enumerate(files):
        print(f'Working on file number {f_idx + 1} out of {len(files)}.')

        ###
        # read one file
        df = pd.read_csv(f, header=None, sep='\t')
        df.columns = colnames

        # filter the file
        df = df.drop(['inchi_key', 'corpus_freq'], axis=1)
        df['id'] = df.id.apply(lambda x: x.split('SCHEMBL')[1]).astype(int)
        df['patent_id'] = df.patent_id.apply(lambda x: ''.join(x.split('-')))
        df = df[df.patent_id.str[0:2] == 'US']
        ###

        ###
        container = dict()
        for tbl, cols in zip([compound, patent, field_freq], [cmp_cols, pat_cols, fld_cols]):
            container[tbl] = df[cols].drop_duplicates(cols[:-1])
        df = None
        ###

        ###
        # upserting compounds
        start = time.time()
        invalid_compounds = []
        i_idx = 0
        for id, smi in zip(container[compound].id, container[compound].smiles):
            if i_idx % 100000 == 0 and i_idx > 0:
                print(f'Processed {i_idx} compounds out of {len(container[compound])}.')
            
            m = Chem.MolFromSmiles(smi)
            if not m: 
                invalid_compounds.append(id)
                continue
            m = remover.StripMol(m)
            can_smi = Chem.MolToSmiles(m, canonical=True)
            if not can_smi: 
                invalid_compounds.append(id)
                continue
            b_mfp = DataStructs.BitVectToBinaryText( AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024) )

            try:
                curs.execute('insert into %(tbl)s values (%(id)s, %(smi)s::mol, bfp_from_binary_text(%(mfp)s)) \
                    on conflict (id) do update set smiles = excluded.smiles, mfp = excluded.mfp',
                    {'tbl': AsIs(compound), 'id': id, 'smi': can_smi, 'mfp': b_mfp})
                
            except Exception as e:
                print(e)
                invalid_compounds.append(id)
            
            m = can_smi = b_mfp = None
            i_idx += 1
        
        del container[compound]

        if len(invalid_compounds) >= 1:
            container[field_freq] = container[field_freq][ ~(container[field_freq].id.isin(invalid_compounds)) ]
            container[patent] = container[patent][ (container[patent].patent_id.isin( container[field_freq].patent_id.unique().tolist() )) ]
        invalid_compounds = None
        
        print(f'Compounds took {time.time() - start} seconds.')
        ###
        
        ###
        # upserting patents
        start = time.time()
        curs.execute(f'create table tmp_{patent} (like {patent} including defaults)')

        buf = io.StringIO(container[patent].to_csv(index=False, header=True, sep=','))
        del container[patent]
        curs.copy_expert(sql=f"copy tmp_{patent} (num, pub_date) from stdin with csv header delimiter as ','", file=buf)
        buf = None

        curs.execute(f'insert into {patent} select * from tmp_{patent} on conflict (num) \
                    do update set pub_date = excluded.pub_date returning num, id')
        
        # need the patent map!
        patent_map = {k: v for k, v in curs.fetchall()}

        curs.execute(f'drop table tmp_{patent}')
        container[field_freq]['patent_id'] = container[field_freq].patent_id.map(patent_map)
        patent_map = None

        print(f'Patents took {time.time() - start} seconds.')
        ###

        ###
        # upserting field frequencies
        start = time.time()

        curs.execute(f'create table tmp_{field_freq} (like {field_freq} including defaults)')

        buf = io.StringIO(container[field_freq].to_csv(index=False, header=True, sep=','))
        del container[field_freq]
        curs.copy_expert(sql=f"copy tmp_{field_freq} from stdin with csv header delimiter as ','", file=buf)
        buf = None

        curs.execute(f'insert into {field_freq} select * from tmp_{field_freq} on conflict (compound_id, patent_id, field_id) \
                    do update set freq = excluded.freq')
        curs.execute(f'drop table tmp_{field_freq}')

        print(f'Field freqs took {time.time() - start} seconds.')
        ###

        print(f'Processed file number {f_idx + 1} out of {len(files)}.')
    ###

    curs.execute('create index smi_idx on compound using gist(smiles)')
    curs.execute('create index mfp_idx on compound using gist(mfp)')

    curs.close()
    conn.close()

if __name__ == '__main__':
    main()
