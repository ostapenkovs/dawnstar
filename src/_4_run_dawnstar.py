# Dawnstar, a patent expiration workflow.
#
# Programmatically find projected expiration dates 
# and other assumed status details for compound patents.
#
#   Requiring Python == 3.11.

from collections.abc import Generator
import argparse
import os
from tqdm import tqdm
import time
from datetime import datetime

import psycopg2
import aiohttp
import asyncio
import json
import bs4

from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*')

from config import db_info

DATA_FOLDER = './data'
BASE = 'https://patents.google.com/patent/'
CHUNK_SIZE = 50
SLEEP_TIME = 3
HEADERS = {'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) \
           AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36'}

def chunks(l: list, n: int) -> Generator:
    '''Break list l into chunks of size n to iterate over.'''
    for i in range(0, len(l), n):
        yield l[i: i+n]

def get_compounds(infile: str, canonical: bool) -> list:
    """Process user input-file and return compound SMILES after sanitizing."""
    compounds = []
    if infile.split('.')[-1] == 'sdf':
        all_mols = Chem.SDMolSupplier(infile)
    else:
        with open(infile, 'r') as f:
            all_mols = [Chem.MolFromSmiles(line.strip()) for line in f]
    for mol in all_mols:
        if mol:
            if canonical:
                smi = Chem.MolToSmiles(mol, canonical=True)
            else:
                smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True)
            if smi:
                compounds.append(smi)
    return compounds

def validate_date(date: str) -> bool:
    """Validate potential before and after date cutoffs.
    Want date to be before NOW and after Jan. 1st, 1962."""
    try:
        if len(date) == 8:
            yr, mo, dy = int(date[0:4]), int(date[4:6]), int(date[6:8])
            date = datetime(year=yr, month=mo, day=dy)
            if date >= datetime(year=1962, month=1, day=1) and date < datetime.now():
                return True
    except Exception as e:
        print(e)
    return False

def safe_select(event: bs4.element.Tag, item: str) -> str | int:
    """Process HTML tags in patent timeline to ensure accurate scraping."""
    try:
        val = event.select_one(item).text.strip()
        val = None if not val else val
    except:
        val = None
    return val

def build_query(compound: str, is_smiles: bool, search_type: str, base_filters: bool, after: None | str, before: None | str) -> tuple:
    """Build two query URLs for a given compound, one for search and one for XHR download.
    Returns tuple of length two with each URL."""
    if is_smiles:
        to_replace = {'%': r'%25', '+': r'%2b', '=': r'%3d', '/': r'%2f', '#': r'%23'}
        for ch in to_replace.keys():
            compound = compound.replace(ch, to_replace[ch])
    
    eq1 = r'%3d'
    sim_sub = ''
    if search_type == 'substructure':
        sim_sub = f'SSS{eq1}'
    elif search_type == 'similarity':
        sim_sub = '~'
    
    base_filter = ''
    if base_filters:
        base_filter = '&country=US&language=ENGLISH&type=PATENT'
    
    after_filter, before_filter = '', ''
    if after:
        after_filter = f'&after=priority:{after}'
    if before:
        before_filter = f'&before=priority:{before}'
    
    query = f'q=CL{eq1}{sim_sub}({compound}){base_filter}{after_filter}{before_filter}'

    search_url = 'https://patents.google.com/?' + query

    eq2, amp, col, pct3d, b1, b2, pct2b = r'%3D', r'%26', r'%3A', r'%253d', r'%5B', r'%5D', r'%252b'
    query = query.replace(eq1, pct3d).replace('=', eq2).replace('&', amp).replace(':', col).replace('[', b1).replace(']', b2).replace(r'%2b', pct2b)
    download_url = f'https://patents.google.com/xhr/query?url=' + query + '&exp=&download=true'

    return search_url, download_url

def parse_timeline(page_text: str) -> tuple:
    """Process webpage text of a single patent for a given compound.
    Extract timeline events, in particular status and expiration."""
    title = date = status = status_detail = None
    soup = bs4.BeautifulSoup(page_text, 'html.parser')
    title, events = soup.find('title').text.split()[0], soup.find_all('dd', {'itemprop': "events"})
    soup = None
    for event in events:
        e_type = safe_select(event, 'span[itemprop="type"]')
        if e_type is not None and e_type == 'legal-status':
            e_title = safe_select(event, 'span[itemprop="title"]')
            e_date = safe_select(event, 'time[itemprop="date"]')
            if e_date is not None and e_date == 'Status':
                status = e_title
            elif e_date is not None:
                status_detail = e_title
                date = e_date
    return title, date, status, status_detail

async def get_url_async(session: aiohttp.ClientSession, url: str, headers: dict, base=None) -> str:
    """Makes async request to get single patent webpage text.
    Raises error if response is not 200."""
    if base is not None: url = base + url
    async with session.get(url=url, headers=headers) as response:
        response.raise_for_status()
        return await response.text()

async def get_timeline_async(session: aiohttp.ClientSession, url: str, headers: dict, base=None, pbar=None) -> tuple | None:
    """Async wrapper for parse_timeline function."""
    page_text = await get_url_async(session=session, url=url, headers=headers, base=base)
    title, date, status, status_detail = parse_timeline(page_text=page_text)
    if pbar is not None: pbar.update(1)
    return (title, (date, status, status_detail))

async def get_all_async(urls: list, headers: dict, base=None, pbar=None) -> list:
    """Collect information for all patents in async fashion."""
    async with aiohttp.ClientSession() as session:
        return await asyncio.gather(*[get_timeline_async(session=session, url=url, headers=headers, base=base, pbar=pbar) for url in urls])

def patents_for_compound(curs, compound: str, search_type: str, threshold=None, after=None, before=None) -> list | None:
    '''Find all patents for a given compound from SureChembl database with optional filters.'''
    data = None
    ###
    if search_type == 'exact':
        search = 'c.smiles @= %s '
    elif search_type == 'substructure':
        search = 'c.smiles @> %s '
    elif search_type == 'similarity':
        m = Chem.MolFromSmiles(compound)
        if not m:
            return
        compound = DataStructs.BitVectToBinaryText(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024))
        m = None
        search = 'c.mfp%%bfp_from_binary_text(%s) '
        if not threshold:
            threshold = 0.75
        curs.execute('set rdkit.tanimoto_threshold = %s', (threshold, ))
    else:
        return
    ###
    sql = 'select p.num from compound c, patent p, field_freq f \
        where c.id = f.compound_id and p.id = f.patent_id and f.field_id = 2 and ' + search
    if after:
        sql += f"and p.pub_date >= {after}::text::timestamp "
    if before:
        sql += f"and p.pub_date < {before}::text::timestamp "
    ###
    curs.execute(sql, (compound, ))
    try:
        data = [x[0] for x in curs.fetchall()]
    except Exception as e:
        print(e)
    return data

def main() -> None:
    """Process command line arguments and execute patent expiration workflow."""
    ### BEGIN USER ARGS ###
    parser = argparse.ArgumentParser(description='Dawnstar patent expiration workflow command line utility.')

    parser.add_argument('-if', '--infile', type=str, help='FilePATH (!) containing compounds in SDF or SMILES format.')
    parser.add_argument('-of', '--outfile', type=str, help='Desired fileNAME (!) of resulting output.')
    parser.add_argument('-st', '--searchtype', type=str, help='What type of search to perform (exact, similarity, substructure).', required=False)
    parser.add_argument('-th', '--threshold', type=float, help='If similarity search, specify the threshold (between 0 and 1).', required=False)
    parser.add_argument('-a', '--after', type=str, help='After (lower-bound) date in yyyymmdd format', required=False)
    parser.add_argument('-b', '--before', type=str, help='Before (upper-bound) date in yyyymmdd format', required=False)
    
    args = parser.parse_args()

    if (not os.path.isfile(args.infile)) or (args.infile.split('.')[-1] not in ['txt', 'smi', 'sdf']):
        print('Unacceptable infile.')
        return
    
    search_type = 'exact'
    if args.searchtype is not None:
        if args.searchtype in ['exact', 'similarity', 'substructure']:
            search_type = args.searchtype
        else:
            print('Unacceptable searchtype.')
            return
    
    threshold = 0.75
    if args.threshold is not None:
        if args.threshold > 0 and args.threshold < 1:
            threshold = args.threshold
        else:
            print('Unacceptable threshold.')
            return
    
    after = None
    if args.after is not None:
        if not validate_date(date=args.after):
            print('Unacceptable after.')
            return
        else:
            after = args.after

    before = None
    if args.before is not None:
        if not validate_date(date=args.before):
            print('Unacceptable before.')
            return
        else:
            before = args.before
    ### END USER ARGS ###

    ### BEGIN GET COMPOUNDS ###
    compounds = get_compounds(infile=args.infile, canonical=True)
    if not compounds:
        print('No compounds found.')
        return
    else:
        print(f'Will process {len(compounds)} compound(s). BEGIN.\n')
    ### END GET COMPOUNDS
    
    ### BEGIN DAWNSTAR WORKFLOW ###
    schema_name = 'sc'

    conn = psycopg2.connect(**db_info)

    curs = conn.cursor()
    curs.execute(f'set search_path to {schema_name}')

    ###
    data = {}
    for compound in tqdm(compounds, desc='Compounds'):
        patents = patents_for_compound(curs=curs, compound=compound, search_type=search_type,
                                       threshold=threshold, after=after, before=before)
        if not patents: continue
        
        inner_pbar = tqdm(total=len(patents), leave=False, desc='Patents')
        results = []
        for pat in chunks(l=patents, n=CHUNK_SIZE):
            res = asyncio.run(get_all_async(urls=pat, headers=HEADERS, base=BASE, pbar=inner_pbar))
            results += res
            res = None
            for _ in tqdm(range(SLEEP_TIME), leave=False, desc='Sleeping'): time.sleep(1)
        inner_pbar.close()
        
        if not results: continue

        data[compound] = {k: v for k, v in results}
        patents = results = None
    ###

    curs.close()
    conn.close()
    ### END DAWNSTAR WORKFLOW ###

    ### BEGIN DATA DUMP ###
    if not data:
        print('Nothing found for given compounds.')
        return

    print('Dumping found data to JSON.')
    if not os.path.exists(DATA_FOLDER):
        os.makedirs(DATA_FOLDER)
    with open(f'{DATA_FOLDER}/{args.outfile}.txt', 'w') as f:
        json.dump(data, f, ensure_ascii=False)
    ### END DATA DUMP ###

if __name__ == '__main__':
    main()
