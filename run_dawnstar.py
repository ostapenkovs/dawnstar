# Dawnstar, a patent expiration workflow.
#
# Programmatically find projected expiration dates 
# and other assumed status details for compound patents.
#
#   Requiring Python == 3.11.

import requests
from io import StringIO
import aiohttp
import asyncio
import os
import sys
import time
from datetime import datetime
import argparse
import json
import pandas as pd
from collections import defaultdict
import bs4
from rdkit import Chem
import random

DATA_FOLDER = './compound_data'
HEADERS = {'user-agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) \
           AppleWebKit/537.36 (KHTML, like Gecko) Chrome/110.0.0.0 Safari/537.36'}

def get_compounds(infile: str) -> list:
    """Process user input-file and return compound SMILES after sanitizing."""
    compounds = []
    if infile.split('.')[-1] == 'sdf':
        all_mols = Chem.SDMolSupplier(infile)
    else:
        with open(infile, 'r') as f:
            all_mols = [Chem.MolFromSmiles(line.strip()) for line in f]
    for mol in all_mols:
        if mol:
            smi = Chem.MolToSmiles(mol, isomericSmiles=True, kekuleSmiles=True)
            if smi:
                compounds.append(smi)
    return compounds

def validate_date(date: str) -> bool:
    """Validate potential before and after date cutoffs.
    Want date to be before NOW and after Jan. 1st, 1970."""
    try:
        if len(date) == 8:
            yr, mo, dy = int(date[0:4]), int(date[4:6]), int(date[6:8])
            date = datetime(year=yr, month=mo, day=dy)
            if date >= datetime(year=1970, month=1, day=1) and date < datetime.now():
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

async def get_url_async(session: aiohttp.ClientSession, url: str, headers: dict) -> str:
    """Makes async request to get single patent webpage text.
    Raises error if response is not 200."""
    async with session.get(url=url, headers=headers) as response:
        response.raise_for_status()
        return await response.text()

async def get_timeline_async(session: aiohttp.ClientSession, url: str, headers: dict) -> tuple | None:
    """Async wrapper for parse_timeline function."""
    page_text = await get_url_async(session=session, url=url, headers=headers)
    title, date, status, status_detail = parse_timeline(page_text=page_text)
    if title and date and (status or status_detail):
        return (title, (date, status, status_detail))
    return None

async def get_all_async(urls: list, headers: dict) -> list:
    """Collect information for all patents in async fashion."""
    l = []
    async with aiohttp.ClientSession() as session:
        for i, url in enumerate(urls):
            l.append(get_timeline_async(session=session, url=url, headers=headers))
            if i % 5 == 0 and i > 0:
                print(f'Processed URL {i}.')
        return await asyncio.gather(*l)

def main() -> None:
    """Process command line arguments and execute patent expiration workflow."""
    ### BEGIN USER ARGS ###
    parser = argparse.ArgumentParser(description='Dawnstar patent expiration workflow command line utility.')

    parser.add_argument('-if', '--infile', type=str, help='FilePATH (!) containing compounds in SDF or SMILES/INCHIKEY format.')
    parser.add_argument('-is', '--issmiles', type=str, help='If SMI format, is it SMILES (versus INCHIKEY) [y/n]?', required=False)
    parser.add_argument('-st', '--searchtype', type=str, help='What type of search to perform (exact, similarity, substructure).', required=False)
    parser.add_argument('-a', '--after', type=str, help='After (lower-bound) date in yyyymmdd format', required=False)
    parser.add_argument('-b', '--before', type=str, help='Before (upper-bound) date in yyyymmdd format', required=False)
    parser.add_argument('-of', '--outfile', type=str, help='Desired fileNAME (!) of resulting output.')
    
    args = parser.parse_args()

    if (not os.path.isfile(args.infile)) or (args.infile.split('.')[-1] not in ['txt', 'smi', 'sdf']):
        print('Unacceptable infile.')
        return
    
    is_smiles = True
    if args.issmiles is not None:
        if args.issmiles == 'n':
            is_smiles = False
        elif args.issmiles != 'y':
            print('Unacceptable issmiles.')
            return
    
    search_type = 'exact'
    if args.searchtype is not None:
        if args.searchtype in ['exact', 'similarity', 'substructure']:
            search_type = args.searchtype
        else:
            print('Unacceptable searchtype.')
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
    compounds = get_compounds(infile=args.infile)
    if not compounds:
        print('No compounds found.')
        return
    else:
        print(f'Will process {len(compounds)} compound(s). BEGIN.\n')
    ### END GET COMPOUNDS
    
    ### BEGIN DAWNSTAR WORKFLOW ###
    reg_session = requests.Session()

    data = defaultdict(dict)
    for i, comp in enumerate(compounds):
        print(f'Working on compound: {comp}...')
        search_url, download_url = build_query(compound=comp, is_smiles=is_smiles, search_type=search_type,
                                               base_filters=True, after=after, before=before)
        
        print(f'Searching URL: {search_url}')

        ### BEGIN GOOGLE XHR ###
        r = reg_session.get('http://google.com')
        r.raise_for_status()
        cookies = reg_session.cookies.get_dict()

        j, max_iter = 0, 10
        start_sec, inc_sec = 5, 10
        while True:
            try:
                r = reg_session.get(download_url, headers=HEADERS, cookies=cookies)
                r.raise_for_status()
                break
            except requests.exceptions.HTTPError:
                if j > max_iter:
                    raise ValueError('Too much sleeping.')
                elif j > 0:
                    start_sec += inc_sec
                print(f'Google XHR error. Sleeping for {start_sec} seconds...')
                time.sleep(start_sec)
                j += 1
        ### END GOOGLE XHR ###

        patent_urls = None
        patent_data = pd.read_csv(StringIO(r.text), skiprows=1)
        if len(patent_data) >= 1:
            patent_urls = patent_data['result link'].dropna(axis=0, how='any').unique().tolist()
        # patent_urls = ['https://patents.google.com/patent/US6421675B1/']
        
        results = None
        if patent_urls:
            print(f'Found {len(patent_urls)} patents for this compound.')
            results = asyncio.run(get_all_async(urls=patent_urls, headers=HEADERS))
            results = [res for res in results if res is not None]
        
        if results:
            data[comp].update(results)
        
        reg_session.cookies.clear()
        
        if len(compounds) > 1:
            print('Sleeping for a random number of seconds before moving on...')
            time.sleep( random.randint(5, 10) )

    reg_session.close()
    ### END DAWNSTAR WORKFLOW ###

    ### BEGIN DATA DUMP ###
    if data:
        print('Dumping found data to JSON.')
        if not os.path.exists(DATA_FOLDER):
            os.makedirs(DATA_FOLDER)
        with open(f'{DATA_FOLDER}/{args.outfile}.txt', 'w') as f:
            json.dump(data, f, ensure_ascii=False)
    ### END DATA DUMP ###

if __name__ == '__main__':
    main()
