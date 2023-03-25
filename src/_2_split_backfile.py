import os
import time
import glob
import gzip

def main() -> None:
    '''Split SureChembl backfile (first, largest file) into multiple for easier reading.'''
    
    start = time.time()
    
    data_folder = '../sure_chembl_data/'
    file = glob.glob(data_folder + '*.txt.gz')[0]

    lpf = 2.5e7
    infile = gzip.open(file, 'rt')
    outfile = None

    j = 0
    for i, line in enumerate(infile):
        if i % lpf == 0:
            j += 1
            if outfile:
                outfile.close()
            
            fname = file[:-7] + f'_{j}.txt.gz'
            outfile = gzip.open(fname, 'wt')
        
        outfile.write(line)

    if outfile:
        outfile.close()

    infile.close()

    os.remove(file)

    print(f'File split took {time.time() - start} seconds.')

if __name__ == '__main__':
    main()
