import numpy as np 
import h5py
from scipy.sparse import csc_matrix
import pandas as pd
import sys

class base_extractor(object):
    def __init__(self, filename, **kwargs):
        self.filename = filename

    def __getitem__(self, i):
        raise NotImplementedError


def process_chrom_sizes(chromsizes_path):
    df = pd.read_table(chromsizes_path).set_index('chrom').sort_index()
    df['cumsum_size'] = df['size'].cumsum().shift(1).fillna(0)
    return df


class FiberSeqExtractor(base_extractor):

    def __init__(self, filename, chromsizes_path, **kwargs):
        super(FiberSeqExtractor, self).__init__(filename, **kwargs)
        self.chromsizes_df = process_chrom_sizes(chromsizes_path)
        self.h5 = h5py.File(filename, 'r')
    
    def __getitem__(self, interval):
        ids_range = self.interval_to_ids(interval)
        # implement any post processing steps here
        return self.read_sparse_matrix(ids_range)

    def __del__(self):
        self.close()

    def close(self):
        if self.h5:
            del self.h5

    def read_sparse_matrix(self, col_indices_range=None):
        """
        Read a sparse matrix from h5 file and return a subset of the columns.

        Parameters:
            col_indices (list): List of column indices to select.

        returns:
            csc_matrix: Sparse matrix with the specified columns.
        """
        if col_indices_range is None:
            raise NotImplementedError("Please specify the column indices to read.")

        col_start, col_end = col_indices_range
        
        indptr_slice = self.h5['indptr'][col_start:col_end + 1]
        data_start = indptr_slice[0]
        data_end = indptr_slice[-1]

        # Adjust indices and data slices according to the new indptr
        new_indptr =  indptr_slice - data_start
        new_data = self.h5['data'][data_start:data_end]
        new_indices = self.h5['indices'][data_start:data_end]
        new_indices -= new_indices.min()
        new_shape = (max(new_indices) + 1, len(new_indptr) - 1)
    
        return csc_matrix((new_data, new_indices, new_indptr), shape=new_shape)

    
    def interval_to_ids(self, interval):
        """
        Convert genomic interval to index in the genome.

        Parameters:
            interval (GenomicInterval): Genomic interval (chr, start, end)
        
        Returns:
            tuple: Start and end index of the interval in the genome.
        """

        chrom_start_index = self.chromsizes_df.loc[interval.chrom, 'cumsum_size']
        start_index = chrom_start_index + interval.start
        end_index = chrom_start_index + interval.end
        return start_index, end_index


from genome_tools.data.extractors import fasta_extractor as FastaExtractor
from genome_tools.genomic_interval import genomic_interval as GenomicInterval
import re
from scipy.sparse import coo_matrix, save_npz
from tqdm import tqdm
import argparse
import pandas as pd
import numpy as np
from save_to_h5 import save_to_h5py

def process_chrom_sizes(chromsizes_path):
    df = pd.read_table(chromsizes_path, names=['chrom', 'size']).set_index('chrom').sort_index()
    df['cumsum_size'] = df['size'].cumsum().shift(1).fillna(0)
    return df


def create_coo_from_bed(bed_df, chromsizes_df):
    """
    Create a h5 file from a bed file.

    Parameters:
        bed_df (pd.DataFrame): Bed file with the following columns: chrom, start, end, ref_m6a, strand.
        chromsizes_path (str): Path to the chromsizes file.
        h5_path (str): Path to the h5 file.
    """
    
    incorrect_positions = set(["", "-1", "."])
    
    bed_df['chrom_index'] = bed_df['chrom'].map(chromsizes_df['cumsum_size'])
    bed_df['start_index'] = bed_df['chrom_index'] + bed_df['start']
    
    coo_data = []
    coo_row = []
    coo_col = []
    print('Starting processing')
    for i, row in tqdm(bed_df.iterrows(), total=len(bed_df.index)):
        assert row['strand'] in {'+', '-'}
        fwd = -1 if row['strand'] == '+' else 1
        for j, mark in enumerate(['ref_m6a', 'ref_5mC'], 1):
            meth_pos = np.array([int(x) for x in row[mark].strip(',').split(',') if x not in incorrect_positions]) - row['start']
            meth_data = np.full_like(meth_pos, j) * fwd
            coo_data.append(meth_data)
            coo_row.append(np.full_like(meth_data, i))
            coo_col.append(meth_pos + row['start_index'])

    print('Finished processing. Creating sparse matrix.')
    coo_data = np.concatenate(coo_data).astype(np.int8)
    coo_row = np.concatenate(coo_row)
    coo_col = np.concatenate(coo_col)
    return coo_matrix((coo_data, (coo_row, coo_col)), shape=(len(bed_df), chromsizes_df['size'].sum()))




if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts output of fibertools-rs extract command to sparse matrix')
    parser.add_argument('bed', help='Path to bed file')
    parser.add_argument('output', help='Path to output npz file')
    parser.add_argument('--chromsizes', help='Path to chromsizes file', default=None)
    parser.add_argument('--chrom', help='Path to genome fasta file', default=None)
    args = parser.parse_args()

    chromsizes_df = process_chrom_sizes(args.chromsizes)
    bed_df = pd.read_table(args.bed).rename(columns=dict(zip(['#ct','st','en'], ['chrom', 'start', 'end'])))
    if args.chrom is not None:
        bed_df.query(f'chrom == "{args.chrom}"', inplace=True)

    coo = create_coo_from_bed(bed_df, args.chromsizes, chromsizes_df)
    print('Converting to csc')
    save_to_h5py(coo.tocsc(), 
        args.output, 
        chrom_start_index=chromsizes_df['cumsum_size'].to_numpy(),
        chrom_name=chromsizes_df.index.to_numpy().astype('S'),
        chrom_size=chromsizes_df['size'].to_numpy().astype(int)
    )
    #save_npz(args.output, coo)