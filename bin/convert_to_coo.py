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

def process_chrom_sizes(chromsizes_path):
    df = pd.read_table(chromsizes_path, names=['chrom', 'size']).set_index('chrom').sort_index()
    df['cumsum_size'] = df['size'].cumsum().shift(1).fillna(0)
    return df


def create_h5_from_bed(bed_path, chromsizes_path, fasta_path):
    """
    Create a h5 file from a bed file.

    Parameters:
        bed_path (str): Path to the bed file.
        chromsizes_path (str): Path to the chromsizes file.
        h5_path (str): Path to the h5 file.
    """
    fasta_extractor = FastaExtractor(fasta_path)
    chromsizes_df = process_chrom_sizes(chromsizes_path)
    incorrect_positions = set(["", "-1", "."])
    bed_df = pd.read_table(bed_path).rename(columns=dict(zip(['#ct','st','en'], ['chrom', 'start', 'end'])))
    bed_df['chrom_index'] = bed_df['chrom'].map(chromsizes_df['cumsum_size'])
    bed_df['start_index'] = bed_df['chrom_index'] + bed_df['start']
    
    coo_data = []
    coo_row = []
    coo_col = []
    for i, row in bed_df.iterrows():
        data = set(np.array([int(x) for x in row['ref_m6a'].strip(',').split(',') if x not in incorrect_positions]) - row['start'])
        sequence = fasta_extractor[GenomicInterval(row['chrom'], row['start'], row['end'])]
        if row['strand'] == '-':
            idxs = [x.start() for x in re.finditer('T', sequence)]
        else:
            idxs = [x.start() for x in re.finditer('A', sequence)]
#         print(min(idxs), min(data, default=0))
        for idx in idxs:
            coo_data.append(2 if idx in data else 1)
            coo_col.append(row['start_index'] + idx)
            coo_row.append(i)
    return coo_matrix((coo_data, (coo_row, coo_col)), shape=(len(bed_df), chromsizes_df['size'].sum()))
#     csc = coo_matrix((coo_data, (coo_row, coo_col)), shape=(len(bed_df), chromsizes_df['size'].sum())).tocsc()
#     with h5py.File(h5_path, 'w') as h5:
#         h5.create_dataset('indptr', data=csc.indptr)
#         h5.create_dataset('indices', data=csc.indices)
#         h5.create_dataset('data', data=csc.data)




if __name__ == "__main__":
    chromsizes_path = "/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.chrom_sizes"
    fasta_path = '/net/seq/data/genomes/human/GRCh38/noalts/GRCh38_no_alts.fa'
    coo = create_h5_from_bed(sys.argv[1], chromsizes_path, fasta_path)
    save_npz(sys.argv[2], coo)