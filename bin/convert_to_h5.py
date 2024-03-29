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


class FiberSeqExtractor(base_extractor):
    def __init__(self, filename, full_overlap=False, **kwargs):
        """
        Initialize the FiberSeqExtractor.

        Parameters:
            filename (str): Path to the h5 file with methylation data.
            full_overlap (bool): If True, the fibers are filtered to have full overlap with the interval.
        
        """
        super(FiberSeqExtractor, self).__init__(filename, **kwargs)
        self.h5 = h5py.File(filename, 'r')
        self.chrom_indices = self.h5['chrom_names'][:]
        self.full_overlap = full_overlap
    
    def __getitem__(self, interval):
        """
        Get the methylation data for a given GenomicInterval.

        Parameters:
            interval (GenomicInterval): Genomic interval (chr, start, end)
        
        Returns:
            csc_matrix: Sparse matrix with the methylation data.
        """
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
        valid_fibers = self.filter_empty_fibers(
            (new_indices.min(), new_indices.max() + 1),
            col_indices_range
        )
    
        new_indices -= new_indices.min()
        new_shape = (np.max(new_indices) + 1, len(new_indptr) - 1)
        matrix = csc_matrix((new_data, new_indices, new_indptr), shape=new_shape)
        return matrix.tocsr()[valid_fibers]

    def filter_empty_fibers(self, indices, col_indices_range):
        """
        Remove empty fibers from the indices.

        Parameters:
            indices (tuple): Start and end index of the extracted fibers
        """
        fiber_ends = self.h5['fiber_ends'][indices[0]: indices[1]]

        mask_unique = fiber_ends >= col_indices_range[0]
        
        if self.full_overlap:
            fiber_starts = self.h5['fiber_starts'][indices[0]: indices[1]]
            mask_unique *= (fiber_starts <= col_indices_range[0]) & (fiber_ends >= col_indices_range[1])
        return mask_unique

    def interval_to_ids(self, interval):
        """
        Convert genomic interval to index in the genome.

        Parameters:
            interval (GenomicInterval): Genomic interval (chr, start, end)
        
        Returns:
            tuple: Start and end index of the interval in the genome.
        """
        chrom = interval.chrom.encode('utf-8')
        chr_index = np.where(self.chrom_indices == chrom)[0]
        if len(chr_index) == 0:
            raise ValueError(f'Chromosome {chrom} not found in chroms used to create the file.\n', 'Chromlist: ', self.chrom_indices)
        chrom_start_index = self.h5['chrom_start_indices'][chr_index[0]]
        assert self.h5['chrom_sizes'][chr_index[0]] >= interval.end, f"Chromosome size: {self.h5['chrom_sizes'][chr_index[0]]} < {interval.end}"
        start_index = chrom_start_index + interval.start
        end_index = chrom_start_index + interval.end
        return start_index, end_index


from genome_tools.genomic_interval import genomic_interval as GenomicInterval
from scipy.sparse import coo_matrix
from tqdm import tqdm
import argparse
import pandas as pd
import numpy as np


def process_chrom_sizes(chromsizes_path):
    """
    Read chromsizes file and create a dataframe with genomic coordinate conversion.

    Parameters:
        chromsizes_path (str): Path to the chromsizes file. No header is expected.

    Returns:
        pd.DataFrame: Dataframe with the following columns: size, cumsum_size. Index is the chromosome name.
    """
    df = pd.read_table(chromsizes_path, names=['chrom', 'size']).set_index('chrom').sort_index()
    df['cumsum_size'] = df['size'].cumsum().shift(1).fillna(0)
    return df


def create_coo_from_bed(bed_df, genome_length):
    """
    Create a h5 file from a bed file.

    Parameters:
        bed_df (pd.DataFrame): Bed file with the following columns: chrom, start, end, ref_m6a, strand.
        genome_length (int): Length of the genome
    
    Returns:
        coo_matrix: Sparse matrix with the methylation data.
    """
    incorrect_positions = set(["", "-1", "."])
    
    coo_data = []
    coo_row = []
    coo_col = []
    print('Starting processing')
    pbar = tqdm(bed_df.iterrows(), total=len(bed_df.index))
    for i, row in pbar:
        assert row['strand'] in {'+', '-'}
        fwd = 1 if row['strand'] == '+' else -1
        for j, mo in enumerate(['ref_m6a'], 1):
            mo_pos = np.array([int(x) for x in row[mo].strip(',').split(',') if x not in incorrect_positions]) - row['start']
            if (mo_pos.size != np.unique(mo_pos).size):
                print(row['fiber'])
                unique, counts = np.unique(mo_pos, return_counts=True)
                repetitive_numbers = unique[counts > 1]
                print(repetitive_numbers)
                raise AssertionError
            coo_data.append(np.full_like(mo_pos, j) * fwd)
            coo_row.append(np.full_like(mo_pos, i))
            coo_col.append(mo_pos + row['start_index'])

    pbar.close()
    print('Finished processing. Creating sparse matrix.')
    sys.stdout.flush()
    coo_data = np.concatenate(coo_data).astype(np.int8)
    coo_row = np.concatenate(coo_row)
    coo_col = np.concatenate(coo_col)
    return coo_matrix((coo_data, (coo_row, coo_col)), shape=(len(bed_df.index), genome_length))


def save_to_h5(matrix, path, **kwargs):
    """
    Save a sparse matrix and additional data to a h5 file.

    Parameters:
        matrix (csc_matrix): Sparse matrix to save.
        path (str): Path to the output file.
        kwargs (dict): Additional data to save.
    """
    with h5py.File(path, 'w') as f:
        f.create_dataset('data', data=matrix.data, compression='gzip')
        f.create_dataset('indices', data=matrix.indices, compression='gzip')
        f.create_dataset('indptr', data=matrix.indptr, compression='gzip')
        for kwarg in kwargs:
            f.create_dataset(kwarg, data=kwargs[kwarg], compression='gzip')


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Converts output of fibertools-rs extract command to sparse matrix')
    parser.add_argument('bed', help='Path to bed file')
    parser.add_argument('output', help='Path to output npz file')
    parser.add_argument('--chromsizes', help='Path to chromsizes file', default=None)
    parser.add_argument('--chrom', help='Chromosome to process', default=None)
    args = parser.parse_args()

    chromsizes_df = process_chrom_sizes(args.chromsizes)
    bed_df = pd.read_table(args.bed).rename(columns=dict(zip(['#ct','st','en'], ['chrom', 'start', 'end'])))
    if args.chrom is not None:
        bed_df = bed_df.query(f'chrom == "{args.chrom}"').reset_index(drop=True)

    bed_df['chrom_index'] = bed_df['chrom'].map(chromsizes_df['cumsum_size'])
    bed_df['start_index'] = bed_df['chrom_index'] + bed_df['start']
    bed_df['end_index'] = bed_df['chrom_index'] + bed_df['end']
    csc = create_coo_from_bed(bed_df, chromsizes_df['size'].sum()).tocsc()

    print('Saving to h5py...')
    save_to_h5(csc,
        args.output,
        chrom_start_indices=chromsizes_df['cumsum_size'].to_numpy().astype(int),
        chrom_names=chromsizes_df.index.to_numpy().astype('S'),
        chrom_sizes=chromsizes_df['size'].to_numpy().astype(int),
        fiber_starts=bed_df['start_index'].to_numpy(),
        fiber_ends=bed_df['end_index'].to_numpy(),
    )