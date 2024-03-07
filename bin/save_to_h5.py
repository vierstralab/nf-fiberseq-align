import h5py


def save_matrix_to_h5py(matrix, path):
    with h5py.File(path, 'w') as f:
        f['data'] = matrix.data
        f['indices'] = matrix.indices
        f['indptr'] = matrix.indptr
