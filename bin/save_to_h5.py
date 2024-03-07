import h5py


def save_to_h5py(matrix, path, **kwargs):
    with h5py.File(path, 'w') as f:
        f['data'] = matrix.data
        f['indices'] = matrix.indices
        f['indptr'] = matrix.indptr
        for kwarg in kwargs:
            f[kwarg] = kwargs[kwarg]
