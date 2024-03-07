import h5py


def save_to_h5py(matrix, path, **kwargs):
    with h5py.File(path, 'w') as f:
        f.create_dataset('data', data=matrix.data, compression='gzip')
        f.create_dataset('indices', data=matrix.indices, compression='gzip')
        f.create_dataset('indptr', data=matrix.indptr, compression='gzip')
        for kwarg in kwargs:
            f.create_dataset(kwarg, data=kwargs[kwarg], compression='gzip')
