import h5py

import numpy as np


def read_aai_fail_to_rep_h5(path: str):
    with h5py.File(path, 'r') as f:
        fail_keys = tuple(x.decode() for x in f['fail_keys'])
        rep_keys = tuple(x.decode() for x in f['rep_keys'])
        n_aminos = np.array(f['n_aminos'])
        n_correct = np.array(f['n_correct'])
    return fail_keys, rep_keys, n_correct, n_aminos


if __name__ == '__main__':
    read_aai_fail_to_rep_h5('/srv/home/uqamussi/projects/gunc-chimeras/output/marker/aai_fail_to_rep/TIGR00037.h5')
