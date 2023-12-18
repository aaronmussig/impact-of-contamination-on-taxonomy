import os
import tempfile
from typing import List, Union
from typing import Optional, Dict

import dendropy
import luigi
import numpy as np
import pandas as pd
from magna.util.disk import get_file_size_fmt, move_file
from magna.util.pandas import optimise_df

from workflow.config import DEBUG
from workflow.util.fasta import read_fasta
from workflow.util.log import log
from workflow.util.obj import get_class_fqn
from workflow.util.paths import maybe_cache


class LuigiTask(luigi.Task):
    """Extended class for Luigi tasks."""

    @property
    def fqn(self):
        return get_class_fqn(self)

    def save_hdf(self, df: pd.DataFrame, index: Optional[Union[List[str], str]] = None,
                 data_columns: Optional[List[str]] = None, path: Optional[str] = None):
        """Save a Pandas DataFrame to a HDF5 file.

        Args:
            df: DataFrame to save.
            index: List of column names to use as index.
            data_columns: Will use index if None.
            path: The path to save to, otherwise self.output().path
        """
        log(f'Saving dataframe {df.shape}')
        if path is None:
            path = self.output().path

        # Set the index and sort it
        if index is not None:
            df.set_index(index, inplace=True)
            if not df.index.is_unique:
                raise Exception('Duplicate index!')
            df.sort_index(inplace=True)

            # If data_columns hasn't been set, then use the index
            if data_columns is None:
                if isinstance(index, str):
                    data_columns = [index]
                else:
                    data_columns = index

        # Minimise dataframe size
        optimise_df(df)

        # Writing to h5 file
        with tempfile.TemporaryDirectory() as tmpdir:
            tmp_path = os.path.join(tmpdir, 'out.h5')
            df.to_hdf(tmp_path, key='root', format='table', data_columns=data_columns,
                      complevel=9, complib='blosc:lz4hc')

            log(f'Copying {get_file_size_fmt(tmp_path)} to {path}')
            move_file(tmp_path, path, checksum=True)

    def make_output_dirs(self):
        if isinstance(self.output(), dict):
            output_dirs = frozenset(os.path.dirname(x.path) for x in self.output().values())
        else:
            output_dirs = frozenset([os.path.dirname(self.output().path)])
        [os.makedirs(x, exist_ok=True) for x in output_dirs]

    def write_sentinel(self):
        with open(self.output().path, 'w') as f:
            f.write(f'{log("Task completed")}\n')


class BaseLocalTarget(luigi.LocalTarget):

    def path_cached(self) -> str:
        return maybe_cache(self.path)


class LocalTargetTree(BaseLocalTarget):

    @staticmethod
    def read_path(path: str) -> dendropy.Tree:
        return dendropy.Tree.get(path=path, schema='newick', preserve_underscores=True)

    def read(self) -> dendropy.Tree:
        return self.read_path(self.path)


class LocalTargetFasta(BaseLocalTarget):

    @staticmethod
    def _read_path(path: str) -> Dict[str, str]:
        return {k: str(v) for k, v in read_fasta(path).items()}

    def read(self) -> Dict[str, str]:
        return self._read_path(self.path)

    def read_cached(self) -> Dict[str, str]:
        path = maybe_cache(self.path)
        return self._read_path(path)
    def maybe_read_cached(self) -> Dict[str, str]:
        if not DEBUG:
            return self.read()
        else:
            return self.read_cached()


class LocalTargetMask(BaseLocalTarget):

    @staticmethod
    def read_path(path: str) -> np.ndarray:
        arr = list()
        with open(path) as f:
            for data in f.read():
                data = data.strip()
                for char in data:
                    if char == '1':
                        arr.append(True)
                    elif char == '0':
                        arr.append(False)
                    else:
                        raise Exception(f'Invalid character {char} in mask file {path}')
        return np.array(arr)

    def read(self) -> np.ndarray:
        return self.read_path(self.path)

    def read_cached(self) -> np.ndarray:
        path = maybe_cache(self.path)
        return self.read_path(path)


class LocalTargetHdf5(BaseLocalTarget):

    def read(self, key: Optional[str] = 'root',
             where: Optional[str] = None,
             start: Optional[int] = None,
             stop: Optional[int] = None,
             columns: Optional[list] = None,
             iterator: Optional[bool] = None,
             chunksize: Optional[int] = None) -> pd.DataFrame:
        return pd.read_hdf(self.path, key=key, where=where, start=start,
                           stop=stop, columns=columns, iterator=iterator,
                           chunksize=chunksize)

    def read_cached(self, key: Optional[str] = 'root',
                    where: Optional[str] = None,
                    start: Optional[int] = None,
                    stop: Optional[int] = None,
                    columns: Optional[list] = None,
                    iterator: Optional[bool] = None,
                    chunksize: Optional[int] = None) -> pd.DataFrame:
        path = maybe_cache(self.path)
        return pd.read_hdf(path, key=key, where=where, start=start,
                           stop=stop, columns=columns, iterator=iterator,
                           chunksize=chunksize)

    def maybe_read_cached(self, key: Optional[str] = 'root',
                          where: Optional[str] = None,
                          start: Optional[int] = None,
                          stop: Optional[int] = None,
                          columns: Optional[list] = None,
                          iterator: Optional[bool] = None,
                          chunksize: Optional[int] = None) -> pd.DataFrame:
        if not DEBUG:
            return self.read(key=key, where=where, start=start,
                             stop=stop, columns=columns, iterator=iterator,
                             chunksize=chunksize)
        else:
            return self.read_cached(key=key, where=where, start=start,
                                    stop=stop, columns=columns, iterator=iterator,
                                    chunksize=chunksize)


class LocalTargetTsv(BaseLocalTarget):

    @staticmethod
    def read_path(path: str) -> pd.DataFrame:
        return pd.read_csv(path, sep='\t')

    def read(self) -> pd.DataFrame:
        return self.read_path(self.path)

    def read_cached(self) -> pd.DataFrame:
        path = maybe_cache(self.path)
        return self.read_path(path)

    def maybe_read_cached(self) -> pd.DataFrame:
        if not DEBUG:
            return self.read()
        else:
            return self.read_cached()
