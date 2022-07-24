import os

from magna.util.disk import copy_file

import workflow
from workflow.config import R95_ROOT, DIR_OUT_GENOMES, ACE_R95_ROOT, ACE_R207_ROOT, DIR_CACHE, DIR_CODE, NCBI_FTP_ROOT


def get_gid_root(gid: str, root: str = DIR_OUT_GENOMES) -> str:
    """Return the path to the root directory of the genome output folder."""
    return os.path.join(root, gid[0:3], gid[4:7], gid[7:10], gid[10:13])


def get_gid_root_cutoff(gid: str, cutoff: int, root: str = DIR_OUT_GENOMES) -> str:
    """Return the path to the root directory of the genome output folder."""
    return os.path.join(root, gid[0:3], gid[4:7], gid[7:10], gid[10:13], 'cutoff', str(cutoff))

def get_gid_root_pct(gid: str, pct: int, root: str = DIR_OUT_GENOMES) -> str:
    """Return the path to the root directory of the genome output folder."""
    return os.path.join(root, gid[0:3], gid[4:7], gid[7:10], gid[10:13], 'pct', str(pct))


def get_gid_r95_root(gid):
    if gid[0:3] == 'GCA':
        db = 'genbank'
    else:
        db = 'refseq'
    r95_root = os.path.join(R95_ROOT, db, gid[0:3], gid[4:7], gid[7:10], gid[10:13])
    if not os.path.isdir(r95_root):
        raise IOError(f'Cannot find: {gid}')
    top = os.listdir(r95_root)[0]
    return os.path.join(r95_root, top)

def get_gid_r207_root(gid):
    if gid[0:3] == 'GCA':
        db = 'genbank'
    else:
        db = 'refseq'
    r207_root = os.path.join(ACE_R207_ROOT, db, gid[0:3], gid[4:7], gid[7:10], gid[10:13])
    if not os.path.isdir(r207_root):
        raise IOError(f'Cannot find: {gid}')

    top = os.listdir(r207_root)[0]
    return os.path.join(r207_root, top)


def get_ace_gid_root(gid):
    """Required as some UBA genomes are in alternative paths. This is easiest."""
    if gid[0:3] == 'GCA':
        db = 'genbank'
    else:
        db = 'refseq'

    r95_root = os.path.join(ACE_R95_ROOT, db, gid[0:3], gid[4:7], gid[7:10], gid[10:13])
    r202_root = os.path.join(ACE_R207_ROOT, db, gid[0:3], gid[4:7], gid[7:10], gid[10:13])

    use_root = r95_root
    if not os.path.isdir(r95_root):
        use_root = r202_root
        if not os.path.isdir(r202_root):
            raise IOError(f'Cannot find: {gid}')

    top = os.listdir(use_root)[0]
    return os.path.join(use_root, top)


def maybe_cache(source: str) -> str:
    """If the source file is not in the cache, copy it there."""
    if DIR_CACHE is None:
        return source
    new_path = os.path.join(DIR_CACHE, source[1:])
    if not os.path.exists(new_path):
        os.makedirs(os.path.dirname(new_path), exist_ok=True)
        copy_file(source, new_path)
    return new_path


def get_module_remote_path(module) -> str:
    """Return the path to the module."""
    root_path = os.path.dirname(os.path.dirname(workflow.__file__))
    return os.path.join(DIR_CODE, module.__file__[len(root_path) + 1:])
