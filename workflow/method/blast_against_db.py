
import subprocess

from workflow.model.genome import Genome


def blast_fna_against_db(path_db: str, path_fna: str):

    cmd = [
        'blastn',
        '-db', path_db,
        '-query', path_fna,
        '-outfmt', '0',
        '-out', 'stdout',
        '-num_threads', '4',
    ]
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = proc.communicate()

    if proc.returncode != 0:
        raise Exception('?')


    return


def _test():
    from workflow.external.gtdb_r207_blastdb import GtdbR207BlastDb
    import os
    import tempfile

    path_db = GtdbR207BlastDb().output().path

    genome = Genome('GCA_001577135.1')

    with tempfile.TemporaryDirectory() as tmp_dir:
        path_fna = os.path.join(tmp_dir, 'genome.fna')
        with open(path_fna, 'w') as f:
            for contig, seq in genome.d_fna.items():
                f.write(f'>{contig}\n{seq}\n')

        blast_fna_against_db(path_db, path_fna)

    return


if __name__ == '__main__':
    _test()
