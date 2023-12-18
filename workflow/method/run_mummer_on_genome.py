import os
import subprocess

from workflow.util.paths import get_gid_root

COORD_HEADER = ('qry_gid', 's1', 'e1', 's2', 'e2', 'len_1', 'len_2',
                'pct_identity', 'cov_r', 'cov_q', 'contig_r', 'contig_q')


def parse_output(output: str, qry_gid: str):
    lines = output.splitlines()

    rows = list()
    for line in lines[4:]:
        cols = line.split('\t')

        row = [qry_gid]
        row.extend(map(int, cols[0:6]))
        row.extend(map(float, cols[6:9]))
        row.extend(cols[9:])

        rows.append(row)
    return rows


def run_nucmer(qry_gid: str, ref_gid: str, tmp_dir: str, cache_dir: str):
    path_ref = os.path.join(get_gid_root(ref_gid, cache_dir), f'{ref_gid}.fna')
    path_qry = os.path.join(get_gid_root(qry_gid, cache_dir), f'{qry_gid}.fna')

    # Create the delta file
    path_delta = os.path.join(tmp_dir, f'{ref_gid}_{qry_gid}.delta')
    cmd = [
        'nucmer',
        '--threads=1',
        '--delta',
        path_delta,
        path_ref,
        path_qry,
    ]
    subprocess.check_output(cmd)

    # Create the coordinates file
    cmd = [
        'show-coords',
        '-T',
        '-c',
        '-L', '100',
        path_delta,
    ]
    output = subprocess.check_output(cmd, encoding='utf-8')
    rows = parse_output(output, qry_gid)
    os.remove(path_delta)
    return rows
