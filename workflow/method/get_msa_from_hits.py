import os
import subprocess
import tempfile

from workflow.config import R207_AR53_HMM, R207_BAC120_HMM, DIR_PFAM_33, DIR_TIGRFAM_15, \
    R207_MARKER_LENGTHS, PATH_R207_AR53_MASK, PATH_R207_BAC120_MASK
from workflow.method.get_genome_domain_from_markers import get_genome_domain_from_markers
from workflow.method.get_marker_hits_for_gid import get_marker_hits_for_gid
from workflow.model.hmm_aligner import get_aligned_marker, get_aligned_markers
from workflow.model.tophit import TopHitPfamFile, TopHitTigrFile
from workflow.util.fasta import read_fasta
from workflow.util.paths import get_gid_root


def align_marker(marker, path_faa):
    if marker.startswith('PF'):
        hmm_path = os.path.join(DIR_PFAM_33, 'individual_hmms', f'{marker}.hmm')
    elif marker.startswith('TIGR'):
        hmm_path = os.path.join(DIR_TIGRFAM_15, 'individual_hmms', f'{marker}.HMM')
    else:
        raise Exception('Unknown marker')

    proc = subprocess.Popen(["hmmalign", "--outformat", "Pfam", hmm_path, path_faa],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
    stdout, stderr = proc.communicate()
    if proc.returncode != 0:
        raise Exception(f'Error running hmmalign: {stderr}')
    return get_aligned_markers(stdout)


def get_msa_from_hits(results_markers, marker_dir=None, user_domain=None):
    domain = get_genome_domain_from_markers(results_markers)
    if domain and user_domain and domain != user_domain:
        print(f'Domain disagreement! {domain} != {user_domain}')
        domain = user_domain
    single_copy_hits = {**results_markers['muq'], **results_markers['unq']}

    if domain == 'd__Archaea':
        marker_ids = R207_AR53_HMM
        path_mask = PATH_R207_AR53_MASK
    elif domain == 'd__Bacteria':
        marker_ids = R207_BAC120_HMM
        path_mask = PATH_R207_BAC120_MASK
    else:
        raise Exception('Unknown domain')

    with open(path_mask) as f:
        mask = [x == '1' for x in f.read().strip()]

    markers_to_align = set(single_copy_hits).intersection(set(marker_ids))

    d_marker_to_aln = dict()
    with tempfile.TemporaryDirectory(prefix='gunc-chim-msa') as tmp_dir:
        for marker in sorted(markers_to_align):
            hit = single_copy_hits[marker]

            path = os.path.join(tmp_dir, f'{marker}.faa')
            with open(path, 'w') as f:
                f.write(f'>{marker}\n{hit["seq"]}\n')

            if marker_dir is None:
                if marker.startswith('PF'):
                    hmm_path = os.path.join(DIR_PFAM_33, 'individual_hmms', f'{marker}.hmm')
                elif marker.startswith('TIGR'):
                    hmm_path = os.path.join(DIR_TIGRFAM_15, 'individual_hmms', f'{marker}.HMM')
                else:
                    raise Exception('Unknown marker')
            else:
                hmm_path = os.path.join(marker_dir, f'{marker}.hmm')

            proc = subprocess.Popen(["hmmalign", "--outformat", "Pfam", hmm_path, path],
                                    stdout=subprocess.PIPE, stderr=subprocess.PIPE, encoding='utf-8')
            stdout, stderr = proc.communicate()
            d_marker_to_aln[marker] = get_aligned_marker(marker, stdout)

    seqs = list()
    for marker in marker_ids:
        if marker not in d_marker_to_aln:
            seqs.append('-' * R207_MARKER_LENGTHS[marker])
        else:
            seqs.append(d_marker_to_aln[marker])

    msa = ''.join(seqs)
    msa_trim = [x for x, y in zip(msa, mask) if y]
    return ''.join(msa_trim)


def main():
    gid = 'GCF_003052605.1'
    gid_root = get_gid_root(gid)

    # Read the fasta file
    d_faa = dict(read_fasta(os.path.join(gid_root, 'prodigal', f'{gid}.faa')))

    pfam_th = TopHitPfamFile(os.path.join(gid_root, 'prodigal', f'{gid}_pfam_tophit.tsv'))
    pfam_th.read()
    tigr_th = TopHitTigrFile(os.path.join(gid_root, 'prodigal', f'{gid}_tigrfam_tophit.tsv'))
    tigr_th.read()

    # omit_contigs = {'NZ_QASO01000049.1'}
    omit_contigs = None
    results_markers = get_marker_hits_for_gid(d_faa, pfam_th, tigr_th, omit_contigs)

    msa = get_msa_from_hits(results_markers)
    return


if __name__ == '__main__':
    main()
