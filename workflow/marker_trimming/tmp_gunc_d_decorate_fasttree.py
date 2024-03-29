import os
import subprocess
import sys
import tempfile

import dendropy
import luigi
from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DIR_OUT_MARKER_TRIMMING_GUNC
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.marker_trimming.gunc_c_root_fasttree import MarkerTrimmingGuncRootFastTree
from workflow.model.luigi import LuigiTask, LocalTargetTree
from workflow.util.log import log

sys.setrecursionlimit(100000000)


class MarkerTrimmingGuncDecorateFastTreeTmp(LuigiTask):
    target_pct = luigi.FloatParameter()

    @property
    def root_dir(self):
        return DIR_OUT_MARKER_TRIMMING_GUNC

    def requires(self):
        return {
            'tree': MarkerTrimmingGuncRootFastTree(replicate_id=0, target_pct=self.target_pct),
            'meta': GtdbMetadataR207()
        }

    def output(self):
        return LocalTargetTree(
            os.path.join(self.root_dir, f'p{self.target_pct}', 'tmp_phylorank', 'fasttree_undecorated_rooted_decorated.tree'))

    def run(self):
        log(f'MarkerTrimmingGuncDecorateFastTreeTmp()', title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        d_gid_to_taxonomy = df_meta['gtdb_taxonomy'].to_dict()

        output_dir = os.path.dirname(self.output().path)

        log('Reading genomes from tree')
        tree = dendropy.Tree.get_from_path(self.input()['tree'].path, schema='newick', preserve_underscores=True)

        log('Creating taxonomy mapping file for PhyloRank')
        genome_ids = {x.label for x in tree.taxon_namespace}
        log(f'Found {len(genome_ids):,} genomes in tree')

        with tempfile.TemporaryDirectory() as tmp_dir:
            log(f'Working in temporary directory: {tmp_dir}')
            path_tax_file_tmp = os.path.join(tmp_dir, 'taxonomy.tsv')
            path_decorated_tmp = os.path.join(tmp_dir, 'decorated.tree')

            with open(path_tax_file_tmp, 'w') as f:
                for gid in tqdm(genome_ids):
                    if gid.startswith("TEST_"):
                        tax = d_gid_to_taxonomy[gid[5:]]
                    else:
                        tax = d_gid_to_taxonomy[gid]
                    f.write(f'{gid}\t{tax}\n')

            cmd = ['phylorank', 'decorate', '--skip_rd_refine', self.input()['tree'].path, path_tax_file_tmp,
                   path_decorated_tmp]
            log(' '.join(cmd))
            proc = subprocess.Popen(cmd, encoding='utf-8', cwd=tmp_dir, stderr=subprocess.PIPE, stdout=subprocess.PIPE)
            stdout, stderr = proc.communicate()

            log('Copying output files')

            if proc.returncode != 0:
                raise Exception(f'Error running phylorank:\n\n{stderr}\n\n{stdout}\n')

            log(f'Copying result to: {self.output().path}')
            copy_file(path_decorated_tmp, self.output().path, checksum=True)
            copy_file(os.path.join(tmp_dir, 'decorated.tree-summary'),
                      os.path.join(output_dir, 'decorated.tree-summary'), checksum=True)
            copy_file(os.path.join(tmp_dir, 'decorated.tree-table'), os.path.join(output_dir, 'decorated.tree-table'),
                      checksum=True)
            copy_file(os.path.join(tmp_dir, 'decorated.tree-taxonomy'),
                      os.path.join(output_dir, 'decorated.tree-taxonomy'), checksum=True)
            copy_file(os.path.join(tmp_dir, 'phylorank.log'), os.path.join(output_dir, 'phylorank.log'), checksum=True)

        return
