import os
import tempfile

from magna.util.disk import copy_file
from tqdm import tqdm

from workflow.config import DEBUG, DIR_OUT_FASTTREE_FULL_TREE_NON_REPS
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.external.gtdb_msa import GtdbMsaBacNonRepsR207
from workflow.gunc_helper.aggregate_max_css_level_merged import AggregateMaxCssLevelMerged
from workflow.model.luigi import LuigiTask, LocalTargetFasta
from workflow.util.log import log


class FastTreeFullSetNonRepsCreateBatchesStandard(LuigiTask):

    def requires(self):
        return {
            'max_css': AggregateMaxCssLevelMerged(),
            'meta': GtdbMetadataR207(),
            'msa': GtdbMsaBacNonRepsR207(),
        }

    def output(self):
        return LocalTargetFasta(os.path.join(DIR_OUT_FASTTREE_FULL_TREE_NON_REPS, 'standard', 'msa.faa'))

    def load_true_msa(self, keep_gids):
        d_gid_to_msa = self.input()['msa'].maybe_read_cached()

        out = dict()
        for gid, seq in d_gid_to_msa.items():
            if gid[3:] in keep_gids:
                out[gid[3:]] = seq
        if len(out) != len(keep_gids):
            raise Exception(f'Expected {len(keep_gids)} gids, but only found {len(out)}')
        return out

    def run(self):
        log(f'Creating batches for FastTree full tree non reps, bac gids (standard)',
            title=True)
        self.make_output_dirs()

        log('Loading metadata')
        df_meta = self.input()['meta'].maybe_read_cached()
        rep_gids = set(df_meta[(df_meta['gtdb_representative'] == 't') & (df_meta['domain'] == 'd__Bacteria')].index)
        bac_gids = set(df_meta[df_meta['domain'] == 'd__Bacteria'].index)

        log('Loading MaxCSS')
        df_max_css = self.input()['max_css'].maybe_read_cached()
        bac_gids_fail = bac_gids.intersection(set(df_max_css.index))

        log('Loading true MSA')
        d_true_msa = self.load_true_msa(rep_gids.union(bac_gids_fail))

        with tempfile.TemporaryDirectory() as tmp_dir:
            path_msa_tmp = os.path.join(tmp_dir, 'msa.faa')

            log('Writing MSA')
            with open(path_msa_tmp, 'w') as fh:
                for gid, seq in tqdm(sorted(d_true_msa.items())):
                    fh.write(f'>{gid}\n{seq}\n')

                if not DEBUG:
                    log('Copying MSA')
                    copy_file(path_msa_tmp, self.output().path, checksum=True)
        return
