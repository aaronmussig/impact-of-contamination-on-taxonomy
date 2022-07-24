import os

from workflow.bootstrap.remove_random_contigs_for_bootstrap_rep import RemoveRandomContigsForBootstrapRep
from workflow.bootstrap.select_genomes_for_bootstrapping import SelectGenomesForBootstrapping
from workflow.config import DEBUG
from workflow.config import DIR_OUT_BOOTSTRAP
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class RemoveRandomContigsForBootstrapRepAggregate(LuigiTask):

    def requires(self):
        return {
            'bootstrap_gids': SelectGenomesForBootstrapping(),
            '_removed_contigs': RemoveRandomContigsForBootstrapRep(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_BOOTSTRAP, 'bootstrap_contigs_removed.h5'))

    def run(self):
        log('Bootstrapping - randomly removing contigs (aggregating', title=True)
        self.make_output_dirs()

        log('Loading bootstrap genomes')
        df_bootstrap = self.input()['bootstrap_gids'].read() if not DEBUG else self.input()['bootstrap_gids'].read_cached()

        return
