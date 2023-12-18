import os
from collections import defaultdict

import pandas as pd
from luigi import LocalTarget

from workflow.config import DEBUG, DIR_OUT_MARKER, PCT_VALUES
from workflow.config import R207_AR53_HMM, R207_BAC120_HMM
from workflow.external.gtdb_metadata import GtdbMetadataR207
from workflow.marker.get_markers_on_contigs_for_failed import GetMarkersOnContigsForFailed
from workflow.model.luigi import LuigiTask, LocalTargetHdf5
from workflow.util.log import log


class PctRemovedVsPctOfDomainMarkers(LuigiTask):

    def requires(self):
        return {
            'meta': GtdbMetadataR207(),
            'markers': GetMarkersOnContigsForFailed(),
        }

    def output(self):
        return LocalTargetHdf5(os.path.join(DIR_OUT_MARKER, 'pct_removed_vs_pct_of_domain_markers.h5'))

    def get_gid_to_domain(self):
        df = self.input()['meta'].read_cached() if DEBUG else self.input()['meta'].read()
        return df['domain'].to_dict()

    def run(self):
        log('Calculating the percentage of genome removed vs the pct of domain markers present', title=True)
        self.make_output_dirs()

        log('Loading gid to domain mapping')
        d_gid_to_domain = self.get_gid_to_domain()

        log('Loading marker dataframe')
        df_markers = self.input()['markers'].read_cached() if DEBUG else self.input()['markers'].read()

        log('Parsing marker data')
        d_gid_to_pct_data = parse_markers(df_markers, d_gid_to_domain)
        del df_markers

        log('Converting into rows')
        rows = convert_to_rows(d_gid_to_pct_data)
        df = pd.DataFrame(rows)
        df.sort_values(by=['gid', 'pct'], inplace=True)

        if not DEBUG:
            self.save_hdf(df)

        return


def convert_to_rows(d_gid_to_pct_data):
    out = list()
    for gid, d_pct_to_pct_domain_markers in d_gid_to_pct_data.items():
        base_pct = d_pct_to_pct_domain_markers[0]
        out.append({
            'gid': gid,
            'pct': 0,
            'pct_markers': base_pct,
        })

        for pct in PCT_VALUES:
            # bring forward the values
            if pct not in d_pct_to_pct_domain_markers:
                previous_data = out[-1]
                out.append({
                    'gid': gid,
                    'pct': pct,
                    'pct_markers': previous_data['pct_markers']
                })
            else:
                out.append({
                    'gid': gid,
                    'pct': pct,
                    'pct_markers': d_pct_to_pct_domain_markers[pct]
                })

    return out


def parse_markers(df_markers, d_gid_to_domain):
    d_gid_to_pct_data = defaultdict(dict)

    for _, row in df_markers.iterrows():
        domain = d_gid_to_domain[row['gid']]
        if domain == 'd__Archaea':
            marker_set = frozenset(R207_AR53_HMM)
        elif domain == 'd__Bacteria':
            marker_set = frozenset(R207_BAC120_HMM)
        else:
            raise Exception(f'Unknown domain: {domain}')

        unq_markers_present = set()
        for marker in marker_set:
            if row[marker] in {'unq', 'muq'}:
                unq_markers_present.add(marker)

        # Domain markers
        domain_markers = marker_set.intersection(unq_markers_present)
        pct_domain_markers = 100 * len(domain_markers) / len(marker_set)

        d_gid_to_pct_data[row['gid']][row['pct']] = pct_domain_markers

        if DEBUG and len(d_gid_to_pct_data) > 10:
            break

    return d_gid_to_pct_data
