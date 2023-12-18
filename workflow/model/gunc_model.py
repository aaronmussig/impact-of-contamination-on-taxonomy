

from enum import Enum


class GuncRefDb(Enum):
    GTDB = 'gtdb'
    PRO = 'progenomes'

GUNC_RANKS = ('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
