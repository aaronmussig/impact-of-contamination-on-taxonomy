from enum import Enum

from workflow.model.gunc_model import GUNC_RANKS


class TaxDomain(Enum):
    BACTERIA = 'd__Bacteria'
    ARCHAEA = 'd__Archaea'


class Taxonomy:
    __slots__ = ('d', 'p', 'c', 'o', 'f', 'g', 's')

    def __init__(self, d, p=None, c=None, o=None, f=None, g=None, s=None):
        self.d = d
        self.p = p
        self.c = c
        self.o = o
        self.f = f
        self.g = g
        self.s = s

    def __repr__(self):
        return str(self)

    def __str__(self):
        ranks = list()
        if self.d:
            ranks.append(self.d)
        if self.p:
            ranks.append(self.p)
        if self.c:
            ranks.append(self.c)
        if self.o:
            ranks.append(self.o)
        if self.f:
            ranks.append(self.f)
        if self.g:
            ranks.append(self.g)
        if self.s:
            ranks.append(self.s)
        return ';'.join(ranks)

    def from_gunc_rank(self, rank: GUNC_RANKS):
        if rank == 'kingdom':
            return self.d
        elif rank == 'phylum':
            return self.p
        elif rank == 'class':
            return self.c
        elif rank == 'order':
            return self.o
        elif rank == 'family':
            return self.f
        elif rank == 'genus':
            return self.g
        elif rank == 'species':
            return self.s
        else:
            raise ValueError(f'Unknown rank: {rank}')