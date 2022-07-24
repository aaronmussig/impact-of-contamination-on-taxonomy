from typing import Dict

from Bio import SeqIO
from frozendict import frozendict


def read_fasta(path) -> Dict[str, str]:
    """Read and make immutable."""
    out = {}
    with open(path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            out[record.id] = str(record.seq)
    return frozendict(out)


def write_fasta(data, path):
    """Write to file."""
    with open(path, "w") as f:
        for gid, seq in data.items():
            f.write(f">{gid}\n{seq}\n")
