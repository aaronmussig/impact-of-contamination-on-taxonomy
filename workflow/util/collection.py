from typing import Callable


def iter_batches(iterable, n=1):
    """Partition a collection into batches of size n."""
    length = len(iterable)
    for ndx in range(0, length, n):
        yield iterable[ndx:min(ndx + n, length)]


def invert_dict(d, obj: Callable) -> dict:
    out = dict()

    if obj is set or obj is frozenset:
        for k, v in d.items():
            if v not in d:
                out[v] = set()
            out[v].add(k)

    elif obj is list or obj is tuple:
        for k, v in d.items():
            if v not in d:
                out[v] = list()
            out[v].append(k)

    if obj is frozenset:
        return {k: frozenset(v) for k, v in out.items()}
    elif obj is tuple:
        return {k: tuple(v) for k, v in out.items()}
    return out
