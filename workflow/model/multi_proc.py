import multiprocessing as mp

from tqdm import tqdm

from workflow.util.collection import iter_batches


def run_mp_batched(fn, cpus, batch_size, queue):
    results = list()
    with tqdm(total=len(queue), smoothing=0.1) as pbar:
        for batch in iter_batches(queue, batch_size):
            with mp.Pool(processes=cpus) as pool:
                for result in pool.imap_unordered(fn, batch):
                    results.append(result)
                    pbar.update()
    return


def _test_fn(job):
    return job ** 2


def _test():
    jobs = list(range(100))
    out = run_mp_batched(_test_fn, 1, 10, jobs)

    print()


if __name__ == '__main__':
    _test()
