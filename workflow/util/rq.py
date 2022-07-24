import time
from collections import Counter
from datetime import timezone
from typing import Callable
from typing import List, Tuple

from redis import Redis
from rq import Queue
from tqdm import tqdm

from workflow.config import REDIS_HOST, REDIS_PASS, RQ_SLEEP_DELAY
from workflow.util.collection import iter_batches
from workflow.util.log import log


def print_job_status(job, delay):
    # Determine the initial progress bar state
    all_statuses = list()
    for dep_job in job.fetch_dependencies():
        cur_status = dep_job.get_status(refresh=False)
        all_statuses.append(cur_status)

    dep_statuses = Counter(all_statuses)
    n_jobs = sum(dep_statuses.values())
    total_done = dep_statuses.get('finished', 0)
    with tqdm(initial=total_done, total=n_jobs, smoothing=0.1, unit='job') as p_bar:

        # Set the initial conditions
        p_bar.start_t = job.created_at.astimezone(timezone.utc).timestamp()

        while job.is_deferred:
            n_done = sum([1 for x in job.fetch_dependencies() if x.get_status(refresh=False) == 'finished'])
            done_delta = n_done - total_done
            p_bar.update(done_delta)
            total_done = n_done
            time.sleep(delay)
    return


def print_queue_status(q: Queue, delay: int):
    # Determine the initial progress bar state
    n_finished = len(q.finished_job_registry)
    n_failed = len(q.failed_job_registry)
    n_deferred = len(q.deferred_job_registry)
    n_started = len(q.started_job_registry)
    n_queued = len(q)
    total_done = n_finished + n_failed

    n_jobs = n_finished + n_failed + n_deferred + n_started + n_queued
    with tqdm(initial=total_done, total=n_jobs, smoothing=0.1, unit='job') as p_bar:
        # Set the initial conditions
        # p_bar.start_t = q.fetch_job(q.job_ids[0]).created_at.astimezone(timezone.utc).timestamp()

        while len(q) > 0:
            n_done = len(q.finished_job_registry) + len(q.failed_job_registry)
            done_delta = n_done - total_done
            p_bar.update(done_delta)
            total_done = n_done
            time.sleep(delay)
    return


def rq_and_wait(job_id: str, fn: Callable, q_args: List[Tuple], queue_name: str, batch_size: int = 1):
    # Prepare the queue
    prep_queue = list()

    # No batching, just take the args
    args_to_enqueue = list()
    if batch_size <= 1:
        args_to_enqueue = q_args
        log(f'No batching, enqueueing {len(args_to_enqueue):,} jobs')
    else:
        for batch in iter_batches(q_args, batch_size):
            args_to_enqueue.append(tuple([tuple(batch, )], ))
        log(f'Batching enabled, enqueuing {len(q_args):,} jobs in batches of {batch_size:,} = {len(args_to_enqueue):,}')

    # Connect to the queue and monitor it
    with Redis(host=REDIS_HOST, password=REDIS_PASS) as conn:
        q = Queue(queue_name, connection=conn, default_timeout='60d')

        # Queue is empty, populate it
        if q.count == 0:
            for arg in args_to_enqueue:
                prep_queue.append(Queue.prepare_data(fn, timeout='60d', ttl='60d',
                                                     failure_ttl='60d', result_ttl='60d',
                                                     args=arg))
            enqueued_jobs = q.enqueue_many(prep_queue)
            log(f'Enqueued {len(enqueued_jobs):,} jobs')

        # Wait for the job to finish
        print(f'Waiting for {RQ_SLEEP_DELAY:,} seconds until jobs are finished...')
        print_queue_status(q, RQ_SLEEP_DELAY)

        # Check if the job didn't succeed
        if not (q.count == 0 and len(q.failed_job_registry) == 0):
            raise Exception(f'Job {job_id} failed!')

    return


def submit_jobs_to_rq(fn: Callable, q_args: List[Tuple], queue_name: str, batch_size: int = 1):
    # Prepare the queue
    prep_queue = list()

    # No batching, just take the args
    args_to_enqueue = list()
    if batch_size <= 1:
        args_to_enqueue = q_args
        log(f'No batching, enqueueing {len(args_to_enqueue):,} jobs')
    else:
        for batch in iter_batches(q_args, batch_size):
            args_to_enqueue.append(tuple([tuple(batch, )], ))
        log(f'Batching enabled, enqueuing {len(q_args):,} jobs in batches of {batch_size:,} = {len(args_to_enqueue):,}')

    # Connect to the queue and monitor it
    with Redis(host=REDIS_HOST, password=REDIS_PASS) as conn:
        q = Queue(queue_name, connection=conn, default_timeout='60d')

        # Queue is empty, populate it
        if q.count == 0:
            for arg in args_to_enqueue:
                prep_queue.append(Queue.prepare_data(fn, timeout='60d', ttl='60d',
                                                     failure_ttl='60d', result_ttl='60d',
                                                     args=arg))
            enqueued_jobs = q.enqueue_many(prep_queue)
            log(f'Enqueued {len(enqueued_jobs):,} jobs')
    return


def rq_wait_for_queue_empty(queue_name: str):
    # Connect to the queue and monitor it
    with Redis(host=REDIS_HOST, password=REDIS_PASS) as conn:
        q = Queue(queue_name, connection=conn, default_timeout='60d')

        # Wait for the job to finish
        print(f'Waiting for {RQ_SLEEP_DELAY:,} seconds until jobs are finished...')
        print_queue_status(q, RQ_SLEEP_DELAY)

        # Check if the job didn't succeed
        if not (q.count == 0 and len(q.failed_job_registry) == 0):
            raise Exception(f'Queue {queue_name} failed!')

    return
