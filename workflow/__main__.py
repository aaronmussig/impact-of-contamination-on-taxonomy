import os
import shutil
import socket
import time

import luigi as lg
import typer
from redis import Redis
from rq import Queue
from rq.command import send_shutdown_command
from rq.worker import Worker

from workflow.config import REDIS_HOST, REDIS_PASS, DIR_CACHE
from workflow.fastani_interspecies.plots import FastAniInterspeciesPlots
from workflow.fastani_interspecies.run_fastani_interspecies_to_genus_ani_for_random_non_reps_agg import \
    FastAniInterspeciesToGenusAniForRandomNonRepsAgg

app = typer.Typer()


class AllTasks(lg.Task):
    """The final node in the Luigi pipeline, creates the DAG."""

    def requires(self):
        return [
            FastAniInterspeciesPlots()
        ]

    def complete(self):
        return False


@app.command()
def run():
    typer.echo('Running workflow...')
    lg.build([AllTasks()], workers=1, local_scheduler=True, log_level="WARNING")
    typer.echo('Done.')


@app.command()
def clear():
    if DIR_CACHE is not None:
        delete = typer.confirm(f'Delete directory? {DIR_CACHE}')
        if delete:
            if os.path.isdir(DIR_CACHE):
                shutil.rmtree(DIR_CACHE)
            typer.echo('Done.')
        else:
            typer.echo('Nothing to do.')
    else:
        typer.echo('Cache is not set, nothing to do.')


@app.command()
def rq(queue_name: str):
    typer.echo(f'Running RQ worker: {queue_name}')
    # os.nice(10)

    # Start a worker with a custom name
    with Redis(host=REDIS_HOST, password=REDIS_PASS) as conn:
        q = Queue(queue_name, connection=conn, default_timeout='60d')
        print(f'Queue size: {len(q)}')
        host_name = socket.gethostname().replace('.ace.uq.edu.au', '')
        worker_name = f'{host_name}-{time.time()}'
        worker = Worker([q], connection=conn, job_monitoring_interval=90, name=worker_name)
        worker.work(burst=True)

    typer.echo('Done.')


@app.command()
def rq_worker_stop(worker_name: str):
    typer.echo(f'Stopping RQ worker: {worker_name}')
    with Redis(host=REDIS_HOST, password=REDIS_PASS) as conn:
        worker = [x for x in Worker.all(conn) if x.name == worker_name]
        if len(worker) == 0:
            typer.echo('No worker found.')
            return
        worker = worker[0]
        send_shutdown_command(conn, worker)
    typer.echo('Done.')


def main():
    app()


if __name__ == '__main__':
    main()
