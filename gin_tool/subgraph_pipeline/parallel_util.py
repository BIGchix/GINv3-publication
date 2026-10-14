"""CPU parallel helpers for context subgraph build (no GPU).

Use ``n_jobs=1`` (default) for the legacy serial path. ``n_jobs>1`` fans out
independent units of work via threads (shared CSR reads) or processes
(GIL-bound Python loops such as per-gene dogma repair).

Process pools use the ``spawn`` start method so they remain safe after earlier
thread-pool steps in the same CLI run.
"""

from __future__ import annotations

import multiprocessing as mp
import os
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor
from typing import Callable, Sequence, TypeVar

T = TypeVar("T")
R = TypeVar("R")


def resolve_n_jobs(n_jobs: int | None) -> int:
    """Normalize ``n_jobs``: ``None``/``0`` → CPU count (capped); ``<0`` → abs."""
    if n_jobs is None or n_jobs == 0:
        cpu = os.cpu_count() or 1
        return max(1, min(cpu, 32))
    if n_jobs < 0:
        cpu = os.cpu_count() or 1
        return max(1, min(cpu, abs(int(n_jobs))))
    return max(1, int(n_jobs))


def map_threads(
    fn: Callable[[T], R],
    items: Sequence[T],
    *,
    n_jobs: int | None = 1,
    chunksize: int = 1,
) -> list[R]:
    """Preserve input order. Serial when ``n_jobs<=1`` or ``len(items)<=1``."""
    n_jobs = resolve_n_jobs(n_jobs)
    if n_jobs <= 1 or len(items) <= 1:
        return [fn(x) for x in items]
    workers = min(n_jobs, len(items))
    with ThreadPoolExecutor(max_workers=workers) as pool:
        return list(pool.map(fn, items, chunksize=max(1, int(chunksize))))


def map_processes(
    fn: Callable[[T], R],
    items: Sequence[T],
    *,
    n_jobs: int | None = 1,
    chunksize: int | None = None,
    initializer: Callable[..., None] | None = None,
    initargs: tuple = (),
) -> list[R]:
    """Preserve input order via ProcessPool (spawn). Serial when ``n_jobs<=1``."""
    n_jobs = resolve_n_jobs(n_jobs)
    if n_jobs <= 1 or len(items) <= 1:
        if initializer is not None:
            initializer(*initargs)
        return [fn(x) for x in items]
    workers = min(n_jobs, len(items))
    cs = chunksize
    if cs is None:
        cs = max(1, len(items) // (workers * 4) or 1)
    ctx = mp.get_context("spawn")
    with ProcessPoolExecutor(
        max_workers=workers,
        mp_context=ctx,
        initializer=initializer,
        initargs=initargs,
    ) as pool:
        return list(pool.map(fn, items, chunksize=cs))


def run_parallel_pair(
    left: Callable[[], R],
    right: Callable[[], R],
    *,
    n_jobs: int | None = 1,
) -> tuple[R, R]:
    """Run two independent callables; uses threads when ``n_jobs>=2``."""
    n_jobs = resolve_n_jobs(n_jobs)
    if n_jobs <= 1:
        return left(), right()
    with ThreadPoolExecutor(max_workers=2) as pool:
        f1 = pool.submit(left)
        f2 = pool.submit(right)
        return f1.result(), f2.result()
