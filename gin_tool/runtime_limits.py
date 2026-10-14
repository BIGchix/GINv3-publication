"""Process-level resource limits for long-running gin_tool scripts."""

from __future__ import annotations

import os
import resource
import shutil
from dataclasses import dataclass
from pathlib import Path


def system_mem_total_bytes() -> int:
    """Read MemTotal from /proc/meminfo (Linux)."""
    with open("/proc/meminfo", encoding="utf-8") as f:
        for line in f:
            if line.startswith("MemTotal:"):
                return int(line.split()[1]) * 1024
    raise RuntimeError("MemTotal not found in /proc/meminfo")


@dataclass(frozen=True)
class MemoryCapInfo:
    fraction: float
    mem_total_bytes: int
    limit_bytes: int
    soft_bytes: int
    hard_bytes: int
    max_gib_clamp: float | None = None

    def to_log_line(self) -> str:
        gb = self.limit_bytes / (1024**3)
        tot = self.mem_total_bytes / (1024**3)
        clamp = (
            f" max_gib={self.max_gib_clamp:g}" if self.max_gib_clamp is not None else ""
        )
        return (
            f"memory_cap fraction={self.fraction:.2f}{clamp} "
            f"limit={gb:.1f}GiB of MemTotal={tot:.1f}GiB "
            f"(RLIMIT_AS soft={self.soft_bytes} hard={self.hard_bytes})"
        )


def apply_memory_cap(
    fraction: float = 0.5,
    *,
    min_bytes: int = 2 * 1024**3,
    max_bytes: int | None = None,
    max_gib: float | None = None,
) -> MemoryCapInfo:
    """
    Cap this process virtual address space (RLIMIT_AS).

    Limit is ``fraction * MemTotal``, then optionally clamped by ``max_bytes`` /
    ``max_gib``, and never below ``min_bytes``. Soft and hard limits are set equal.
    """
    if not (0.0 < fraction <= 1.0):
        raise ValueError(f"fraction must be in (0, 1], got {fraction}")
    total = system_mem_total_bytes()
    limit = max(int(total * fraction), int(min_bytes))
    clamp_gib = max_gib
    if max_gib is not None:
        max_bytes = int(max_gib * (1024**3)) if max_bytes is None else min(max_bytes, int(max_gib * (1024**3)))
    if max_bytes is not None:
        limit = min(limit, int(max_bytes))
        limit = max(limit, int(min_bytes))
    soft, hard = resource.getrlimit(resource.RLIMIT_AS)
    if hard != resource.RLIM_INFINITY:
        limit = min(limit, hard)
    resource.setrlimit(resource.RLIMIT_AS, (limit, limit))
    soft2, hard2 = resource.getrlimit(resource.RLIMIT_AS)
    return MemoryCapInfo(
        fraction=fraction,
        mem_total_bytes=total,
        limit_bytes=limit,
        soft_bytes=soft2,
        hard_bytes=hard2,
        max_gib_clamp=clamp_gib,
    )


_THREAD_ENV_KEYS = (
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "NUMEXPR_NUM_THREADS",
    "BLIS_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def apply_thread_cap(n_threads: int = 2) -> dict[str, str]:
    """Pin BLAS/OpenMP thread env vars (call before heavy numeric imports when possible)."""
    if n_threads < 1:
        raise ValueError(f"n_threads must be >= 1, got {n_threads}")
    applied: dict[str, str] = {}
    val = str(int(n_threads))
    for key in _THREAD_ENV_KEYS:
        os.environ[key] = val
        applied[key] = val
    return applied


def disk_free_bytes(path: os.PathLike[str] | str) -> int:
    return shutil.disk_usage(Path(path)).free


def require_disk_free(
    path: os.PathLike[str] | str,
    *,
    min_free_gib: float = 50.0,
) -> int:
    """Raise RuntimeError if free space on ``path``'s filesystem is below threshold."""
    free = disk_free_bytes(path)
    need = int(min_free_gib * (1024**3))
    if free < need:
        raise RuntimeError(
            f"disk free {free / (1024**3):.1f}GiB < required {min_free_gib:g}GiB "
            f"on filesystem of {Path(path).resolve()}"
        )
    return free


def find_conflicting_batch_pids(
    *,
    patterns: tuple[str, ...] = (
        "batch_kegg_subgraph_spn_v2.py",
        "batch_kegg_gene_met_spn_v2.py",
    ),
    exclude_pid: int | None = None,
) -> list[tuple[int, str]]:
    """Return (pid, cmdline) for other running kegg SPN batch processes."""
    me = exclude_pid if exclude_pid is not None else os.getpid()
    hits: list[tuple[int, str]] = []
    proc = Path("/proc")
    if not proc.is_dir():
        return hits
    for entry in proc.iterdir():
        if not entry.name.isdigit():
            continue
        pid = int(entry.name)
        if pid == me:
            continue
        cmdline_path = entry / "cmdline"
        try:
            raw = cmdline_path.read_bytes()
        except OSError:
            continue
        cmd = raw.replace(b"\x00", b" ").decode("utf-8", errors="replace").strip()
        if any(p in cmd for p in patterns):
            hits.append((pid, cmd))
    return hits


def safe_relpath(path: os.PathLike[str] | str, root: os.PathLike[str] | str) -> str:
    """Return path relative to root when possible; otherwise absolute string.

    Handles symlink roots (e.g. ``/home/...`` vs ``/disk1/...``) and relative
    inputs that would make bare ``Path.relative_to`` raise ``ValueError``.
    """
    p = Path(path)
    r = Path(root)
    candidates = (
        (p.resolve(), r.resolve()),
        (p.absolute(), r.absolute()),
        (p, r),
    )
    for pp, rr in candidates:
        try:
            return str(pp.relative_to(rr))
        except ValueError:
            continue
    try:
        return os.path.relpath(str(p), start=str(r))
    except ValueError:
        return str(p.resolve() if p.exists() else p)
