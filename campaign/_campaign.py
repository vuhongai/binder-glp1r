"""Shared campaign glue: config loading, path resolution, notifications, stop."""
from __future__ import annotations

import os
import signal
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
BINDER = Path(os.environ.get("BINDER", "/DATA/binder"))

PATH_FIELDS = ("structure", "assembly", "conservation_csv")

_notifier = None
_stopping = False


def import_binder():
    if not BINDER.is_dir():
        raise SystemExit(f"binder package not found at {BINDER}; set $BINDER")
    if str(BINDER) not in sys.path:
        sys.path.insert(0, str(BINDER))


def log(msg: str) -> None:
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S')}] {msg}", flush=True)


def announce(title: str, body: str, tags: str = "dna",
             priority: str = "default") -> None:
    log(f"ntfy: {title} -- {body}")
    if _notifier:
        try:
            _notifier.send(title, body, tags=tags, priority=priority)
        except Exception as exc:                  # a phone failure is not a run failure
            log(f"ntfy failed: {exc}")


def load(config: str):
    """Load a config and make its repo-relative paths absolute."""
    import_binder()
    from binder.config import load_config
    cfg = load_config(config)
    t = cfg.target
    for field in PATH_FIELDS:
        value = getattr(t, field, None)
        if value and not os.path.isabs(value):
            setattr(t, field, str(ROOT / value))
    for species, path in list(t.species_fasta.items()):
        if path and not os.path.isabs(path):
            t.species_fasta[species] = str(ROOT / path)
    if cfg.results_root and not os.path.isabs(cfg.results_root):
        cfg.results_root = str(ROOT / cfg.results_root)
    return cfg


def missing_inputs(cfg) -> list[str]:
    t = cfg.target
    missing = [getattr(t, f) for f in PATH_FIELDS
               if getattr(t, f, None) and not os.path.isfile(getattr(t, f))]
    missing += [p for p in t.species_fasta.values() if not os.path.isfile(p)]
    if t.structure.endswith(".pdb"):
        contig = t.structure[:-4] + ".contig"
        if not os.path.isfile(contig):
            missing.append(contig + "  (without it the contig is the FIRST CHAIN "
                                    "ONLY -- diffusion would see a monomer)")
    return missing


def require_inputs(cfg) -> None:
    missing = missing_inputs(cfg)
    if missing:
        print("MISSING INPUTS:")
        for path in missing:
            print(f"  {path}")
        raise SystemExit(1)


def round1_stages(cfg) -> list[str]:
    """Round-1 stage directories that got as far as a shortlist."""
    import glob
    base = os.path.join(cfg.run_root, "round1")
    return sorted(d for d in glob.glob(os.path.join(base, "*"))
                  if os.path.isdir(os.path.join(d, "shortlist")))


def pick_stage(cfg, stage: str | None) -> str:
    """A named round-1 stage, or the first one with a shortlist."""
    if stage:
        path = os.path.join(cfg.run_root, "round1", stage)
        if not os.path.isdir(os.path.join(path, "shortlist")):
            raise SystemExit(f"no shortlist under {path}")
        return path
    stages = round1_stages(cfg)
    if not stages:
        raise SystemExit(f"no round-1 stage with a shortlist under "
                         f"{cfg.run_root}/round1 -- run campaign/run_round1.py "
                         f"first")
    return stages[0]


def run_paths(cfg, stage: str):
    """RunPaths rooted at a round-1 stage directory."""
    import_binder()
    from binder.paths import RunPaths
    return RunPaths(root=stage, round_name=os.path.basename(stage),
                    target_name=cfg.target.name)


def group_pids(exclude: int) -> list[int]:
    """PIDs in this process group except `exclude`: the tool subprocesses."""
    pids = []
    for entry in os.listdir("/proc"):
        if not entry.isdigit():
            continue
        try:
            with open(f"/proc/{entry}/stat") as fh:
                fields = fh.read().rsplit(")", 1)[1].split()
            if int(fields[2]) == os.getpgid(0) and int(entry) != exclude:
                pids.append(int(entry))
        except (OSError, IndexError, ValueError):
            continue
    return pids


def _handle_stop(signum, _frame) -> None:
    global _stopping
    if _stopping:
        os._exit(130)
    _stopping = True
    announce("stopped", f"signal {signum} -- killed the tool subprocesses; "
                        f"relaunch to resume (finished work is kept)",
             tags="stopwatch", priority="high")
    me = os.getpid()
    for sig in (signal.SIGTERM, signal.SIGKILL):
        for pid in group_pids(exclude=me):
            try:
                os.kill(pid, sig)
            except OSError:
                pass
        time.sleep(5)
    os._exit(130)


def tee_to(path: str) -> None:
    """Send this process's AND every subprocess's stdout/stderr to `path`.

    os.dup2 on fd 1/2 rather than replacing sys.stdout, because the package
    shells out for every tool: a Python-level wrapper would capture the driver's
    own prints and lose RFdiffusion, MPNN and AF2 entirely. The terminal goes
    quiet, which is what a multi-day detached run wants; follow it with tail -f.
    Appended, so a resumed run adds to the same log.
    """
    os.makedirs(os.path.dirname(path), exist_ok=True)
    fd = os.open(path, os.O_WRONLY | os.O_CREAT | os.O_APPEND, 0o644)
    sys.stdout.flush()
    sys.stderr.flush()
    os.dup2(fd, 1)
    os.dup2(fd, 2)
    os.close(fd)


def setup(cfg, log: str | None = None):
    """Log file, notifications, the GPU-hold notice, and the stop handler.

    Returns the Notifier; wrap the run in `notifier.arm(label, log)` so a
    failure notification carries the tail of that log.
    """
    global _notifier
    from binder import runner
    from binder.notify import Notifier
    if log:
        print(f"log: {log}", flush=True)
        tee_to(log)
        log_start(cfg)
    signal.signal(signal.SIGTERM, _handle_stop)
    signal.signal(signal.SIGINT, _handle_stop)
    _notifier = Notifier.from_config(cfg)
    runner.gpu_notify = lambda gpu, free: announce(
        f"{cfg.target.name} GPU{gpu} holding",
        f"waiting for GPU {gpu} ({free} MiB free) -- stop the LLM server",
        tags="hourglass_flowing_sand")
    return _notifier


def log_start(cfg) -> None:
    log("=" * 72)
    log(f"{cfg.target.name}  pid {os.getpid()}  cwd {os.getcwd()}")
    log(f"binder {BINDER}")
