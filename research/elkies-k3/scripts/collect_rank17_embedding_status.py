from pathlib import Path
import re
import numpy as np

BASE = Path(__file__).resolve().parents[1]
R = BASE / "results"

JOBS = {
    "E29": {
        "log": R / "E29-short-pool-v2.log",
        "vectors": R / "E29-short-pool-v2-vectors.npy",
        "norms": R / "E29-short-pool-v2-norms.npy",
        "embedding": R / "rank17-E29-embedding-v2.log",
    },
    "rank19": {
        "log": R / "Nagao-rank19-short-pool-v2.log",
        "vectors": R / "E29-short-pool-v2-vectors.npy",  # fallback only
        "norms": R / "E29-short-pool-v2-norms.npy",
        "embedding": R / "rank17-Nagao-rank19-embedding.log",
    },
    "rank21": {
        "log": R / "Nagao-rank21-short-pool-v2.log",
        "vectors": R / "E29-short-pool-v2-vectors.npy",  # fallback only
        "norms": R / "E29-short-pool-v2-norms.npy",
        "embedding": R / "rank17-Nagao-rank21-embedding.log",
    },
}

# Because the generator currently writes generic E29-named files,
# try job-specific aliases first if you create/copy them later.
ALIASES = {
    "E29": (
        R / "E29-short-pool-v2-vectors.npy",
        R / "E29-short-pool-v2-norms.npy",
    ),
    "rank19": (
        R / "Nagao-rank19-short-pool-v2-vectors.npy",
        R / "Nagao-rank19-short-pool-v2-norms.npy",
    ),
    "rank21": (
        R / "Nagao-rank21-short-pool-v2-vectors.npy",
        R / "Nagao-rank21-short-pool-v2-norms.npy",
    ),
}


def last_matching(lines, prefix):
    for line in reversed(lines):
        if line.startswith(prefix):
            return line
    return None


def field(line, name, cast=float):
    if not line:
        return None

    m = re.search(
        rf"(?:^|\|){re.escape(name)}=([^|]+)",
        line
    )

    if not m:
        return None

    try:
        return cast(m.group(1))
    except Exception:
        return m.group(1)


def parse_log(path):
    if not path.exists():
        return {}

    lines = [
        x.strip()
        for x in path.read_text(errors="replace").splitlines()
        if x.strip()
    ]

    step = last_matching(lines, "STEP|")
    prune = last_matching(lines, "PRUNE_DONE|")
    checkpoint = last_matching(lines, "CHECKPOINT|")

    return {
        "step": field(step, "n", int),
        "unique": field(step, "unique", int),
        "accepted": field(step, "accepted", int),
        "rate": field(step, "rate", float),
        "walker_min": field(step, "walker_min", float),
        "walker_med": field(step, "walker_med", float),
        "walker_max": field(step, "walker_max", float),

        "prune_step": field(prune, "step", int),
        "prune_unique": field(prune, "unique", int),
        "active_cutoff": field(prune, "active_cutoff", float),

        "checkpoint_step": field(checkpoint, "step", int),
        "checkpoint_saved": field(checkpoint, "saved", int),
        "checkpoint_min": field(checkpoint, "min", float),
        "checkpoint_max": field(checkpoint, "max", float),
    }


def pool_stats(vpath, npath):
    if not vpath.exists() or not npath.exists():
        return None

    norms = np.load(
        npath,
        mmap_mode="r"
    )

    vectors = np.load(
        vpath,
        mmap_mode="r"
    )

    return {
        "count": len(norms),
        "dim": vectors.shape[1],
        "min": float(np.min(norms)),
        "q10": float(np.quantile(norms, 0.10)),
        "median": float(np.median(norms)),
        "q90": float(np.quantile(norms, 0.90)),
        "max": float(np.max(norms)),
    }


def embedding_stats(path):
    if not path.exists():
        return {}

    text = path.read_text(errors="replace")

    depths = [
        int(x)
        for x in re.findall(
            r"(?:BEST|DEPTH)\|[^\n]*?\|depth=(\d+)",
            text
        )
    ]

    m = re.search(
        r"best_depth\s*=\s*(\d+)",
        text
    )

    if m:
        depths.append(int(m.group(1)))

    full = "FOUND FULL EMBEDDING" in text

    return {
        "best_depth": max(depths) if depths else None,
        "full": full,
    }


print()
print("=" * 105)
print(
    f"{'job':<10}"
    f"{'step':>12}"
    f"{'unique':>12}"
    f"{'cutoff':>11}"
    f"{'rate/s':>11}"
    f"{'pool':>10}"
    f"{'norm min':>11}"
    f"{'median':>11}"
    f"{'norm max':>11}"
    f"{'depth':>7}"
)
print("=" * 105)

for name, job in JOBS.items():

    log = parse_log(job["log"])

    vpath, npath = ALIASES[name]

    pool = pool_stats(
        vpath,
        npath
    )

    emb = embedding_stats(
        job["embedding"]
    )

    def fmt(v, spec=".1f"):
        if v is None:
            return "-"
        if isinstance(v, int):
            return str(v)
        return format(v, spec)

    print(
        f"{name:<10}"
        f"{fmt(log.get('step')):>12}"
        f"{fmt(log.get('unique')):>12}"
        f"{fmt(log.get('active_cutoff'), '.3f'):>11}"
        f"{fmt(log.get('rate'), '.0f'):>11}"
        f"{fmt(pool['count'] if pool else None):>10}"
        f"{fmt(pool['min'] if pool else None, '.3f'):>11}"
        f"{fmt(pool['median'] if pool else None, '.3f'):>11}"
        f"{fmt(pool['max'] if pool else None, '.3f'):>11}"
        f"{fmt(emb.get('best_depth')):>7}"
    )

print("=" * 105)

print()
print("DETAILS")

for name, job in JOBS.items():

    print()
    print(f"[{name}]")

    log = parse_log(job["log"])

    for k in (
        "step",
        "unique",
        "accepted",
        "rate",
        "walker_min",
        "walker_med",
        "walker_max",
        "prune_step",
        "prune_unique",
        "active_cutoff",
        "checkpoint_step",
        "checkpoint_saved",
        "checkpoint_min",
        "checkpoint_max",
    ):
        print(f"  {k:20s} {log.get(k, '-')}")

    emb = embedding_stats(
        job["embedding"]
    )

    print(
        f"  {'embedding_depth':20s} "
        f"{emb.get('best_depth', '-')}"
    )

    print(
        f"  {'full_embedding':20s} "
        f"{emb.get('full', False)}"
    )

print()
