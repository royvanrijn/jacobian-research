"""Exact adapter from historical Curve302 V3 raw transcripts to exposure ledger.

The historical cloud JSON stores chart *definitions* (centre/index/mapping/search),
not quotient words.  This module therefore joins the chart schedule to recorded
rational-point results, reconstructs the seeded displayed-D basis from MW-state
snapshots, and recognizes every relevant recorded point in that basis with a
high-precision Neron--Tate solve followed by exact elliptic-curve group-law
verification.

Nothing here treats an absent record as a negative until the originating chart
search has explicit completion evidence.  Unsupported joins/point schemas remain
UNKNOWN rather than being converted into zero exposure.
"""
from __future__ import annotations

from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction as F
from hashlib import sha256
import json
import math
from pathlib import Path
import re
from typing import Any, Iterable, Sequence

try:
    from curve302_short_core_controls import integer, primitive, require
except Exception:  # tests may import this module standalone
    class InvalidEvidence(ValueError):
        pass
    def require(cond, message):
        if not cond:
            raise InvalidEvidence(message)
    def integer(v):
        require(not isinstance(v, (bool, float)), "expected exact integer")
        if isinstance(v, int): return v
        q = F(str(v)); require(q.denominator == 1, "expected integer"); return q.numerator
    def primitive(v):
        vals = tuple(int(x) for x in v); g = 0
        for x in vals: g = math.gcd(g, abs(x))
        require(g > 0, "zero vector")
        vals = tuple(x//g for x in vals)
        for x in vals:
            if x:
                return vals if x > 0 else tuple(-y for y in vals)

SCHEMA_MW = "elliptic-curves.mw-state.v1"
SCHEMA_POINT = "elliptic-curves.recorded-point-mod2-rank.v3"
SCHEMA_CLOUD = "elliptic-curves.retained-cloud-modl.v1"
LEDGER_SCHEMA = "curve302-chart-exposure-ledger.v1"

POINT_CONTAINER_HINTS = (
    "basis", "generators", "generator", "points", "point", "mw", "known",
    "new", "accepted", "recorded", "candidate",
)
INDEX_KEYS = (
    "chart_index", "cloud_index", "call_index", "search_index", "source_index",
    "index", "chart", "call", "attempt",
)
COMPLETE_KEYS = (
    "exhaustive", "complete", "search_complete", "replay_complete", "finished",
)
STATUS_KEYS = ("status", "state", "result")
GOOD_STATUS = {"pass", "passed", "success", "successful", "complete", "completed", "done", "ok"}
BAD_STATUS = {"timeout", "timed_out", "censored", "failed", "error", "killed", "interrupted"}
PATH_INDEX_PATTERNS = (
    r"(?:chart|cloud|call|attempt|index)[-_]?([0-9]+)",
)


def read_json(path: Path):
    def unique(pairs):
        out = {}
        for k, v in pairs:
            require(k not in out, f"duplicate JSON key {k} in {path}")
            out[k] = v
        return out
    def reject(value):
        raise ValueError(f"nonfinite JSON value {value} in {path}")
    return json.loads(Path(path).read_text(), object_pairs_hook=unique, parse_constant=reject)


def canonical_json(obj) -> str:
    return json.dumps(obj, sort_keys=True, separators=(",", ":"), ensure_ascii=True)


def hash_obj(obj) -> str:
    return sha256(canonical_json(obj).encode()).hexdigest()


def as_fraction(value) -> F | None:
    if isinstance(value, bool) or value is None:
        return None
    if isinstance(value, F):
        return value
    if isinstance(value, int):
        return F(value)
    if isinstance(value, str):
        text = value.strip()
        try:
            return F(text)
        except Exception:
            return None
    return None


def context_seed(path: Path, names: Sequence[str]) -> str | None:
    text = str(path).lower().replace("_", "-")
    found = [name for name in names if name.lower() in text]
    return found[0] if len(found) == 1 else None


def context_epoch(path: Path) -> int | None:
    """Bind epoch from *any* path component, not only the filename."""
    text = str(path).lower()
    found = []
    for pat in (r"(?:^|[/\\])epoch[-_]?([0-9]+)(?:[/\\]|$)",
                r"(?:^|[/\\])stage[-_]?([0-9]+)(?:[/\\]|$)",
                r"(?:^|[/\\])round[-_]?([0-9]+)(?:[/\\]|$)"):
        found += [int(x) for x in re.findall(pat, text)]
    vals = sorted(set(found))
    return vals[0] if len(vals) == 1 else None


def path_index(path: Path) -> int | None:
    found = []
    for part in Path(path).parts:
        low = part.lower()
        for pat in PATH_INDEX_PATTERNS:
            found += [int(x) for x in re.findall(pat, low)]
    vals = sorted(set(found))
    return vals[0] if len(vals) == 1 else None


def _parse_point_string(text: str):
    s = text.strip().strip("()[]")
    sep = ":" if ":" in s else "," if "," in s else None
    if sep is None:
        return None
    parts = [x.strip() for x in s.split(sep)]
    if len(parts) not in (2, 3):
        return None
    vals = [as_fraction(x) for x in parts]
    if any(v is None for v in vals):
        return None
    if len(vals) == 3:
        # Historical rational-point files use affine/projective z=1.  Do not
        # guess weighted projective conventions for any other z.
        if vals[2] != 1:
            return None
    return (vals[0], vals[1])


def parse_point(value):
    """Parse an exact affine rational EC point from common historical encodings."""
    if isinstance(value, str):
        return _parse_point_string(value)
    if isinstance(value, (list, tuple)) and len(value) in (2, 3):
        vals = [as_fraction(v) for v in value]
        if all(v is not None for v in vals):
            if len(vals) == 3 and vals[2] != 1:
                return None
            return (vals[0], vals[1])
    if isinstance(value, dict):
        # Direct affine forms.
        for xk, yk in (("x", "y"), ("X", "Y")):
            if xk in value and yk in value:
                x, y = as_fraction(value[xk]), as_fraction(value[yk])
                if x is not None and y is not None:
                    z = value.get("z", value.get("Z", 1))
                    zq = as_fraction(z)
                    if zq == 1:
                        return (x, y)
        # Explicit wrapper forms.
        for key in ("point", "elliptic_point", "affine_point", "coordinates", "xy"):
            if key in value:
                p = parse_point(value[key])
                if p is not None:
                    return p
    return None


def point_text(point) -> tuple[str, str]:
    x, y = point
    return (str(x), str(y))


def _walk(obj, path="$", depth=0):
    if depth > 16:
        return
    yield path, obj
    if isinstance(obj, dict):
        for k, v in obj.items():
            if isinstance(v, (dict, list, tuple)):
                yield from _walk(v, f"{path}.{k}", depth+1)
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            if isinstance(v, (dict, list, tuple)):
                yield from _walk(v, f"{path}[{i}]", depth+1)


def collect_points(obj, *, prefer_paths: Iterable[str] = ()):  # [(path, point)]
    found = []
    preferred = tuple(x.lower() for x in prefer_paths)
    for path, value in _walk(obj):
        p = parse_point(value)
        if p is None:
            continue
        score = sum(h in path.lower() for h in preferred)
        found.append((score, path, p))
    if preferred and any(score > 0 for score, _, _ in found):
        found = [row for row in found if row[0] > 0]
    # Avoid duplicates created by parsing both wrapper and child coordinate nodes.
    best = {}
    for score, path, p in found:
        key = point_text(p)
        row = best.get(key)
        if row is None or (score, -len(path)) > (row[0], -len(row[1])):
            best[key] = (score, path, p)
    return sorted(best.values(), key=lambda r: (-r[0], r[1], point_text(r[2])))


def point_lists(obj):
    """Return plausible MW basis/generator lists without interpreting arbitrary pairs."""
    rows = []
    for path, value in _walk(obj):
        if not isinstance(value, list) or not (10 <= len(value) <= 40):
            continue
        low = path.lower()
        if not any(h in low for h in POINT_CONTAINER_HINTS):
            continue
        parsed = [parse_point(v) for v in value]
        if all(p is not None for p in parsed):
            uniq = []
            seen = set()
            for p in parsed:
                k = point_text(p)
                if k not in seen:
                    seen.add(k); uniq.append(p)
            if len(uniq) >= 10:
                rows.append({"path": path, "points": tuple(uniq), "count": len(uniq)})
    rows.sort(key=lambda r: (0 if "basis" in r["path"].lower() else 1,
                             0 if "generator" in r["path"].lower() else 1,
                             r["count"], r["path"]))
    return rows


def collect_integer_fields(obj, keys=INDEX_KEYS):
    out = []
    wanted = {k.lower() for k in keys}
    for path, value in _walk(obj):
        if not isinstance(value, dict):
            continue
        for k, v in value.items():
            if k.lower() not in wanted:
                continue
            try:
                iv = integer(v)
            except Exception:
                continue
            if iv >= 0:
                out.append((f"{path}.{k}", iv))
    return out


def completion_evidence(obj):
    """Return True/False/None with conservative, explicit completion semantics."""
    if not isinstance(obj, dict):
        return None
    votes = []
    for path, value in _walk(obj):
        if not isinstance(value, dict):
            continue
        for key in COMPLETE_KEYS:
            if key in value and isinstance(value[key], bool):
                votes.append(value[key])
        for key in STATUS_KEYS:
            raw = value.get(key)
            if isinstance(raw, str):
                s = raw.strip().lower()
                if s in GOOD_STATUS: votes.append(True)
                if s in BAD_STATUS: votes.append(False)
        for key in ("timeout", "timed_out", "censored", "interrupted"):
            if isinstance(value.get(key), bool) and value[key]:
                votes.append(False)
        if "returncode" in value:
            try:
                rc = integer(value["returncode"])
                votes.append(rc == 0)
            except Exception:
                pass
    if not votes:
        return None
    if False in votes:
        return False
    return True


def infer_ainvariants(points: Sequence[tuple[F,F]]):
    """Infer [a1,a2,a3,a4,a6] exactly from rational points, then verify all."""
    from sympy import Matrix, Rational
    pts = list(points)
    require(len(pts) >= 5, "need five rational points to infer Weierstrass model")
    A, b = [], []
    for x, y in pts:
        # a1*x*y - a2*x^2 + a3*y - a4*x - a6 = x^3-y^2
        A.append([x*y, -x*x, y, -x, -1])
        b.append(x*x*x - y*y)
    Ms = Matrix([[Rational(q.numerator, q.denominator) for q in row] for row in A])
    bs = Matrix([Rational(q.numerator, q.denominator) for q in b])
    require(Ms.rank() >= 5, "point set does not determine a unique Weierstrass model")
    # Use exact least subsystem via gauss_jordan_solve on the overdetermined system.
    sol = Ms.gauss_jordan_solve(bs)
    vec, params = sol
    require(params.rows == 0, "Weierstrass model is not uniquely determined")
    ainv = tuple(F(int(v.p), int(v.q)) for v in vec)
    a1,a2,a3,a4,a6 = ainv
    for x,y in pts:
        require(y*y + a1*x*y + a3*y == x*x*x + a2*x*x + a4*x + a6,
                "inferred curve does not contain every basis point")
    return ainv


def on_curve(point, ainv):
    x,y = point; a1,a2,a3,a4,a6 = ainv
    return y*y + a1*x*y + a3*y == x*x*x + a2*x*x + a4*x + a6


def negate_point(point, ainv):
    x,y = point; a1,_,a3,_,_ = ainv
    return (x, -y-a1*x-a3)


def canonical_point_sign(point, ainv):
    neg = negate_point(point, ainv)
    return min(point_text(point), point_text(neg))


def inspect_raw_tree(raw_root: Path, names: Sequence[str]):
    """One-pass inventory specialized to the actual replay adapter inputs."""
    raw_root = Path(raw_root).resolve()
    files = sorted(p for p in raw_root.rglob("*.json") if p.is_file())
    require(files, "raw transcript tree contains no JSON")
    out = {"json_files": len(files), "schemas": Counter(), "mw_files": [], "point_files": [],
           "cloud_files": [], "chart_snapshots": defaultdict(list), "parse_errors": []}
    for file in files:
        try:
            obj = read_json(file)
        except Exception as exc:
            if len(out["parse_errors"]) < 30:
                out["parse_errors"].append({"file": str(file), "error": repr(exc)})
            continue
        schema = obj.get("schema") if isinstance(obj, dict) else None
        if schema is not None: out["schemas"][str(schema)] += 1
        rel = file.relative_to(raw_root)
        seed, epoch = context_seed(rel, names), context_epoch(rel)
        row = {"file": str(file), "seed": seed, "epoch": epoch, "schema": schema}
        if schema == SCHEMA_MW:
            pls = point_lists(obj)
            row["point_lists"] = [{"path": r["path"], "count": r["count"]} for r in pls[:8]]
            out["mw_files"].append(row)
        elif schema == SCHEMA_POINT:
            pts = collect_points(obj, prefer_paths=("point","candidate","recorded"))
            row["point_count"] = len(pts)
            row["indices"] = collect_integer_fields(obj)[:20]
            out["point_files"].append(row)
        elif schema == SCHEMA_CLOUD:
            out["cloud_files"].append(row)
        # Independently discover cumulative {charts:[...]} snapshots.
        if isinstance(obj, dict) and isinstance(obj.get("charts"), list) and obj["charts"] and all(isinstance(c, dict) for c in obj["charts"]):
            if seed is not None and epoch is not None:
                out["chart_snapshots"][(seed, epoch)].append({"file": str(file), "charts": len(obj["charts"])})
    out["schemas"] = dict(out["schemas"].most_common())
    out["chart_snapshots"] = {f"{s}:{e}": sorted(v, key=lambda r:(r["charts"],r["file"]))
                              for (s,e),v in sorted(out["chart_snapshots"].items())}
    return out


def load_stage_charts(raw_root: Path, names: Sequence[str], expected_epochs: dict[str,set[int]]):
    """Select the maximal cumulative chart snapshot for every historical stage."""
    raw_root = Path(raw_root).resolve()
    candidates = {}
    for file in sorted(raw_root.rglob("*.json")):
        if not file.is_file(): continue
        rel = file.relative_to(raw_root)
        seed, epoch = context_seed(rel, names), context_epoch(rel)
        if seed not in names or epoch is None or epoch not in expected_epochs.get(seed,set()):
            continue
        try: obj = read_json(file)
        except Exception: continue
        charts = obj.get("charts") if isinstance(obj, dict) else None
        if not (isinstance(charts, list) and charts and all(isinstance(c,dict) for c in charts)):
            continue
        key = (seed, epoch)
        row = candidates.get(key)
        count = len(charts)
        if row is None or count > row["count"]:
            candidates[key] = {"count": count, "file": file, "charts": charts,
                               "canonical": None}
        elif count == row["count"]:
            current = canonical_json(charts)
            if row["canonical"] is None:
                row["canonical"] = canonical_json(row["charts"])
            require(current == row["canonical"],
                    f"ambiguous maximal chart snapshots for {seed}/epoch-{epoch:02d}")
    selected = {}
    for seed, epochs in expected_epochs.items():
        for epoch in sorted(epochs):
            row = candidates.get((seed,epoch))
            require(row is not None, f"no chart snapshot bound to {seed}/epoch-{epoch:02d}")
            file, charts = row["file"], row["charts"]
            normalized = []
            for order, chart in enumerate(charts):
                idx_fields = collect_integer_fields(chart)
                direct = chart.get("index")
                try: idx = integer(direct)
                except Exception: idx = order
                normalized.append({"chart_id": f"{seed}:{epoch}:{idx}", "order": order, "index": idx,
                                   "centre": chart.get("centre"), "mapping": chart.get("mapping"),
                                   "search": chart.get("search"), "complete_flag": completion_evidence(chart.get("search", chart)),
                                   "source_file": str(file), "raw": chart})
            # Index must identify a chart inside the stage; duplicate index would make hit joins ambiguous.
            ids = [c["index"] for c in normalized]
            require(len(set(ids)) == len(ids), f"duplicate chart index in {seed}/epoch-{epoch:02d}")
            selected[(seed,epoch)] = normalized
    return selected


def derive_seeded_basis(raw_root: Path, names: Sequence[str]):
    """Recover common M17 and one exact seed point per exceptional direction from MW-state snapshots."""
    raw_root = Path(raw_root).resolve()
    by_seed = defaultdict(list)
    provenance = []
    for file in sorted(raw_root.rglob("*.json")):
        if not file.is_file(): continue
        rel = file.relative_to(raw_root); seed = context_seed(rel,names)
        if seed not in names: continue
        try: obj = read_json(file)
        except Exception: continue
        if not (isinstance(obj,dict) and obj.get("schema") == SCHEMA_MW): continue
        for row in point_lists(obj):
            if row["count"] in (17,18,19):
                s = frozenset(point_text(p) for p in row["points"])
                by_seed[seed].append({"file": str(file), "path": row["path"], "points": row["points"], "set": s, "count": row["count"],
                                      "epoch": context_epoch(file.relative_to(raw_root))})
                provenance.append({"seed":seed,"file":str(file),"path":row["path"],"count":row["count"]})
    missing = [s for s in names if not by_seed[s]]
    require(not missing, f"MW-state basis lists unavailable for seeds {missing[:5]}")
    # Prefer exact rank-18 seeded states.  The common 17 points must appear in at least one
    # candidate state for every seed.
    rank18 = {}
    for s in names:
        rows=[r for r in by_seed[s] if r["count"]==18]
        epochs=[r["epoch"] for r in rows if r["epoch"] is not None]
        if epochs:
            first=min(epochs); rows=[r for r in rows if r["epoch"]==first]
        rank18[s]=rows
    require(all(rank18[s] for s in names), "need one rank-18 MW-state snapshot per seed to derive M17+seed basis")
    # Frequency over seeds, not files, avoids cumulative snapshot multiplicity.
    point_seed_count = Counter()
    for s in names:
        union = set().union(*(r["set"] for r in rank18[s]))
        for p in union: point_seed_count[p] += 1
    common = {p for p,c in point_seed_count.items() if c == len(names)}
    require(len(common) == 17, f"expected exactly 17 common M17 points, found {len(common)}")
    chosen = {}; extras = {}
    for s in names:
        eligible = [r for r in rank18[s] if common <= r["set"] and len(r["set"]-common)==1]
        require(eligible, f"no rank-18 state for {s} equals common M17 plus one seed")
        # Multiple equivalent snapshots are allowed only if they have the same extra point.
        extra_sets = {next(iter(r["set"]-common)) for r in eligible}
        require(len(extra_sets)==1, f"ambiguous seeded exceptional point for {s}")
        extras[s] = next(iter(extra_sets)); chosen[s] = eligible[0]
    # Recover Fraction objects from selected snapshots.
    lookup = {}
    for r in [chosen[s] for s in names]:
        for p in r["points"]: lookup[point_text(p)] = p
    m17 = tuple(lookup[k] for k in sorted(common))
    exceptional = tuple(lookup[extras[s]] for s in names)
    ainv = infer_ainvariants(m17 + exceptional)
    require(all(on_curve(p,ainv) for p in m17+exceptional), "derived D basis contains off-curve point")
    return {"ainvariants": ainv, "m17": m17, "exceptional": exceptional,
            "chosen_states": {s:{"file":chosen[s]["file"],"path":chosen[s]["path"]} for s in names}}


def collect_strings(obj):
    out=[]
    for path,value in _walk(obj):
        if isinstance(value,dict):
            for k,v in value.items():
                if isinstance(v,str): out.append((f"{path}.{k}",v))
        elif isinstance(value,list):
            for i,v in enumerate(value):
                if isinstance(v,str): out.append((f"{path}[{i}]",v))
    return out

def resolve_json_reference(text: str, *, source_file: Path, raw_root: Path):
    raw=text.strip()
    if not raw.lower().endswith('.json'): return None
    candidates=[Path(raw)]
    if not Path(raw).is_absolute():
        candidates += [source_file.parent/raw, raw_root/raw]
    for p in candidates:
        try: q=p.resolve()
        except Exception: continue
        try: q.relative_to(raw_root.resolve())
        except Exception: continue
        if q.is_file(): return q
    # basename fallback is intentionally unique-only.
    matches=list(raw_root.rglob(Path(raw).name))
    return matches[0].resolve() if len(matches)==1 else None

def _candidate_hit_files(raw_root: Path, names: Sequence[str]):
    rows = []
    for file in sorted(Path(raw_root).rglob("*.json")):
        if not file.is_file(): continue
        rel = file.relative_to(raw_root); seed, epoch = context_seed(rel,names), context_epoch(rel)
        if seed not in names or epoch is None: continue
        try: obj = read_json(file)
        except Exception: continue
        if not (isinstance(obj,dict) and obj.get("schema") == SCHEMA_POINT): continue
        pts = collect_points(obj, prefer_paths=("point","candidate","recorded"))
        if not pts: continue
        idxs=[]
        for key in INDEX_KEYS:
            if key in obj:
                try:
                    iv=integer(obj[key])
                    if iv >= 0: idxs.append((f"$.{key}",iv))
                except Exception: pass
        pi = path_index(rel)
        if pi is not None: idxs.append(("$path",pi))
        if not idxs: idxs=collect_integer_fields(obj)
        # One file may contain wrapper+child aliases of the same point; collect_points dedupes.
        rows.append({"file":str(file),"seed":seed,"epoch":epoch,"payload":obj,
                     "points":[p for _,_,p in pts],"indices":idxs,
                     "complete_flag":completion_evidence(obj)})
    return rows


def join_hits_to_charts(raw_root: Path, names: Sequence[str], charts_by_stage):
    """Join recorded points to charts using independent exact identifiers/references."""
    raw_root=Path(raw_root).resolve()
    hits = _candidate_hit_files(raw_root,names)
    joined = defaultdict(list); unresolved = []
    # First consume points embedded directly in each chart search payload and exact
    # JSON result references.  The original schema probe intentionally did not
    # descend through `search`, so this is new evidence rather than reinterpretation
    # of its zero point-only count.
    loaded_refs={}
    for charts in charts_by_stage.values():
        for chart in charts:
            search=chart.get('search')
            for _,_,p in collect_points(search, prefer_paths=('point','found','hit','result','recorded')):
                joined[chart['chart_id']].append({'point':p,'source_file':chart['source_file']+'#search',
                                                  'record_complete_flag':completion_evidence(search)})
            for _,text in collect_strings(search):
                ref=resolve_json_reference(text,source_file=Path(chart['source_file']),raw_root=raw_root)
                if ref is None: continue
                if ref not in loaded_refs:
                    try: loaded_refs[ref]=read_json(ref)
                    except Exception: loaded_refs[ref]=None
                obj=loaded_refs[ref]
                if obj is None: continue
                for _,_,p in collect_points(obj, prefer_paths=('point','found','hit','result','recorded','candidate')):
                    joined[chart['chart_id']].append({'point':p,'source_file':str(ref),
                                                      'record_complete_flag':completion_evidence(obj)})
    for hit in hits:
        key=(hit["seed"],hit["epoch"]); charts=charts_by_stage.get(key)
        if not charts:
            unresolved.append({"file":hit["file"],"reason":"stage-not-in-roster"}); continue
        index_values = sorted(set(v for _,v in hit["indices"]))
        matches=[]
        for chart in charts:
            if chart["index"] in index_values or chart["order"] in index_values:
                matches.append(chart)
        # Stronger fallback: exact search/mapping object embedded in hit payload.
        if len(matches)!=1:
            embedded=[]
            for _,obj in _walk(hit["payload"]):
                if not isinstance(obj,dict): continue
                for chart in charts:
                    if obj == chart.get("search") or obj == chart.get("mapping"):
                        embedded.append(chart)
            uniq={c["chart_id"]:c for c in embedded}
            if len(uniq)==1: matches=list(uniq.values())
        if len(matches)!=1:
            unresolved.append({"file":hit["file"],"seed":hit["seed"],"epoch":hit["epoch"],
                               "index_values":index_values,"candidate_matches":[c["chart_id"] for c in matches[:20]],
                               "reason":"ambiguous-or-missing-chart-join"})
            continue
        chart=matches[0]
        for p in hit["points"]:
            joined[chart["chart_id"]].append({"point":p,"source_file":hit["file"],
                                               "record_complete_flag":hit["complete_flag"]})
    # Deduplicate wrapper aliases/repeated cumulative references per chart.
    for cid,rows in list(joined.items()):
        dedup={}
        for row in rows:
            key=point_text(row['point'])
            prev=dedup.get(key)
            if prev is None: dedup[key]=row
            elif prev.get('record_complete_flag') is not True and row.get('record_complete_flag') is True:
                dedup[key]=row
        joined[cid]=list(dedup.values())
    return joined, unresolved


@dataclass
class Recognition:
    point: tuple[F,F]
    coefficients: tuple[int,...]
    quotient_word: tuple[int,...]
    precision: int


class SageHeightRecognizer:
    """Recover D coordinates numerically, then prove them by exact Sage group equality."""
    def __init__(self, ainvariants, m17, exceptional, *, precisions=(128,256,512,1024)):
        from sage.all import EllipticCurve, QQ, RealField, matrix, vector, ZZ  # type: ignore
        self.EllipticCurve=EllipticCurve; self.QQ=QQ; self.RealField=RealField
        self.matrix=matrix; self.vector=vector; self.ZZ=ZZ
        self.E = EllipticCurve(QQ, [QQ(a.numerator)/a.denominator for a in ainvariants])
        self.ainvariants = ainvariants
        self.basis_points = tuple(self.E(QQ(x.numerator)/x.denominator, QQ(y.numerator)/y.denominator)
                                  for x,y in tuple(m17)+tuple(exceptional))
        require(len(self.basis_points)==31, "displayed D basis must have rank 31")
        self.precisions=tuple(int(p) for p in precisions)
        self._gram={}; self._hb={}

    def _height(self,P,prec):
        fn=getattr(P,"height",None)
        require(callable(fn), "Sage EC point has no canonical height method")
        try: return fn(precision=prec)
        except TypeError:
            try: return fn(prec=prec)
            except TypeError: return fn()

    def _gram_at(self,prec):
        if prec in self._gram: return self._gram[prec], self._hb[prec]
        RR=self.RealField(prec); B=self.basis_points
        hb=[RR(self._height(P,prec)) for P in B]
        G=self.matrix(RR,31,31)
        for i in range(31):
            G[i,i]=hb[i]
            for j in range(i):
                val=(RR(self._height(B[i]+B[j],prec))-hb[i]-hb[j])/2
                G[i,j]=G[j,i]=val
        self._gram[prec]=G; self._hb[prec]=hb
        return G,hb

    def _sage_point(self, point):
        x,y=point; QQ=self.QQ
        return self.E(QQ(x.numerator)/x.denominator, QQ(y.numerator)/y.denominator)

    def _sum(self, coeffs):
        P=self.E(0)
        for c,B in zip(coeffs,self.basis_points):
            if c: P += int(c)*B
        return P

    def recognize(self, point) -> Recognition | None:
        P=self._sage_point(point)
        for prec in self.precisions:
            RR=self.RealField(prec); G,hb=self._gram_at(prec)
            hp=RR(self._height(P,prec))
            b=[]
            for i,B in enumerate(self.basis_points):
                b.append((RR(self._height(P+B,prec))-hp-hb[i])/2)
            try: sol=G.solve_right(self.vector(RR,b))
            except Exception: continue
            coeffs=tuple(int(x.round()) for x in sol)
            exact=self._sum(coeffs)
            if exact == P:
                q=tuple(coeffs[17:])
                if any(q): q=primitive(q)
                return Recognition(point,coeffs,q,prec)
            if exact == -P:
                coeffs=tuple(-x for x in coeffs)
                q=tuple(coeffs[17:])
                if any(q): q=primitive(q)
                return Recognition(point,coeffs,q,prec)
        return None


def build_replay_ledger(raw_root: Path, names: Sequence[str], expected_epochs: dict[str,set[int]], *,
                        recognizer_factory=SageHeightRecognizer, max_unresolved_hits=None, progress=None):
    """Construct a normalized chart-exposure ledger from actual raw transcript schemas."""
    raw_root=Path(raw_root).resolve()
    charts=load_stage_charts(raw_root,names,expected_epochs)
    basis=derive_seeded_basis(raw_root,names)
    hit_map, unresolved=join_hits_to_charts(raw_root,names,charts)
    if max_unresolved_hits is not None:
        require(len(unresolved) <= max_unresolved_hits,
                f"{len(unresolved)} recorded-point files could not be bound uniquely to charts; first {unresolved[:5]}")
    unresolved_stages={(r.get('seed'),r.get('epoch')) for r in unresolved if r.get('seed') in names and r.get('epoch') is not None}
    global_unresolved=any(r.get('seed') not in names or r.get('epoch') is None for r in unresolved)
    rec=recognizer_factory(basis["ainvariants"],basis["m17"],basis["exceptional"])
    cache={}; outside=[]; recognized=0; unique_points=set()
    runs=[]; source_files=set()
    for seed in names:
        stages=[]
        for epoch in sorted(expected_epochs[seed]):
            normcharts=[]
            for chart in charts[(seed,epoch)]:
                exposures=[]; seen=set()
                for item in hit_map.get(chart["chart_id"],[]):
                    p=item["point"]; key=point_text(p); unique_points.add(key)
                    if key not in cache:
                        cache[key]=rec.recognize(p)
                        if progress is not None:
                            progress(len(cache), key, cache[key])
                    rr=cache[key]
                    if rr is None:
                        outside.append({"point":key,"source_file":item["source_file"],"chart_id":chart["chart_id"]})
                        continue
                    if not rr.quotient_word or rr.quotient_word in seen:
                        continue
                    seen.add(rr.quotient_word); recognized += 1
                    exposures.append({"quotient_word":list(rr.quotient_word),"primitive_quotient_word":list(rr.quotient_word),
                                      "parameter_height":None,"within_bound":True,
                                      "recognition_precision":rr.precision,"source_file":item["source_file"]})
                    source_files.add(item["source_file"])
                normcharts.append({"chart_id":chart["chart_id"],"order":chart["order"],"index":chart["index"],
                                   "score_band":None,"quartic_coefficients":chart.get("centre"),
                                   "search_bound":None,"complete": (None if global_unresolved or (seed,epoch) in unresolved_stages else chart.get("complete_flag")),
                                   "exposures":exposures})
                source_files.add(chart["source_file"])
            stages.append({"epoch":epoch,"charts":normcharts})
        runs.append({"seed":seed,"stages":stages})
    ledger={"schema":LEDGER_SCHEMA,"direction_ids":list(names),"runs":runs,
            "adapter": {"schema":"curve302-chart-replay-adapter.v1",
                        "basis": {"ainvariants":[str(x) for x in basis["ainvariants"]],
                                  "m17_points":len(basis["m17"]),"exceptional_points":len(basis["exceptional"]),
                                  "chosen_states":basis["chosen_states"]},
                        "recorded_hit_files":sum(len(v) for v in hit_map.values()),
                        "unique_recorded_points":len(unique_points),
                        "recognized_exposures":recognized,"outside_D_records":outside[:200],
                        "unresolved_hit_files":unresolved[:200]},
            "source_files":sorted(source_files)}
    return ledger


def replay_schema_audit(raw_root: Path, names: Sequence[str], expected_epochs: dict[str,set[int]]):
    """Cheap, non-promoting audit used before expensive height recognition."""
    raw_root=Path(raw_root).resolve()
    charts=load_stage_charts(raw_root,names,expected_epochs)
    basis=derive_seeded_basis(raw_root,names)
    hits, unresolved=join_hits_to_charts(raw_root,names,charts)
    status="PASS_REPLAY_SCHEMA_ADAPTER" if not unresolved else "PASS_REPLAY_SCHEMA_WITH_UNBOUND_POINT_FILES"
    return {"schema":"curve302-chart-replay-schema-audit.v1","status":status,
            "stages":len(charts),"charts":sum(len(v) for v in charts.values()),
            "bound_hit_records":sum(len(v) for v in hits.values()),"charts_with_hits":sum(bool(v) for v in hits.values()),
            "unresolved_hit_files":unresolved[:200],
            "basis":{"m17":len(basis["m17"]),"exceptional":len(basis["exceptional"]),
                     "ainvariants":[str(x) for x in basis["ainvariants"]],"chosen_states":basis["chosen_states"]},
            "boundary":"This audit proves path binding, MW-basis reconstruction, and recorded-hit/chart joins only. No quotient exposure is claimed until Sage height recognition and exact group equality succeed."}

# ---------------------------------------------------------------------------
# One-pass raw-tree index used by the real controller.  The earlier public
# helpers above remain independently testable and are useful for small fixtures.

def scan_replay_index(raw_root: Path, names: Sequence[str], expected_epochs: dict[str,set[int]]):
    raw_root=Path(raw_root).resolve(); require(raw_root.is_dir(),f'raw root missing: {raw_root}')
    files=sorted(p for p in raw_root.rglob('*.json') if p.is_file())
    require(files,'raw transcript tree contains no JSON')
    best_charts={}; mw=defaultdict(list); hits=[]; schemas=Counter(); parse_errors=[]
    examined=0
    for file in files:
        rel=file.relative_to(raw_root); seed,epoch=context_seed(rel,names),context_epoch(rel)
        try: obj=read_json(file)
        except Exception as exc:
            if len(parse_errors)<50: parse_errors.append({'file':str(file),'error':repr(exc)})
            continue
        examined += 1
        schema=obj.get('schema') if isinstance(obj,dict) else None
        if schema is not None: schemas[str(schema)] += 1
        # Cumulative stage chart snapshots: retain only the maximal list.
        if seed in names and epoch is not None and epoch in expected_epochs.get(seed,set()) and isinstance(obj,dict):
            charts=obj.get('charts')
            if isinstance(charts,list) and charts and all(isinstance(c,dict) for c in charts):
                key=(seed,epoch); row=best_charts.get(key); n=len(charts)
                canon=None
                if row is None or n > row['count']:
                    best_charts[key]={'count':n,'file':str(file),'charts':charts,'canonical':None}
                elif n == row['count']:
                    canon=canonical_json(charts)
                    if row['canonical'] is None: row['canonical']=canonical_json(row['charts'])
                    require(canon==row['canonical'],f'ambiguous maximal chart snapshots for {seed}/epoch-{epoch:02d}')
        if seed in names and schema == SCHEMA_MW:
            for row in point_lists(obj):
                if row['count'] in (17,18,19):
                    mw[seed].append({'file':str(file),'path':row['path'],'points':row['points'],
                                     'set':frozenset(point_text(p) for p in row['points']),
                                     'count':row['count'],'epoch':epoch})
        if seed in names and epoch is not None and schema == SCHEMA_POINT:
            pts=collect_points(obj,prefer_paths=('point','candidate','recorded'))
            if pts:
                idxs=[]
                for key in INDEX_KEYS:
                    if key in obj:
                        try:
                            iv=integer(obj[key])
                            if iv>=0: idxs.append((f'$.{key}',iv))
                        except Exception: pass
                pi=path_index(rel)
                if pi is not None: idxs.append(('$path',pi))
                if not idxs: idxs=collect_integer_fields(obj)
                hits.append({'file':str(file),'seed':seed,'epoch':epoch,'payload':obj,
                             'points':[p for _,_,p in pts],'indices':idxs,
                             'complete_flag':completion_evidence(obj)})
    missing=[(s,e) for s,epochs in expected_epochs.items() for e in sorted(epochs) if (s,e) not in best_charts]
    require(not missing,f'one-pass index missing chart snapshots for {missing[:10]}')
    return {'raw_root':raw_root,'files':len(files),'examined':examined,'schemas':dict(schemas.most_common()),
            'best_charts':best_charts,'mw':mw,'hits':hits,'parse_errors':parse_errors}


def charts_from_index(index, names: Sequence[str], expected_epochs: dict[str,set[int]]):
    out={}
    for seed,epochs in expected_epochs.items():
        for epoch in sorted(epochs):
            row=index['best_charts'][(seed,epoch)]; charts=row['charts']; normalized=[]
            for order,chart in enumerate(charts):
                try: idx=integer(chart.get('index',order))
                except Exception: idx=order
                normalized.append({'chart_id':f'{seed}:{epoch}:{idx}','order':order,'index':idx,
                                   'centre':chart.get('centre'),'mapping':chart.get('mapping'),'search':chart.get('search'),
                                   'complete_flag':completion_evidence(chart.get('search',chart)),
                                   'source_file':row['file'],'raw':chart})
            ids=[c['index'] for c in normalized]; require(len(set(ids))==len(ids),f'duplicate chart index in {seed}/{epoch}')
            out[(seed,epoch)]=normalized
    return out


def basis_from_index(index, names: Sequence[str]):
    by_seed=index['mw']; missing=[s for s in names if not by_seed.get(s)]
    require(not missing,f'MW-state basis lists unavailable for seeds {missing[:5]}')
    rank18={}
    for s in names:
        rows=[r for r in by_seed[s] if r['count']==18]
        epochs=[r['epoch'] for r in rows if r['epoch'] is not None]
        if epochs:
            first=min(epochs); rows=[r for r in rows if r['epoch']==first]
        rank18[s]=rows
    require(all(rank18[s] for s in names),'need one rank-18 MW-state snapshot per seed to derive M17+seed basis')
    freq=Counter()
    for s in names:
        union=set().union(*(r['set'] for r in rank18[s]))
        for p in union: freq[p]+=1
    common={p for p,c in freq.items() if c==len(names)}
    require(len(common)==17,f'expected exactly 17 common M17 points, found {len(common)}')
    chosen={}; extras={}; lookup={}
    for s in names:
        eligible=[r for r in rank18[s] if common<=r['set'] and len(r['set']-common)==1]
        require(eligible,f'no rank-18 state for {s} equals common M17 plus seed')
        es={next(iter(r['set']-common)) for r in eligible}; require(len(es)==1,f'ambiguous seed point for {s}')
        extras[s]=next(iter(es)); chosen[s]=eligible[0]
        for p in eligible[0]['points']: lookup[point_text(p)]=p
    m17=tuple(lookup[k] for k in sorted(common)); exceptional=tuple(lookup[extras[s]] for s in names)
    ainv=infer_ainvariants(m17+exceptional)
    return {'ainvariants':ainv,'m17':m17,'exceptional':exceptional,
            'chosen_states':{s:{'file':chosen[s]['file'],'path':chosen[s]['path']} for s in names}}


def join_hits_from_index(index, charts_by_stage, names: Sequence[str]):
    raw_root=index['raw_root']; joined=defaultdict(list); unresolved=[]; loaded_refs={}
    # Embedded/search-referenced results.
    for charts in charts_by_stage.values():
        for chart in charts:
            search=chart.get('search')
            for _,_,p in collect_points(search,prefer_paths=('point','found','hit','result','recorded')):
                joined[chart['chart_id']].append({'point':p,'source_file':chart['source_file']+'#search',
                                                  'record_complete_flag':completion_evidence(search)})
            for _,text in collect_strings(search):
                ref=resolve_json_reference(text,source_file=Path(chart['source_file']),raw_root=raw_root)
                if ref is None: continue
                if ref not in loaded_refs:
                    try: loaded_refs[ref]=read_json(ref)
                    except Exception: loaded_refs[ref]=None
                obj=loaded_refs[ref]
                if obj is None: continue
                for _,_,p in collect_points(obj,prefer_paths=('point','found','hit','result','recorded','candidate')):
                    joined[chart['chart_id']].append({'point':p,'source_file':str(ref),'record_complete_flag':completion_evidence(obj)})
    # Standalone recorded-point files.
    for hit in index['hits']:
        key=(hit['seed'],hit['epoch']); charts=charts_by_stage.get(key)
        if not charts:
            unresolved.append({'file':hit['file'],'reason':'stage-not-in-roster'}); continue
        vals=sorted(set(v for _,v in hit['indices'])); matches=[c for c in charts if c['index'] in vals or c['order'] in vals]
        if len(matches)!=1:
            embedded=[]
            for _,obj in _walk(hit['payload']):
                if not isinstance(obj,dict): continue
                for c in charts:
                    if obj==c.get('search') or obj==c.get('mapping'): embedded.append(c)
            uniq={c['chart_id']:c for c in embedded}
            if len(uniq)==1: matches=list(uniq.values())
        if len(matches)!=1:
            unresolved.append({'file':hit['file'],'seed':hit['seed'],'epoch':hit['epoch'],'index_values':vals,
                               'candidate_matches':[c['chart_id'] for c in matches[:20]],'reason':'ambiguous-or-missing-chart-join'})
            continue
        for p in hit['points']:
            joined[matches[0]['chart_id']].append({'point':p,'source_file':hit['file'],'record_complete_flag':hit['complete_flag']})
    for cid,rows in list(joined.items()):
        dedup={}
        for r in rows:
            k=point_text(r['point']); prev=dedup.get(k)
            if prev is None or (prev.get('record_complete_flag') is not True and r.get('record_complete_flag') is True): dedup[k]=r
        joined[cid]=list(dedup.values())
    return joined,unresolved


def replay_schema_audit_index(index, names: Sequence[str], expected_epochs: dict[str,set[int]]):
    charts=charts_from_index(index,names,expected_epochs); basis=basis_from_index(index,names); hits,unresolved=join_hits_from_index(index,charts,names)
    status='PASS_REPLAY_SCHEMA_ADAPTER' if not unresolved else 'PASS_REPLAY_SCHEMA_WITH_UNBOUND_POINT_FILES'
    return {'schema':'curve302-chart-replay-schema-audit.v1','status':status,'raw_json_files':index['files'],'files_examined':index['examined'],
            'schemas':index['schemas'],'stages':len(charts),'charts':sum(len(v) for v in charts.values()),
            'recorded_point_files':len(index['hits']),'bound_hit_records':sum(len(v) for v in hits.values()),
            'charts_with_hits':sum(bool(v) for v in hits.values()),'unresolved_hit_files':unresolved[:200],
            'basis':{'m17':17,'exceptional':14,'ainvariants':[str(x) for x in basis['ainvariants']],
                     'chosen_states':basis['chosen_states']},
            'boundary':'One-pass schema audit only; no quotient exposure is promoted until exact D-coordinate recognition succeeds.'}


def build_replay_ledger_index(index, names: Sequence[str], expected_epochs: dict[str,set[int]], *,
                              recognizer_factory=SageHeightRecognizer, max_unresolved_hits=None, progress=None):
    charts=charts_from_index(index,names,expected_epochs); basis=basis_from_index(index,names); hit_map,unresolved=join_hits_from_index(index,charts,names)
    if max_unresolved_hits is not None: require(len(unresolved)<=max_unresolved_hits,f'{len(unresolved)} unbound hit files exceed limit')
    unresolved_stages={(r.get('seed'),r.get('epoch')) for r in unresolved if r.get('seed') in names and r.get('epoch') is not None}
    global_unresolved=any(r.get('seed') not in names or r.get('epoch') is None for r in unresolved)
    rec=recognizer_factory(basis['ainvariants'],basis['m17'],basis['exceptional']); cache={}; outside=[]; recognized=0; unique=set(); runs=[]; source_files=set()
    for seed in names:
        stages=[]
        for epoch in sorted(expected_epochs[seed]):
            norm=[]
            for chart in charts[(seed,epoch)]:
                ex=[]; seen=set()
                for item in hit_map.get(chart['chart_id'],[]):
                    p=item['point']; k=point_text(p); unique.add(k)
                    if k not in cache:
                        cache[k]=rec.recognize(p)
                        if progress: progress(len(cache),k,cache[k])
                    rr=cache[k]
                    if rr is None or not rr.quotient_word or rr.quotient_word in seen: continue
                    seen.add(rr.quotient_word); recognized+=1
                    ex.append({'quotient_word':list(rr.quotient_word),'primitive_quotient_word':list(rr.quotient_word),
                               'parameter_height':None,'within_bound':True,'recognition_precision':rr.precision,
                               'source_file':item['source_file']}); source_files.add(item['source_file'])
                norm.append({'chart_id':chart['chart_id'],'order':chart['order'],'index':chart['index'],'score_band':None,
                             'search_bound':None,'complete':None if global_unresolved or (seed,epoch) in unresolved_stages else chart.get('complete_flag'),
                             'exposures':ex}); source_files.add(chart['source_file'])
            stages.append({'epoch':epoch,'charts':norm})
        runs.append({'seed':seed,'stages':stages})
    return {'schema':LEDGER_SCHEMA,'direction_ids':list(names),'runs':runs,
            'adapter':{'schema':'curve302-chart-replay-adapter.v1','basis':{'ainvariants':[str(x) for x in basis['ainvariants']],
                        'm17_points':17,'exceptional_points':14,'chosen_states':basis['chosen_states']},
                       'raw_files_examined':index['examined'],'recorded_point_files':len(index['hits']),
                       'unique_recorded_points':len(unique),'recognized_exposures':recognized,'outside_D_records':outside[:200],
                       'unresolved_hit_files':unresolved[:200]},'source_files':sorted(source_files)}
