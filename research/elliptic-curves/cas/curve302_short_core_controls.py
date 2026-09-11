"""Exact lattice checks and a conditional short-vector trajectory control.

Integer lattices use rows. A QuotientMap is a surjective homomorphism Z^n ->
Z^(n-r) with kernel the *actual* primitive prefix, not an implicit saturation.
Only the already-enumerated, sign-canonical primitive directions are sampled.
"""
from __future__ import annotations

from array import array
from bisect import bisect_right
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import math
from pathlib import Path
import random
from typing import Iterable, Sequence

from sympy import Matrix, ZZ
from sympy.matrices.normalforms import hermite_normal_form, smith_normal_form
from sympy.polys.matrices import DomainMatrix
from sympy.polys.matrices.normalforms import smith_normal_decomp


class InvalidEvidence(ValueError):
    pass


def require(condition, message):
    if not condition:
        raise InvalidEvidence(message)


def integer(value):
    """Accept integer JSON values/text, never silently truncate a float/rational."""
    require(not isinstance(value, (bool, float)), "expected exact integer, not bool/float")
    if isinstance(value, int):
        return value
    q = F(str(value))
    require(q.denominator == 1, "nonintegral coefficient")
    return q.numerator


def primitive(values):
    v = tuple(integer(x) for x in values)
    g = math.gcd(*v) if v else 0
    if not g:
        return v
    sign = 1 if next(x for x in v if x) > 0 else -1
    return tuple(sign * x // g for x in v)


def qnorm(form, v):
    """Exact norm; sparse summation is important for the supplied basis words."""
    support = [(i, integer(x)) for i, x in enumerate(v) if x]
    return sum((form[i][i] * x*x + sum(2*form[i][j]*x*y for j, y in support if j > i)
                for i, x in support), F(0))


def dot(a, b):
    return sum(x*y for x, y in zip(a, b))


def row_basis(rows, n):
    """Small exact rational echelon basis, stored as primitive integer rows."""
    pivots = {}
    for raw in rows:
        v = tuple(integer(x) for x in raw)
        require(len(v) == n, "row width mismatch")
        for p in sorted(pivots):
            if v[p]:
                b = pivots[p]
                v = primitive(tuple(b[p]*x - v[p]*y for x, y in zip(v, b)))
        if any(v):
            v = primitive(v)
            pivots[next(i for i, x in enumerate(v) if x)] = v
    return tuple(pivots[p] for p in sorted(pivots))


def rational_rank(rows, n):
    return len(row_basis(rows, n))


def row_hnf(rows, n):
    if not rows:
        return ()
    h = hermite_normal_form(Matrix(rows).T).T
    require(h.cols == n, "HNF width mismatch")
    return tuple(tuple(int(x) for x in h.row(i)) for i in range(h.rows))


def smith_packet(rows, n, *, witnesses=True):
    """Smith invariants plus independently checked integer generation witnesses."""
    require(rows and all(len(r) == n for r in rows), "empty or malformed Smith input")
    a = Matrix([[integer(v) for v in row] for row in rows])
    if witnesses:
        dm = DomainMatrix.from_Matrix(a).convert_to(ZZ)
        dd, uu, vv = smith_normal_decomp(dm)
        d, u, v = dd.to_Matrix(), uu.to_Matrix(), vv.to_Matrix()
        require(u*a*v == d, "Smith transformation identity failed")
        require(abs(u.det()) == 1 and abs(v.det()) == 1, "Smith transforms not unimodular")
    else:
        # The tied shell may be large. Never allocate a K-by-K transform for
        # that secondary index check; Smith invariants alone suffice.
        d = smith_normal_form(a, domain=ZZ)
    invariants = [abs(int(d[i, i])) for i in range(min(d.shape)) if d[i, i]]
    rank = len(invariants)
    require(all(b % a == 0 for a, b in zip(invariants, invariants[1:])), "Smith divisibility failed")
    index = math.prod(invariants)
    h = hermite_normal_form(a.T)
    if rank == n:
        require(h.rows == h.cols == n and abs(int(h.det())) == index, "HNF/Smith index disagreement")
    result = {"rows": len(rows), "rank": rank, "smith_invariants": invariants,
              "index_in_saturation": str(index), "index_in_ambient": str(index) if rank == n else None,
              "generates_ambient_integrally": rank == n and index == 1,
              "hnf_rows": row_hnf(rows, n), "axis_multiple_witnesses": []}
    if witnesses and rank == n:
        # c*A=m*e_j iff (m*e_j*V)_i is divisible by d_i. U records the
        # integer combination in the ORIGINAL ordered generator list.
        for j in range(n):
            m = math.lcm(*(abs(int(d[i, i])) // math.gcd(abs(int(d[i, i])), abs(int(v[j, i])))
                           for i in range(n)))
            t = Matrix(1, a.rows, lambda _, i: m*v[j, i]/d[i, i] if i < n else 0)
            require(all(x.q == 1 for x in t), "nonintegral Smith witness")
            coeff = t*u
            expected = Matrix(1, n, lambda _, k: m if k == j else 0)
            require(coeff*a == expected, "axis reconstruction witness failed")
            result["axis_multiple_witnesses"].append({"axis": j, "smallest_positive_multiple": str(m),
                                                       "coefficients": [int(x) for x in coeff]})
    return result


def egcd(a, b):
    old_r, r, old_s, s, old_t, t = abs(a), abs(b), 1, 0, 0, 1
    while r:
        q = old_r // r
        old_r, r = r, old_r-q*r
        old_s, s = s, old_s-q*s
        old_t, t = t, old_t-q*t
    return old_r, old_s if a >= 0 else -old_s, old_t if b >= 0 else -old_t


class QuotientMap:
    """Exact primitive-extension tracker; no lattice is silently saturated."""
    def __init__(self, n, equations=None):
        self.n = n
        self.equations = (tuple(tuple(int(i == j) for j in range(n)) for i in range(n))
                          if equations is None else tuple(tuple(r) for r in equations))

    @property
    def rank(self):
        return self.n - len(self.equations)

    def image(self, v):
        require(len(v) == self.n, "candidate width mismatch")
        return tuple(dot(a, v) for a in self.equations)

    def contains(self, rows):
        return all(not any(self.image(row)) for row in rows)

    def admission(self, v):
        z = self.image(v)
        if not any(z):
            return "DEPENDENT", z
        return ("PRIMITIVE" if math.gcd(*z) == 1 else "NONPRIMITIVE"), z

    def extend(self, v):
        kind, z0 = self.admission(v)
        require(kind == "PRIMITIVE", "extension is " + kind)
        b = list(self.equations)
        z = list(z0)
        p = next(i for i, x in enumerate(z) if x)
        b[0], b[p] = b[p], b[0]
        z[0], z[p] = z[p], z[0]
        for i in range(1, len(b)):
            if not z[i]:
                continue
            a, c = z[0], z[i]
            g, s, t = egcd(a, c)
            first, other = b[0], b[i]
            b[0] = tuple(s*x+t*y for x, y in zip(first, other))
            b[i] = tuple((-c//g)*x+(a//g)*y for x, y in zip(first, other))
            z[0], z[i] = g, 0
        require(abs(z[0]) == 1 and all(x == 0 for x in z[1:]), "primitive quotient reduction failed")
        result = QuotientMap(self.n, row_hnf(b[1:], self.n))
        require(result.rank == self.rank+1 and result.contains((v,)), "quotient update failed")
        return result


def common_rank(maps, n):
    require(maps, "no maps for common intersection")
    return n-rational_rank([r for q in maps for r in q.equations], n)


def verify_common_core(maps, basis, n):
    """For saturated inputs, containment plus rank proves exact integral equality."""
    require(all(q.contains(basis) for q in maps), "reported core not in every prefix")
    require(common_rank(maps, n) == len(basis) == rational_rank(basis, n), "common-core rank mismatch")
    if basis:
        snf = smith_normal_form(Matrix(basis), domain=ZZ)
        require(math.prod(abs(int(snf[i, i])) for i in range(len(basis))) == 1, "reported common core is not saturated")


class Vocabulary:
    """Compact read-only vocabulary: int64 storage, unbounded-integer arithmetic.

    Norms are streamed, not retained for 1.3 million rows. Only shell boundaries,
    requested-vector positions, and the first 29 plus boundary shell are kept.
    """
    def __init__(self, n):
        self.n = n
        self.values = array("q")
        self.shells = array("Q")
        self.positions = {}
        self.wanted_norms = {}
        self.first = []
        self.last_norm = None

    def __len__(self):
        return len(self.values)//self.n

    def vector(self, i):
        require(0 <= i < len(self), "vocabulary index outside range")
        return tuple(self.values[i*self.n:(i+1)*self.n])

    def shell(self, index):
        k = bisect_right(self.shells, index)-1
        require(0 <= index < len(self), "shell index outside range")
        return int(self.shells[k]), int(self.shells[k+1])

    def band(self, index):
        """Decades of 1-based worst tied rank, expanded to complete norm shells."""
        _, hi = self.shell(index)
        upper = 10
        lower = 1
        while hi > upper:
            lower, upper = upper+1, upper*10
        lo0, hi0 = lower-1, min(upper, len(self))
        lo, _ = self.shell(lo0)
        _, end = self.shell(hi0-1)
        return {"nominal_rank_low": lower, "nominal_rank_high": upper,
                "index_start": lo, "index_stop": end, "size": end-lo,
                "tie_expanded": (lo, end) != (lo0, hi0),
                "ceiling_truncated": upper > len(self)}

    @classmethod
    def load(cls, path, n, count, bound, wanted=(), form=None):
        obj = cls(n)
        wanted = set(map(tuple, wanted))
        previous = None
        shell29_norm = None
        with Path(path).open() as stream:
            require(stream.readline().rstrip("\n") == "norm\tvector", "wrong enumeration TSV header")
            for i, line in enumerate(stream):
                fields = line.rstrip("\n").split("\t")
                require(len(fields) == 2, "malformed enumeration row")
                norm = F(fields[0])
                v = tuple(integer(x) for x in fields[1].split(","))
                require(len(v) == n and any(v) and primitive(v) == v, "invalid primitive direction")
                require(norm > 0 and norm <= bound, "norm outside certified ceiling")
                require(previous is None or (norm, v) > previous, "enumeration reordered or duplicated")
                if previous is None or norm != previous[0]:
                    obj.shells.append(i)
                try:
                    obj.values.extend(v)
                except OverflowError as exc:
                    raise InvalidEvidence("vocabulary coordinate exceeds exact int64 STORAGE; no truncation allowed") from exc
                if i == 28:
                    shell29_norm = norm
                if i < 29 or norm == shell29_norm:
                    obj.first.append((norm, v))
                if v in wanted:
                    require(v not in obj.positions, "duplicate requested direction")
                    obj.positions[v] = i
                    obj.wanted_norms[v] = norm
                if form is not None and (i < 29 or norm == shell29_norm or v in wanted):
                    require(qnorm(form, v) == norm, "witness norm disagrees with supplied quotient form")
                previous = norm, v
        obj.shells.append(len(obj))
        require(len(obj) == count and count > 0, "enumeration direction count mismatch")
        require(wanted <= obj.positions.keys(), "acquisition missing from certified vocabulary")
        obj.last_norm = previous[0]
        require(obj.last_norm == bound, "enumeration did not reach declared observed ceiling")
        return obj


def band_for_rank(vocab, index):
    lo, hi = vocab.shell(index)
    return {"observed_rank_best": lo+1, "observed_rank_worst": hi, **vocab.band(index)}


def first_direction_index(vocab):
    count = min(29, len(vocab))
    literal = smith_packet([v for _, v in vocab.first[:count]], vocab.n)
    complete_shell = smith_packet([v for _, v in vocab.first], vocab.n, witnesses=False)
    return {"status": "PASS_EXACT_INDEX" if count == 29 else "INSUFFICIENT_VOCABULARY_FOR_29",
            "literal_first_29": literal, "whole_29th_norm_shell": complete_shell,
            "requested_directions": 29, "available_literal_directions": count,
            "boundary_norm": str(vocab.first[count-1][0]),
            "tie_boundary": len(vocab.first) != count,
            "boundary": "Index concerns the literal integer span; no saturation is applied. Tied 29th shell is separately included."}


def corrected_basins(rows):
    """Separate automatic hits; aggregate ONLY deficit-one cases, with counts."""
    rows = [dict(row) for row in rows]
    seen = set()
    for row in rows:
        for key in ("landmark_dimension", "deficit_before", "hits_at_actual_norm", "eligible_at_actual_norm",
                    "hits_at_complete_bound", "eligible_at_complete_bound"):
            row[key] = integer(row[key])
        key = row["landmark_dimension"], row["seed"]
        require(key not in seen, "duplicate basin row")
        seen.add(key)
        require(integer(row["deficit_before"]) in (0, 1), "unexpected basin deficit")
        for suffix in ("actual_norm", "complete_bound"):
            num = integer(row["hits_at_"+suffix]); den = integer(row["eligible_at_"+suffix])
            require(0 <= num <= den and den > 0, "invalid basin counts")
            if row["deficit_before"] == 0:
                require(num == den, "already-contained row is not automatic")
        require(row["hits_at_actual_norm"] <= row["hits_at_complete_bound"] and
                row["eligible_at_actual_norm"] <= row["eligible_at_complete_bound"], "basin bound counts decrease")

    def summarize(block):
        good = [r for r in block if r["deficit_before"] == 1]
        stats = {"total_cases": len(block), "already_contained_excluded": len(block)-len(good),
                 "nontrivial_cases": len(good), "bounds": {}}
        for suffix in ("actual_norm", "complete_bound"):
            counts = [(integer(r["hits_at_"+suffix]), integer(r["eligible_at_"+suffix])) for r in good]
            num, den = sum(a for a, _ in counts), sum(b for _, b in counts)
            stats["bounds"][suffix] = {"hits": num, "eligible": den,
                                       "pooled_fraction": str(F(num, den)) if den else None,
                                       "mean_case_fraction": str(sum((F(a, b) for a, b in counts), F(0))/len(counts)) if counts else None}
        return stats
    return {"status": "PASS_NONTRIVIAL_BASIN_REAGGREGATION", "overall": summarize(rows),
            "by_landmark": [{"landmark_dimension": d, **summarize([r for r in rows if r["landmark_dimension"] == d])}
                            for d in sorted({integer(r["landmark_dimension"]) for r in rows})],
            "nontrivial_rows": [r for r in rows if r["deficit_before"] == 1],
            "excluded_already_contained": [{"seed": r["seed"], "landmark_dimension": r["landmark_dimension"]}
                                           for r in rows if r["deficit_before"] == 0],
            "boundary": "Exact reaggregation of hash-verified source counts, not a new basin enumeration. Denominators include independent directions with potentially nonprimitive extensions; these are saturated-span counts, not the primitive-only sampling law or V3 probabilities."}


def seeded_rng(master, panel, seed_name):
    label = f"curve302-short-core-control.v1|{master}|{panel}|{seed_name}".encode()
    return random.Random(int.from_bytes(sha256(label).digest(), "big"))


def sample_extension(vocab, state, band, rng, *, max_draws=4096, rescue_after=64, rescue_scan_limit=4096):
    """Uniform conditional primitive extension, without widening the frozen band.

    Exact rescue scanning after failed proposals is also uniform on admissible
    choices, so it does not change the conditional sampling law. Empty bands
    and exhausted large-band draw budgets are *different censored outcomes*.
    """
    lo, hi = band["index_start"], band["index_stop"]
    require(0 <= lo < hi <= len(vocab), "empty/malformed frozen band")
    rejected = Counter()
    for draw in range(1, max_draws+1):
        index = rng.randrange(lo, hi)
        kind, _ = state.admission(vocab.vector(index))
        if kind == "PRIMITIVE":
            return {"status": "ACCEPTED", "index": index, "draws": draw, "rejected": dict(sorted(rejected.items())),
                    "selection": "uniform_rejection", "band_size": hi-lo}
        rejected[kind] += 1
        if draw == rescue_after and hi-lo <= rescue_scan_limit:
            accepted, choice = 0, None
            for index in range(lo, hi):
                if state.admission(vocab.vector(index))[0] == "PRIMITIVE":
                    accepted += 1
                    if rng.randrange(accepted) == 0:
                        choice = index
            return {"status": "ACCEPTED" if accepted else "CENSORED_EMPTY_ADMISSIBLE_BAND",
                    "index": choice, "draws": draw, "rejected": dict(sorted(rejected.items())),
                    "selection": "uniform_reservoir_rescue", "exact_admissible_count": accepted, "band_size": hi-lo}
    return {"status": "CENSORED_DRAW_CAP", "index": None, "draws": max_draws,
            "rejected": dict(sorted(rejected.items())), "selection": "bounded_rejection", "band_size": hi-lo}


def track_observed(runs, n, source_cores):
    maps = {}
    prefix_count = 0
    for run in runs:
        seed = tuple(int(i == run["seed_index"]) for i in range(n))
        state = QuotientMap(n).extend(seed)
        states = {1: state}
        for event in run["events"]:
            state = state.extend(event["word"])
            states[state.rank] = state
        require(state.rank == run["final_dimension"], "observed final rank mismatch")
        maps[run["seed"]] = states
        prefix_count += len(states)
    require(len({r["seed"] for r in runs}) == len(runs), "duplicate trajectory seed")
    for row in source_cores:
        d = integer(row["quotient_dimension"])
        participating = [maps[r["seed"]][d] for r in runs if r["final_dimension"] >= d]
        require(len(participating) == integer(row["run_count"]), "core participant count mismatch")
        verify_common_core(participating, row["basis"], n)
    return maps, prefix_count


def panel_metrics(states_by_run, runs, source_cores, n, axes):
    """Keep censored panels in denominators; late 13-run cohorts are explicit."""
    cores = []
    for core in source_cores:
        d = integer(core["quotient_dimension"])
        participants = [r["seed"] for r in runs if r["final_dimension"] >= d]
        present = [states_by_run[s].get(d) for s in participants]
        complete = all(s is not None for s in present)
        known = [s for s in present if s is not None]
        # A known counterexample settles non-containment even if another run
        # is censored. Equality/rank is left unknown without the full cohort.
        hit = False if any(not s.contains(core["basis"]) for s in known) else True if complete else None
        rank = common_rank(present, n) if complete else None
        cores.append({"dimension": d, "participants": len(participants), "observed_core_rank": len(core["basis"]),
                      "complete_cohort": complete, "common_rank": rank, "contains_observed_core": hit,
                      "rank_at_least_observed": rank >= len(core["basis"]) if rank is not None else None,
                      "equals_observed_core": (hit and rank == len(core["basis"])) if rank is not None else None,
                      "primary": d <= n-2 and len(core["basis"]) > 0,
                      "trivial_target": len(core["basis"]) == 0 or d == n})
    axis_rows = []
    for name, index in axes:
        target = (tuple(int(j == index) for j in range(n)),)
        per_run = []
        for run in runs:
            states = states_by_run[run["seed"]]
            first = next((d for d in sorted(states) if states[d].contains(target)), None)
            per_run.append({"seed": run["seed"], "first_dimension": first,
                            "attained_dimension": max(states), "planned_final_dimension": run["final_dimension"]})
        first_all = max(r["first_dimension"] for r in per_run) if all(r["first_dimension"] is not None for r in per_run) else None
        axis_rows.append({"axis": name, "all_runs_first_containment": first_all, "per_run": per_run})
    return {"cores": cores, "axes": axis_rows}


def simulate_panel(vocab, runs, source_cores, n, axes, master, panel, policy):
    states_by_run = {}
    receipts = []
    for run in runs:
        rng = seeded_rng(master, panel, run["seed"])
        state = QuotientMap(n).extend(tuple(int(i == run["seed_index"]) for i in range(n)))
        states = {1: state}
        events = []
        for step, original in enumerate(run["events"]):
            record = sample_extension(vocab, state, original["band"], rng,
                                      max_draws=policy["max_draws"], rescue_after=policy["rescue_after"],
                                      rescue_scan_limit=policy["rescue_scan_limit"])
            record.update({"step": step, "before": state.rank})
            events.append(record)
            if record["status"] != "ACCEPTED":
                break
            v = vocab.vector(record["index"])
            state = state.extend(v)
            states[state.rank] = state
        states_by_run[run["seed"]] = states
        receipts.append({"seed": run["seed"], "planned_final_dimension": run["final_dimension"],
                         "attained_dimension": state.rank, "complete": state.rank == run["final_dimension"], "events": events})
    return {"panel": panel, "all_runs_complete": all(r["complete"] for r in receipts), "runs": receipts,
            "metrics": panel_metrics(states_by_run, runs, source_cores, n, axes)}


def boolean_reference(values):
    n = len(values)
    yes = sum(v is True for v in values); no = sum(v is False for v in values); unknown = sum(v is None for v in values)
    require(n > 0 and yes+no+unknown == n, "invalid reference observations")
    return {"panels": n, "true": yes, "false": no, "unknown": unknown,
            "empirical_fraction_lower": str(F(yes, n)), "empirical_fraction_upper": str(F(yes+unknown, n)),
            "meaning": "Censoring bounds on this finite ensemble; NOT confidence bounds or a p-value."}


def summarize_panels(panels, observed):
    require(panels, "no panels")
    core_rows = []
    for i, expected in enumerate(observed["cores"]):
        block = [p["metrics"]["cores"][i] for p in panels]
        require(all(x["dimension"] == expected["dimension"] for x in block), "panel core alignment failed")
        core_rows.append({"dimension": expected["dimension"], "participants": expected["participants"],
                          "observed_common_rank": expected["common_rank"], "primary": expected["primary"],
                          "trivial_target": expected["trivial_target"],
                          "common_rank_histogram": dict(sorted(Counter(str(x["common_rank"]) for x in block).items())),
                          "contains_observed_core": boolean_reference([x["contains_observed_core"] for x in block]),
                          "equals_observed_core": boolean_reference([x["equals_observed_core"] for x in block]),
                          "rank_at_least_observed": boolean_reference([x["rank_at_least_observed"] for x in block])})
    axes = []
    for i, expected in enumerate(observed["axes"]):
        deadline = expected["all_runs_first_containment"]
        vals = []
        for panel in panels:
            rows = panel["metrics"]["axes"][i]["per_run"]
            if deadline is None:
                vals.append(None)
                continue
            if any(r["attained_dimension"] >= deadline and (r["first_dimension"] is None or r["first_dimension"] > deadline) for r in rows):
                vals.append(False)
            elif all(r["first_dimension"] is not None and r["first_dimension"] <= deadline for r in rows):
                vals.append(True)
            else:
                vals.append(None)
        axes.append({"axis": expected["axis"], "observed_all_runs_dimension": deadline,
                     "all_runs_by_observed_dimension": boolean_reference(vals)})
    return {"panels": len(panels), "fully_completed_panels": sum(p["all_runs_complete"] for p in panels),
            "censored_runs": sum(not r["complete"] for p in panels for r in p["runs"]),
            "cores": core_rows, "axes": axes,
            "boundary": "A conditional short-vector-only reference ensemble on a known quotient lattice, matched to retrospective height bands. Observed cores were selected from these same histories. No formal significance claim, independence of acquisitions, prospective selector, or new rank theorem follows."}
