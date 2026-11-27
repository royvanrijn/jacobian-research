"""Directed slab with full simplex clipping and stronger annular slack."""

import argparse
import hashlib
import json
import math
from fractions import Fraction as Q
from pathlib import Path
import sys
import time

import numpy as np


LOWER = Q(221, 500)
UPPER = Q(8859, 20000)
ORIGINAL_SHA256 = "c974ff60f8c1bf4723cd22d1a3074f932437b2303eb8219bcd29741a950aba06"
LOWER_FLOAT = math.nextafter(float(LOWER), -math.inf)
UPPER_FLOAT = math.nextafter(float(UPPER), math.inf)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_original():
    location = Path(__file__).resolve()
    original = location.with_name("check_three_cell_finite_amplitude.py")
    if not original.is_file():
        original = location.parents[1] / "research/plurality-stability" / original.name
    assert digest(original) == ORIGINAL_SHA256
    content = original.read_text()
    # Retarget angular integration and widen its trial-norm denominator.
    replacements = (
        ("rho = I(math.nextafter(.44, -math.inf), math.nextafter(.442, math.inf))",
         f"rho = I({LOWER_FLOAT!r}, {UPPER_FLOAT!r})"),
        ("norm = I.arb(arb(271)/1000)", "norm = I.arb(arb(272)/1000)"),
    )
    for old, new in replacements:
        assert content.count(old) == 1
        content = content.replace(old, new)
    namespace = {"__file__": str(original), "__name__": "_finite_amplitude_original"}
    exec(compile(content, str(original), "exec"), namespace)
    return namespace


def projected_assembly(original):
    """FA8-FA9, full-vector clipping AP3, and final quadratic determinant."""
    arb = original["arb"]
    pi = original["pi"]
    j1 = original["j1"]
    cs, ls = original["cs"], original["ls"]
    mb = original["M"]
    base = original["s"].base
    rat = lambda n, d=1: arb(n) / d
    rho = base.ball(LOWER_FLOAT, UPPER_FLOAT)
    d = 1-rho*rho
    alpha = rho/d.sqrt()
    kappa = rho*rho/d
    x0 = rat(1, 10)
    q = alpha*x0
    special = lambda z: pi.sqrt()*(1+2*z*z)*(z*z).exp()*z.erfc()-2*z
    weighted_j = j1*alpha*(-(1+kappa)*x0*x0).exp()*(special(q/2)+special(q)/2)
    slope = sum((c*l for c, l in zip(cs, ls) if c > 0), arb(0))
    pole = (sum(cs, arb(0))-slope*x0*x0)/rat(272, 1000).sqrt()
    constant = abs(mb)/rat(268, 1000).sqrt()
    A0, c, C, j = rat(1, 100), rat(15, 8), rat(14, 25), rat(107, 1000)
    gap, e, R2, P = rat(27, 20), rat(3, 500), rat(1, 15000), rat(97, 20)
    H2 = ((R2.sqrt()+(R2+gap*e).sqrt())/gap)**2
    radial = rat(2, 3).sqrt()+A0*C
    X = A0*c/radial
    Phi = 2*j*c**3*A0/(3*radial)
    zeta = rat(1, 6)-rat(2, 3).sqrt()*A0*P
    anchor = zeta*zeta/20-(2*P*A0)**2/240
    rest = rat(13, 20)**2-(rat(22, 5)-rat(399, 200)**2)
    final = (rat(1, 160)-e)*(gap-rat(20, 27))-R2
    assert weighted_j > j and pole > c and constant < C
    assert X < x0 and Phi > H2 and zeta > 0 and anchor > 0
    assert j1 > rat(9, 100) and rest > 0 and final > 0
    values = dict(weighted_J_local=weighted_j, pole_lower=pole,
                  pole_constant_upper=constant, integration_endpoint=X,
                  feasibility_Phi=Phi, feasibility_H_squared=H2,
                  feasibility_margin=Phi-H2, outside_cube_margin=anchor,
                  orthogonal_block_margin=rest, final_determinant=final)
    return {name: list(base.bounds(value)) for name, value in values.items()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    assert Q.from_float(LOWER_FLOAT) < LOWER < UPPER < Q.from_float(UPPER_FLOAT)
    start = time.monotonic()
    original = load_original()
    core = original["s"]
    base = core.base
    positive = core.positive
    constants = projected_assembly(original)
    print("Analytic assembly", constants, flush=True)

    print("Rebuilding validated one-variable tables", flush=True)
    table = core.InterpolatedTable(4., 8192)
    exp_table = base.build_exp_table(64., 131072, 2.)
    controls = original["numerical_controls"](table)
    gamma = original["gamma_cover"]()
    assert gamma["status"] == "pass"
    assert gamma["attempts"] == 2*(gamma["accepted"]+gamma["skipped"])-1
    angular_band = original["angular_band"](exp_table)
    print("Angular band", angular_band, flush=True)

    mesh = original["Mesh"](16384)
    parents = []
    all_children = []
    edges = np.linspace(LOWER_FLOAT, UPPER_FLOAT, 3)
    for lo, hi in zip(edges[:-1], edges[1:]):
        lo, hi = float(lo), float(hi)
        parent = original["cover"](lo, hi, mesh, table, exp_table, 38)
        parent["rho"] = [lo, hi]
        assert original["interval_covers"](parent["slabs"], (Q.from_float(lo), Q.from_float(hi)))
        for child in parent["slabs"]:
            assert lo <= child["rho"][0] < child["rho"][1] <= hi
            lo_d, hi_d = map(Q.from_float, child["D"])
            assert Q(268, 1000) <= lo_d <= hi_d <= Q(272, 1000), child
            lo_a, hi_a = map(Q.from_float, child["a"])
            assert Q(399, 200) <= lo_a <= hi_a <= Q(1003, 500), child
            assert Q.from_float(child["r2"]) <= Q(1, 15000), child
        parents.append(parent)
        all_children.extend(parent["slabs"])
        print("Spectral parent", {k: v for k, v in parent.items() if k != "slabs"}, flush=True)
    assert original["interval_covers"](all_children, (LOWER, UPPER))
    assert not original["interval_covers"](all_children[1:], (LOWER, UPPER))

    qt = positive.build_positive_q_table(4., 1024, 64, 2.)
    at = positive.build_positive_a_table(64., 4096, 64, 2.)
    upstream_constants = positive.build_constants()
    hs = []
    edges = np.linspace(LOWER_FLOAT, UPPER_FLOAT, 5)
    for lo, hi in zip(edges[:-1], edges[1:]):
        lo, hi = float(lo), float(hi)
        value = positive.certify_slab(lo, hi, 7., 6400, 48, at, qt, exp_table, upstream_constants)
        assert Q.from_float(value.total) <= Q(22, 5), value
        hs.append(dict(rho=[lo, hi], hs_squared=value.total, finite=value.finite,
                       tail=value.tail, kernel_tail=value.kernel_tail,
                       additive_tail=value.additive_tail))
        print("Hilbert--Schmidt parent", hs[-1], flush=True)
    assert original["interval_covers"](hs, (LOWER, UPPER))
    assert not original["interval_covers"](hs[1:], (LOWER, UPPER))
    controls.update(exact_rational_endpoint_cover=True, child_parent_containment=True,
                    incomplete_spectral_cover_rejected=True,
                    incomplete_hilbert_schmidt_cover_rejected=True,
                    all_spectral_caps=True, all_hilbert_schmidt_caps=True,
                    gamma_full_bisection_cover=True)
    sources = {src: digest(original["SOURCE_DIR"] / src)
               for src in ("plur_cert.py", "plur_cert_positive.py",
                           "check_positive_spectral_extension.py")}
    summary = dict(minimum_trial_norm=min(v["minD"] for v in parents),
                   maximum_trial_norm=max(v["maxD"] for v in parents),
                   minimum_rayleigh=min(v["mina"] for v in parents),
                   maximum_rayleigh=max(v["maxa"] for v in parents),
                   maximum_raw_residual_squared=max(v["maxr2"] for v in parents),
                   maximum_hilbert_schmidt_squared=max(v["hs_squared"] for v in hs))
    result = dict(status="pass", exact_correlation_range=[str(LOWER), str(UPPER)],
                  original_checker_sha256=ORIGINAL_SHA256,
                  extension_checker_sha256=digest(Path(__file__)),
                  spectral_cover=parents, hilbert_schmidt_cover=hs,
                  angular_cube=gamma, angular_band=angular_band,
                  analytic_assembly=constants, controls=controls,
                  input_sources=sources, summary=summary,
                  environment=dict(python=sys.version.split()[0],
                                   python_flint=__import__("flint").__version__,
                                   numpy=np.__version__),
                  elapsed_seconds=time.monotonic()-start)
    args.output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(dict(status="pass", elapsed_seconds=result["elapsed_seconds"],
                          summary=summary)), flush=True)


if __name__ == "__main__":
    main()
