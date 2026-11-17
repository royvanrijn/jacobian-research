"""Exact finite controls; all-degree transfers are proved in the companion note."""
from pathlib import Path
from math import prod
import hashlib
import json
import sympy as s

root = Path(__file__).resolve().parent
receipts = []
for sizes in [(1, 2), (1, 2, 3), (2, 2, 2)]:
    r = len(sizes)
    alpha = [s.Rational(k, 2) for k in sizes]
    A = sum(alpha)
    t, y = s.symbols("t:" + str(r)), s.symbols("y:" + str(r))
    X = s.symbols("x:" + str(sum(sizes)))
    R = sum(y)
    p = t[0] + s.I * t[1] ** 2 - (s.Rational(1, 3) if r == 2 else t[2] / 3)
    # Homogenize term by term so no uncancelled rational expression is used.
    F = s.expand(sum(c * prod(v ** a for v, a in zip(y, exps))
                     * R ** (2 - sum(exps))
                     for exps, c in s.Poly(p, *t).terms()))
    blocks, j = [], 0
    for k in sizes:
        blocks.append(X[j:j+k])
        j += k
    Fx = s.expand(F.subs({v: sum(z*z for z in block) / 2
                         for v, block in zip(y, blocks)}))
    markers = [(0,) * r, (1,) + (0,) * (r-1), (0, 2) + (0,) * (r-2)]
    if r == 3:
        markers.append((1, 0, 1))

    def expectation(poly, variables, kind):
        ans = 0
        for exps, coefficient in s.Poly(poly, *variables).terms():
            factor = prod(s.rf(a, e) for a, e in zip(alpha, exps))
            if kind == "dir":
                factor /= s.rf(A, sum(exps))
            ans += coefficient * factor
        return s.expand(ans)

    checks = 0
    for m in range(4):
        gaussian_terms = s.Poly(s.expand(Fx ** m), *X).terms()
        for b in markers:
            q = prod(v ** a for v, a in zip(t, b))
            G = prod(v ** a for v, a in zip(y, b))
            e = sum(b)
            dm = expectation(s.expand(q * p ** m), t, "dir")
            gm = expectation(s.expand(G * F ** m), y, "gamma")
            assert s.expand(gm - s.rf(A, e + 2*m) * dm) == 0

            # Original Gaussian monomial, not its radial average.
            qexps = []
            for block, bi in zip(blocks, b):
                qexps.extend([2*bi] + [0] * (len(block)-1))
            wick = 0
            for exps, coefficient in gaussian_terms:
                ee = [u+v for u, v in zip(exps, qexps)]
                if any(a % 2 for a in ee):
                    continue
                wick += coefficient * prod(s.factorial2(a-1) for a in ee)
            kappa = prod(2**bi * s.rf(s.Rational(1, 2), bi) / s.rf(ai, bi)
                         for ai, bi in zip(alpha, b))
            assert s.expand(wick - kappa * s.rf(A, e+2*m) * dm) == 0
            checks += 2
    receipts.append({"block_sizes": sizes, "dirichlet_parameters": [str(a) for a in alpha],
                     "exact_independent_identities": checks})

note = root / "GAUSSIAN_GAMMA_DIRICHLET_AND_OUTWARD_TRANSFERS_2026_11_16.tex"
receipt = {
    "date": "2026-11-16",
    "status": "exact_controls_pass",
    "scope": "Finite exact Dirichlet/Gamma monomial formulas compared independently to expanded Gaussian Wick moments; universal result rests on written proof.",
    "controls": receipts,
    "total_exact_identities": sum(x["exact_independent_identities"] for x in receipts),
    "proof_sha256": hashlib.sha256(note.read_bytes()).hexdigest(),
    "checker_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    "registry_promoted": False,
    "new_agent_messages": 0,
    "novelty_claim": False,
    "formal_verification": False,
    "external_review": False
}
out = root / "gaussian-gamma-dirichlet-outward-controls-2026-11-16.json"
out.write_text(json.dumps(receipt, indent=2) + "\n")
print(json.dumps(receipt, indent=2))
