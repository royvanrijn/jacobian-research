# Degree 5–20 pilot scan

The unnormalized pilot evaluated 96 exact sparse systems: six fixed-`B`
Newton-edge families at each degree from 5 through 20. Results:

- 64 systems are inconsistent over `Q`;
- each has an exact RREF contradiction row `0=1`;
- ranks agree modulo `1000003`, `1000033`, and `1000037`;
- 32 systems survive, exactly the coordinate and triangular positive controls;
- all 32 survivors are birational automorphisms.

The full machine-readable records are generated at
`results/scan_2d_5_20.json` and deliberately ignored by Git because they are
generated output.

Run the stricter pass with two colliding normalized points and identity
derivative at the origin:

```bash
.venv/bin/python scripts/scan_2d_charts.py \
  --collision-normalized \
  --output results/scan_2d_5_20_collision.json
```

That stricter pass also evaluated 96 systems and rejected all 96 by exact
linear inconsistency. The modular ranks again agree at all three primes. This
is the expected outcome for these deliberately small support families.

These are exhaustive statements only for the enumerated support families.
They are not a degree-5–20 proof of the plane Jacobian conjecture.

Characteristic-zero F4 certificates for the newer normalized Stage A support
families live in `results/certificates`. Each `.input` is a complete msolve
input and the matching `.msolve` output has reduced Gröbner basis `[1]`. Recheck
all of them with `scripts/verify_certificates.py`.
