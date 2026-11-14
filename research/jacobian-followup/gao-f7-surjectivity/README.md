# Gao F7 is surjective

**GAOF7SURJ1**, 14 November 2026: every complex target of Gao's explicit
five-dimensional F7 has a preimage. This closes the empty-fibre question
in the preceding F7 follow-up. The map remains noninjective of generic
degree twelve and nonproper.

An empty eliminant would have to be `Q2^3*T3^2` or `Q4^3`. The quartic-cube
equations have an explicit identity equal to one. Centering the quadratic
reduces the other case to two variables; an exact identity forces its
discriminant to vanish, contradicting the required squarefreeness.

[Written proof](PROOF.md) · [checker](verify.py) · [receipt](receipt.json) ·
[two sparse rational identities](certificates.json) · [source pins](sources/manifest.json).

The identities were generated in Singular and independently replayed by
FLINT. SymPy rederives their necessary input equations. Verification needs
no Groebner calculation or search:

```sh
timeout 60s .venv/bin/python research/jacobian-followup/gao-f7-surjectivity/verify.py
```

Only the global cardinality three remains undecided for F7; full strata,
F5, JC2 and unrestricted HC4 remain open. Independent review and novelty
confirmation are outstanding. This is a result about this specific F7,
not a claim to the first surjective Keller example.
