# Q80 orbit 1222 characteristic-zero pipeline

These scripts are the active/reproducible pipeline for the certified Q80
orbit-1222 characteristic-zero Jacobian.  They resolve repository paths from
their own permanent location, so the commands may be launched from any
working directory.

Permanent exact inputs and outputs are in
[`../../data/fibrations/q80-orbit1222-char0/`](../../data/fibrations/q80-orbit1222-char0/).
Cached modular packets and Frobenius-character audits are in
[`../../../artifacts/generated-results/q80-orbit1222-char0/modular/`](../../../artifacts/generated-results/q80-orbit1222-char0/modular/).
Superseded experiments and older script copies are retained under the sibling
`debug/` artifact directory.

## Fast certificate replay

From the repository root:

```bash
sage elkies-k3/data/fibrations/q80-orbit1222-char0/q80_char0_orbit1222_jacobian_normalized.sage
sage elkies-k3/data/fibrations/q80-orbit1222-char0/q80_char0_orbit1222_P1_P3_normalized.sage
sage elkies-k3/scripts/q80-orbit1222-char0/finalize_orbit1222.sage
```

The final command must end with these statuses:

```text
PASS_EXACT_SQUARECLASS_REDUCTION
PASS_NORMALIZED_EXACT_MODEL_AND_SECTIONS
PASS_ALL_FROBENIUS_CHARACTERS
PASS_EXACT_ORBIT1222
```

It checks the pinned 16 rational primes (32 split places), rewrites the
normalized exact model and self-contained section file, and refreshes
`Q80_CHAR0_ORBIT1222_FINAL_CERTIFICATE.md` without performing modular
reconstruction.

## Reconstruction stages

The cached-data stages are, in dependency order:

```bash
# Generate/reuse the four modular embeddings for one split prime.
python3 elkies-k3/scripts/q80-orbit1222-char0/modular_packet.py --prime 79

# Recover the certified exact monic C and S using all cached packets and the
# independent p=73 holdout.
sage elkies-k3/scripts/q80-orbit1222-char0/export_exact_cs.sage

# Optional fast modular mu packet for an additional split prime.
sage elkies-k3/scripts/q80-orbit1222-char0/mu_prime_fast.sage \
  --prime 1009 \
  --out artifacts/generated-results/q80-orbit1222-char0/modular/q80_mu_fast_p1009.json

# Recover exact mu from the critical factors, then build the exact base model.
sage elkies-k3/scripts/q80-orbit1222-char0/mu_from_critical_factors.sage
sage elkies-k3/scripts/q80-orbit1222-char0/build_j_fiber_model.sage

# Optional Frobenius-character packet for an additional split prime.
sage elkies-k3/scripts/q80-orbit1222-char0/twist_character_probe.sage \
  --prime 79 \
  --out artifacts/generated-results/q80-orbit1222-char0/modular/q80_twchar_p79.json

# Expensive exact section/twist lift, followed by the fast final certificate.
sage elkies-k3/scripts/q80-orbit1222-char0/lift_exact_twist_P1_P3.sage
sage elkies-k3/scripts/q80-orbit1222-char0/finalize_orbit1222.sage
```

`simultaneous_reconstruct_p73.sage` retains the successful simultaneous
`C,S,mu` reconstruction route.  `mu_reconstruct_mixed.py` is an independent
mixed-cache rational-reconstruction cross-check.  Neither is needed for the
fast final certificate replay.

The modular fiber computations use Sage, Singular's `brnoeth.lib`, and the
compatibility patch embedded in the worker scripts.  They are substantially
more expensive than the cached certificate replay.
