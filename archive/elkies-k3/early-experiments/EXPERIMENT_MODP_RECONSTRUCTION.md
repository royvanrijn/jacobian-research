# Elkies rank-17 K3: finite-field reconstruction probe

This attack uses the existing three-hub relation graph but eliminates `A(t)` and `B(t)` from the polynomial system. The remaining unknowns are 15 quartic x-coordinates, one sextic root y-coordinate, and 21 quadratic edge slopes.

## Why this experiment

The current 167-variable system is too large to attack over Q directly. The intended route is:

1. reduce the system;
2. work modulo several good primes;
3. determine the actual component dimension after gauge freedom;
4. obtain a zero-dimensional generic slice;
5. solve one mod-p slice;
6. Hensel/Newton lift a smooth solution;
7. rationally reconstruct coefficients or, more likely, identify the 1-dimensional Shimura/moduli parameter.

The expected coordinate-space dimension is roughly 5: three base PGL2 parameters, one Weierstrass scaling, and one genuine K3 moduli parameter. This is a hypothesis to test, not an assumption to bake into a proof.

## Install (Linux compute box)

    sudo apt update && sudo apt install -y msolve

Use your existing Sage installation.

## Generate one system

    sage elkies-k3/scripts/export_three_hub_msolve.sage \
      --p 101 --slices 5 --seed 1 \
      --out artifacts/local/elkies-k3/p101-s5-seed1.ms

Dimension/GB probe:

    msolve -t 32 -v 2 -g 1 \
      -f artifacts/local/elkies-k3/p101-s5-seed1.ms \
      -o artifacts/local/elkies-k3/p101-s5-seed1.dim

If this is zero-dimensional, rerun without `-g 1` to ask msolve for a finite-field parametrization:

    msolve -t 32 -v 2 \
      -f artifacts/local/elkies-k3/p101-s5-seed1.ms \
      -o artifacts/local/elkies-k3/p101-s5-seed1.solve

## Recommended first batch

    python3 elkies-k3/scripts/run_three_hub_msolve_batch.py \
      --primes 101,103,107,109 \
      --slices 4,5,6 \
      --seeds 1,2,3 \
      --jobs 2 --threads 16 --timeout 7200

On a 32-core box, `--jobs 2 --threads 16` is the sensible first run.

## What I want back

The most useful outputs are the lines beginning `PROBE|` plus, for any `GB_DONE` case, the corresponding `.dim`/`.log` files. Do not increase to huge primes yet. Small good primes make the algebra much cheaper and are exactly what we want for model recovery.
