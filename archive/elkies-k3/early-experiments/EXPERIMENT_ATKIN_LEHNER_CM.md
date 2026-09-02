# Atkin-Lehner CM orbit test for X_0^6(79)

Starting from the exact Eichler order O of quaternion discriminant 6 and level 79:

1. Search exact normalizer elements of reduced norm ±2, ±3, ±79.
2. Verify conjugation preserves O with unimodular action.
3. Construct the exact Gross lattice O^T.
4. Compute induced integral isometries W_2, W_3, W_79 on O^T.
5. Take primitive norm-3 and norm-24 vectors.
6. Compute their full and ±-identified Atkin-Lehner orbits and stabilizers.

Run:

    sage elkies-k3/scripts/atkin_lehner_cm_orbits.sage \
      --search-bound 12 --hit-bound 80 --targets 3,24

Key output:

    AL|normalizer|p=...
    AL|gross_action|p=...
    ALORBIT|target=3|full_size=...|pm_size=...
    ALFIX|target=3|p=2|relation=...
    ...
    ALORBIT|target=24|...

If a normalizer is missing, increase only --search-bound. Do not increase
--hit-bound unless a target seed is missing.

Caveat: this computes the action for one concrete Eichler order. Relating the
resulting embedding orbit exactly to Elkies's coordinate t on the Shimura curve
may still require matching his chosen Atkin-Lehner quotient, but this is the
correct exact arithmetic invariant to compute next.
