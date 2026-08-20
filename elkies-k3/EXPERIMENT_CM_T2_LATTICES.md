# Elkies rank-17: transcendental / CM endpoint search

This replaces the failed Gröbner attacks.

Facts used:

- generic K3 has NS rank 19;
- `disc(NS)=948`, hence the transcendental lattice has rank 3, signature `(2,1)`, determinant `-948`;
- the discriminant group in the recovered rank-17 MW lattice is cyclic of order 948;
- the Shimura curve is `X(6,79)`: quaternion discriminant 6, level 79. Therefore the relevant rational ternary quadratic space should have even Clifford algebra ramified at 2 and 3.

The script searches small even integral ternary forms satisfying those constraints. For each candidate it enumerates primitive negative vectors. Their orthogonal complements are positive-definite binary even lattices: candidate transcendental lattices of CM/singular K3 specializations.

Run:

    sage elkies-k3/scripts/search_cm_t2_lattices.sage \
      --coeff-bound 8 \
      --c-bound 100 \
      --top 20 \
      --vbound 12 \
      --cm-disc-max 20000

Expected runtime is minutes, not hours. If useful candidates stabilize, increase only `--vbound`, e.g. 20 or 30.

Important: matching determinant/signature/cyclic discriminant/Clifford type is a strong genus-level filter, but this script does NOT yet prove integral discriminant-form isometry with the recovered NS lattice. The next verifier should compare the finite quadratic discriminant modules exactly once a short list exists.
