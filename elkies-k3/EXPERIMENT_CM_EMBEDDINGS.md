# Identify CM orders for X_0^6(79)

Uses Eichler's optimal-embedding formula for a squarefree Eichler order:

    e_{D,N}(Delta)
      = h(Delta)
        prod_{p|D} (1 - {Delta/p})
        prod_{q|N} (1 + {Delta/q})

with D=6 and N=79.

The modified Eichler symbol is 1 when p divides the conductor of the quadratic order,
otherwise it is the Kronecker symbol of the fundamental field discriminant.

The script computes ring class numbers directly by enumerating primitive reduced
positive-definite binary quadratic forms, avoiding Sage-version-specific APIs.

First run:

    python3 elkies-k3/scripts/identify_cm_embeddings.py --disc-max 20000 --embedding-count 4

Then broaden if needed:

    python3 elkies-k3/scripts/identify_cm_embeddings.py --disc-max 200000 --embedding-count 4 --top 200

Interpretation:
- e=0: no optimal embedding into the (6,79) Eichler order.
- e=4: especially interesting because Elkies reports four CM points at |t|=2, |u|=32.
- e is before any additional quotient by Atkin-Lehner involutions, so e=4 is evidence, not identification.

The next step is to connect surviving Delta values to negative-vector orbits in the
rank-3 transcendental lattice. Do NOT identify Delta with det(v^\perp); those are
different discriminants in the K3/Shimura correspondence.
