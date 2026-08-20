# Gross/Clifford companion test

Literature bridge:

For a quaternion order O, define the Gross lattice

    O^T = { 2x - Tr(x) : x in O }.

An optimal embedding of a quadratic order of discriminant Delta<0 into O is
equivalent to a primitive element iota(sqrt(Delta)) in O^T of reduced norm -Delta.

Our recovered rank-3 K3 transcendental lattice T is known only as an integral
quadratic lattice in the same rational ternary/quaternionic world. It is NOT safe
to assert that T itself is O^T with identical scale.

For a ternary quadratic form q with matrix A, the trace-zero norm on the even
Clifford algebra is naturally represented rationally by the adjugate/exterior-square
form adj(A). This experiment therefore compares:

- q_T;
- -q_T;
- adj(q_T);
- -adj(q_T);
- q_T^{-1} and its sign;

and their primitive integral normalizations.

It tests representations of 3 and 24, corresponding to the two CM discriminants
surviving the Eichler embedding-count fingerprint.

Run:

    sage elkies-k3/scripts/test_gross_companion.sage --targets 3,24 --bound 80

The `SQUARECLASS` lines separately test whether negative vectors in T have
rational norm square-class -3 or -6 (the CM fields of discriminants -3 and -24).

Do not yet claim that a hit identifies t=2. A hit tells us which normalization of
the Clifford/Gross construction is compatible and what exact order lattice needs
to be reconstructed next.
