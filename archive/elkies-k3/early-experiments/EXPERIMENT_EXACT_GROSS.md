# Exact Gross lattice for the (6,79) Eichler order

Sage supports rational quaternion algebras and `order_with_level`.

For quaternion discriminant D=6 and Eichler level N=79 we construct:

    A = QuaternionAlgebra(6)
    O = A.order_with_level(6*79)

and assert

    [O_max : O] = 79.

Then compute the exact Gross lattice

    O^T = { 2x - Tr(x) : x in O }

from the order basis and the reduced norm form.

Primitive vectors beta in O^T with norm d correspond to optimal embeddings of
quadratic orders of discriminant -d.

Run:

    sage elkies-k3/scripts/construct_exact_gross_lattice.sage --targets 3,24 --bound 80

This is now an exact quaternion-order experiment; it no longer depends on guessing
the normalization relating the K3 transcendental lattice to a Clifford companion.

Key outputs:

    EXACTGROSS|ramified_primes=(2,3)
    EXACTGROSS|index_in_maximal=79
    EXACTGROSS|norm_gram_det=...
    EXACTGROSS|target=3|primitive_count_box=...
    EXACTGROSS|target=24|primitive_count_box=...
    GROSSHIT|...

If both 3 and 24 occur, the next discriminator is the action of the Atkin-Lehner
normalizer on these primitive vectors / optimal embeddings. If only one occurs,
that identifies the viable CM order immediately.
