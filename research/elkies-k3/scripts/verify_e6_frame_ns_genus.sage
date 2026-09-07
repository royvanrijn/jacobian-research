from sage.all import *
from sage.quadratic_forms.genera.genus import Genus
from pathlib import Path


def load_gram(path):
    return matrix(ZZ, [
        [ZZ(value) for value in line.split()]
        for line in Path(path).read_text().splitlines()
        if line.strip() and not line.startswith("#")
    ])


original_mw = load_gram("elkies-k3/data/lattice/rank17_gram.txt")
e6_frame = load_gram(
    "elkies-k3/data/fibrations/mw3_e6_a3a3_a1a1_frame.txt"
)
U = matrix(ZZ, [[0, 1], [1, 0]])
original_ns = block_diagonal_matrix(U, -original_mw)
e6_ns = block_diagonal_matrix(U, -e6_frame)

assert original_mw.nrows() == e6_frame.nrows() == 17
assert original_mw.det() == e6_frame.det() == 948
assert all(value % 2 == 0 for value in original_mw.diagonal())
assert all(value % 2 == 0 for value in e6_frame.diagonal())

# The positive-definite frames, and hence their hyperbolic U-extensions, have
# exactly the same local genus symbols.
assert Genus(original_mw) == Genus(e6_frame)
assert Genus(original_ns) == Genus(e6_ns)

# Both discriminant groups are cyclic of order 948, so every p-primary length
# is one.  Rank 19 is far above l(A_p)+2=3.  The exceptional equality clause at
# p=2 in Nikulin, Theorem 1.14.2, is therefore not invoked.  That theorem says
# this indefinite even genus contains a single integral isometry class.
def nontrivial_smith_factors(gram):
    smith = gram.smith_form()[0]
    return tuple(
        abs(ZZ(smith[index, index]))
        for index in range(smith.nrows())
        if abs(ZZ(smith[index, index])) != 1
    )


assert nontrivial_smith_factors(original_ns) == (948,)
assert nontrivial_smith_factors(e6_ns) == (948,)
rank = original_ns.nrows()
discriminant_length = 1
assert rank == 19 and rank > discriminant_length+2

# Signature is (1,18): U contributes (+1,-1), and the frame is positive
# definite of rank 17 before negation.
assert original_mw.is_positive_definite() and e6_frame.is_positive_definite()

print("E6NS|original_frame_det=948|e6_frame_det=948", flush=True)
print("E6NS|signature=1,18|discriminant_group=Z/948|length=1", flush=True)
print("E6NS|local_genus_match=1|rank=19|length_plus_2=3", flush=True)
print("E6NS|nikulin_1.14.2=APPLIES|indefinite_genus_classes=1", flush=True)
print("E6NS|same_integral_NS_isometry_class=1", flush=True)
print("E6NS|note=existence_certificate_not_explicit_isometry_matrix", flush=True)
print("E6NS|PASS", flush=True)
