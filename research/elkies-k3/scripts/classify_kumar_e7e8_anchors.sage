from sage.all import *
from sage.quadratic_forms.genera.genus import Genus
from pathlib import Path


BASE = Path(__file__).resolve().parents[1]


def load_gram(path):
    return matrix(
        ZZ,
        [
            [ZZ(value) for value in line.split()]
            for line in path.read_text().splitlines()
            if line.strip() and not line.startswith("#")
        ],
    )


target = load_gram(BASE / "data/lattice/rank17_gram.txt")
assert target.nrows() == 17 and target.det() == 948

E7 = CartanMatrix(["E", 7])
E8 = CartanMatrix(["E", 8])
roots = block_diagonal_matrix(E7, E8)

# The nonzero E7 discriminant class is represented by the minuscule
# fundamental weight.  In this numbering it is the last fundamental weight,
# of norm 3/2.  A section meeting the nonidentity E7 component has this root
# part; a section meeting the identity component has root part zero.
weight_pairing = vector(ZZ, [0, 0, 0, 0, 0, 0, 1] + [0] * 8)
weight = vector(QQ, [0, 0, 0, 0, 0, 0, 1]) * E7.inverse()
assert weight * E7 * weight == QQ(3) / 2


def kumar_frame(a, b, c, epsilon):
    """Glue height Gram (1/2)[a,b;b,c] to E7+E8."""
    e1, e2 = epsilon
    height = matrix(
        QQ,
        [[QQ(a) / 2, QQ(b) / 2], [QQ(b) / 2, QQ(c) / 2]],
    )
    component = matrix(
        ZZ,
        15,
        2,
        lambda i, j: weight_pairing[i] * epsilon[j],
    )
    section = height + QQ(3) / 2 * matrix(
        ZZ,
        [[e1 * e1, e1 * e2], [e1 * e2, e2 * e2]],
    )
    assert all(value.denominator() == 1 for value in section.list())
    frame = block_matrix(
        [[roots, component], [component.transpose(), matrix(ZZ, section)]]
    )
    assert frame.det() == 948
    assert all(value % 2 == 0 for value in frame.diagonal())
    return height, frame


def root_invariants(frame):
    form = QuadraticForm(ZZ, frame)
    short = form.short_vector_list_up_to_length(2, up_to_sign_flag=True)
    half_roots = short[1] if len(short) > 1 else []
    if not half_roots:
        return 0, 0, 1
    root_module = matrix(ZZ, [list(row) for row in half_roots]).row_module()
    basis = root_module.basis_matrix()
    return basis.rank(), 2 * len(half_roots), abs(
        (basis * frame * basis.transpose()).det()
    )


def smith_factors(frame):
    diagonal = frame.smith_form()[0]
    return tuple(
        abs(ZZ(diagonal[i, i]))
        for i in range(diagonal.nrows())
        if abs(ZZ(diagonal[i, i])) != 1
    )


# If H is a rank-2 height Gram of determinant 474, then G=2H is integral
# with determinant 4*474=1896.  Every positive binary form has a reduced
# representative 0 <= 2b <= a <= c, and then a <= sqrt(4*1896/3).
scaled_determinant = ZZ(4 * 474)
reduced_forms = []
bound = ZZ(floor(sqrt(QQ(4) * scaled_determinant / 3)))
for a in range(1, bound + 1):
    for b in range(0, a // 2 + 1):
        if (scaled_determinant + b * b) % a:
            continue
        c = (scaled_determinant + b * b) // a
        if 2 * b <= a <= c:
            reduced_forms.append((ZZ(a), ZZ(b), ZZ(c)))

# Integrality and evenness of the glued frame determine epsilon from each
# diagonal modulo 4 and impose the corresponding off-diagonal parity.
admissible = []
for a, b, c in reduced_forms:
    for e1 in (0, 1):
        for e2 in (0, 1):
            if a % 4 != (1 if e1 else 0):
                continue
            if c % 4 != (1 if e2 else 0):
                continue
            if (b + 3 * e1 * e2) % 2:
                continue
            height, frame = kumar_frame(a, b, c, (e1, e2))
            admissible.append(((a, b, c), (e1, e2), height, frame))

assert len(reduced_forms) == 36
assert len(admissible) == 12

genus_hits = []
for scaled, epsilon, height, frame in admissible:
    if Genus(frame) == Genus(target):
        genus_hits.append((scaled, epsilon, height, frame))

expected = [
    ((5, 2, 380), (1, 0)),
    ((8, 0, 237), (0, 1)),
    ((21, 6, 92), (1, 0)),
]
assert [(tuple(scaled), epsilon) for scaled, epsilon, _, _ in genus_hits] == expected

pinned = [
    BASE / "data/fibrations/kumar_e7e8_mw2_frame_1.txt",
    BASE / "data/fibrations/kumar_e7e8_mw2_frame_2.txt",
    BASE / "data/fibrations/kumar_e7e8_mw2_frame_3.txt",
]

print(
    "KUMAR17|stage=enumeration|reduced_binary_forms={}|admissible_even_glues={}".format(
        len(reduced_forms), len(admissible)
    ),
    flush=True,
)
for index, ((scaled, epsilon, height, frame), path) in enumerate(
    zip(genus_hits, pinned), 1
):
    stored = load_gram(path)
    assert stored == frame
    assert height.det() == 474
    assert smith_factors(frame) == (948,)
    root_rank, root_count, root_determinant = root_invariants(frame)
    assert (root_rank, root_count, root_determinant) == (15, 366, 2)
    # Clear denominators without changing the integral automorphism group.
    height_form = QuadraticForm(ZZ, matrix(ZZ, 4 * height))
    automorphisms = height_form.number_of_automorphisms()
    assert automorphisms == (4 if index == 2 else 2)
    print(
        "KUMAR17|anchor={}|scaled_height={}|epsilon={}|height={}|"
        "height_aut_order={}|root_rank={}|roots={}|root_det={}|"
        "disc_group=Z/948|genus_match=1".format(
            index,
            scaled,
            epsilon,
            height,
            automorphisms,
            root_rank,
            root_count,
            root_determinant,
        ),
        flush=True,
    )

# U plus the negative of any hit has signature (1,18), cyclic discriminant
# group of length one, and the same local genus as the recovered target.
# Nikulin's rank-versus-length uniqueness theorem therefore identifies their
# stable Neron--Severi lattices.  This is an existence result; the explicit
# U-embedding/neighbor transport is the next calculation.
print(
    "KUMAR17|stable_NS|signature=1,18|rank=19|disc_length=1|"
    "nikulin_unique_genus=1|same_integral_NS_class=1",
    flush=True,
)
print("KUMAR17|status=PASS", flush=True)
