from sage.all import ZZ, block_diagonal_matrix, gcd, matrix, vector, xgcd
from pathlib import Path
import hashlib


ROOT = Path(__file__).resolve().parents[2]
U = matrix(ZZ, [[0, 1], [1, 0]])


def load_matrix(path):
    rows = []
    for line in path.read_text().splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            rows.append([ZZ(x) for x in line.split()])
    return matrix(ZZ, rows)


def bezout_vector(pairings):
    """Reproduce the deterministic iterative xgcd used by the searches."""
    current = ZZ(0)
    coeffs = [ZZ(0)] * len(pairings)
    for i, pairing in enumerate(pairings):
        if pairing == 0:
            continue
        new_gcd, s, t = xgcd(current, ZZ(pairing))
        coeffs = [s * x for x in coeffs]
        coeffs[i] += t
        current = new_gcd
    assert abs(current) == 1
    if current == -1:
        coeffs = [-x for x in coeffs]
    return vector(ZZ, coeffs)


def reconstruct_step(parent, a, b, v):
    ns_parent = block_diagonal_matrix(U, -parent)
    f = vector(ZZ, [ZZ(a), ZZ(b)] + list(v))
    assert f * ns_parent * f == 0
    assert gcd([abs(ZZ(x)) for x in ns_parent * f]) == 1

    g = bezout_vector(list(ns_parent * f))
    assert f * ns_parent * g == 1
    g_square = ZZ(g * ns_parent * g)
    assert g_square % 2 == 0
    g0 = g - (g_square // 2) * f
    assert g0 * ns_parent * g0 == 0
    assert f * ns_parent * g0 == 1

    orthogonal = matrix(
        ZZ, [list(f * ns_parent), list(g0 * ns_parent)]
    ).right_kernel_matrix()
    child = -(orthogonal * ns_parent * orthogonal.transpose())
    transport = matrix(ZZ, [list(f), list(g0)] + orthogonal.rows())

    assert transport.det() in (-1, 1)
    assert transport * ns_parent * transport.transpose() == block_diagonal_matrix(
        U, -child
    )
    return child, transport, f, g0


DATA = ROOT / "elkies-k3" / "data"
STEPS = [
    {
        "name": "rank17_to_q90_mw7",
        "parent": DATA / "lattice" / "rank17_gram.txt",
        "child": DATA / "fibrations" / "q90_mw7_frame.txt",
        "q": 90,
        "a": 9,
        "b": 10,
        "v": (0, 0, 0, 0, 0, 0, -2, -1, 0, 0, 6, -5, 1, 0, 0, 0, 0),
    },
    {
        "name": "q90_mw7_to_q90_mw4",
        "parent": DATA / "fibrations" / "q90_mw7_frame.txt",
        "child": DATA / "fibrations" / "q90_mw4_frame.txt",
        "q": 4,
        "a": 2,
        "b": 2,
        "v": (-1, -1, -1, -3, 0, 1, 1, -1, 0, 2, 1, -1, -2, 0, 1, 0, 0),
    },
    {
        "name": "q90_mw4_to_e6_mw3",
        "parent": DATA / "fibrations" / "q90_mw4_frame.txt",
        "child": DATA / "fibrations" / "mw3_e6_a3a3_a1a1_frame.txt",
        "q": 4,
        "a": 2,
        "b": 2,
        "v": (0, 0, -2, 0, -2, 3, -2, 2, 1, -1, -1, 0, 0, -1, -1, 0, 0),
    },
]


transports = []
for index, step in enumerate(STEPS, 1):
    parent = load_matrix(step["parent"])
    expected_child = load_matrix(step["child"])
    assert parent.nrows() == 17 and parent.det() == 948
    assert vector(ZZ, step["v"]) * parent * vector(ZZ, step["v"]) == 2 * step["q"]

    child, transport, f, g0 = reconstruct_step(
        parent, step["a"], step["b"], step["v"]
    )
    assert child == expected_child
    transports.append(transport)
    print(
        "E6CHAIN|step={}|name={}|q={}|ab={},{}|child_exact=1|transport_det={}".format(
            index,
            step["name"],
            step["q"],
            step["a"],
            step["b"],
            transport.det(),
        )
    )
    print("E6CHAIN|step={}|f={}".format(index, tuple(f)))
    print("E6CHAIN|step={}|g0={}".format(index, tuple(g0)))


composite = transports[2] * transports[1] * transports[0]
original = load_matrix(STEPS[0]["parent"])
e6 = load_matrix(STEPS[-1]["child"])
ns_original = block_diagonal_matrix(U, -original)
ns_e6 = block_diagonal_matrix(U, -e6)
assert composite.det() == 1
assert composite * ns_original * composite.transpose() == ns_e6

certificate_path = DATA / "fibrations" / "e6_ns_transport_from_rank17.txt"
certificate = load_matrix(certificate_path)
assert certificate == composite
digest = hashlib.sha256(certificate_path.read_bytes()).hexdigest()

print("E6CHAIN|composite_exact=1|det={}|max_entry={}".format(
    composite.det(), max(abs(x) for x in composite.list())
))
print("E6CHAIN|fiber_in_rank17={}".format(tuple(composite.row(0))))
print("E6CHAIN|mate_in_rank17={}".format(tuple(composite.row(1))))
print("E6CHAIN|certificate={}".format(certificate_path.relative_to(ROOT)))
print("E6CHAIN|sha256={}".format(digest))
print("E6CHAIN|status=PASS")
