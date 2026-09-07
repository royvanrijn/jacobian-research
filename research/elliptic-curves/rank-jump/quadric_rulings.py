#!/usr/bin/env python3
"""Exact shared-quadric and ruling-line analysis; no parameter/point search."""
import argparse
from itertools import permutations
from pathlib import Path
import retrospective as r
from cubic_bridge import Cubic
from cover_experiment import mul, sub, evaluate, sqrtq

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE / "QUADRIC_RULING_PROTOCOL.json"
ANCHOR = r.OUT / "rank_jump_local_collision_inputs_v1.json"
BRIDGE = r.OUT / "rank_jump_cubic_bridge_v1.json"
OUTPUT = r.OUT / "rank_jump_quadric_rulings_v1.json"
CHECKPOINT = r.ROOT / "artifacts/local/rank-jump-quadric-rulings-v1"


class Quad:
    """Exact a+b sqrt(B); B is set once for the frozen anchor."""
    B = None

    def __init__(self, a=0, b=0):
        if isinstance(a, Quad):
            self.a, self.b = a.a, a.b
        else:
            self.a, self.b = r.F(a), r.F(b)

    def __add__(self, other):
        other = Quad(other)
        return Quad(self.a+other.a, self.b+other.b)

    __radd__ = __add__

    def __neg__(self):
        return Quad(-self.a, -self.b)

    def __sub__(self, other):
        return self + -Quad(other)

    def __rsub__(self, other):
        return Quad(other) + -self

    def __mul__(self, other):
        other = Quad(other)
        return Quad(self.a*other.a+self.B*self.b*other.b,
                    self.a*other.b+self.b*other.a)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = Quad(other)
        d = other.a*other.a-self.B*other.b*other.b
        assert d
        return self * Quad(other.a/d, -other.b/d)

    def __rtruediv__(self, other):
        return Quad(other)/self

    def __pow__(self, n):
        if n < 0:
            return (Quad(1)/self)**(-n)
        out = Quad(1)
        for _ in range(n):
            out *= self
        return out

    def __eq__(self, other):
        other = Quad(other)
        return self.a == other.a and self.b == other.b

    def __bool__(self):
        return bool(self.a or self.b)

    def encode(self):
        return [str(self.a), str(self.b)]


def determinant(M):
    n = len(M)
    return sum((-1)**sum(p[i] > p[j] for i in range(n) for j in range(i+1,n))
               * product(M[i][p[i]] for i in range(n))
               for p in permutations(range(n)))


def product(values):
    result = 1
    for value in values:
        result *= value
    return result


def trace(K, z):
    return 3*z[0]-2*K.A*z[2]


def scale(poly, a):
    return [a*x for x in poly]


def add(a,b):
    return sub(a, scale(b,-1))


def form_matrix(K, beta, coordinate):
    e = [tuple(r.F(i == j) for i in range(3)) for j in range(3)]
    return [[K.mul(beta, K.mul(x,y))[coordinate] for y in e] for x in e]


def quadric_isomorphism(K,p,q):
    A,B = K.A,K.B
    H = p*p+A
    T = [[H,-H*p,-H*A,H], [1,-p,-A,-1], [0,H,B,0], [0,0,q,0]]
    standard = [[0,r.F(1,2),0,0], [r.F(1,2),0,0,0],
                [0,0,-1,0], [0,0,0,B]]
    source = [row+[r.F(0)] for row in form_matrix(K,(p,-1,0),1)]
    source.append([0,0,0,1])
    assert determinant(T)
    assert all(sum(T[a][i]*standard[a][b]*T[b][j] for a in range(4) for b in range(4))
               == -H*source[i][j] for i in range(4) for j in range(4))
    return [[str(x) for x in row] for row in T]


def section(K, p, q, sign):
    A,B = K.A,K.B
    h = Quad(0,1)
    b = (-B+sign*q*h)/(p*p+A)
    v = (p*b+A, b, Quad(1))
    beta = (p,r.F(-1),r.F(0))
    bv = K.mul(beta,v)
    bvv = K.mul(beta,K.square(v))
    assert bv[1] == 0 and bvv[1] == 0
    u = [Quad(0), -2*bv[2], -bvv[2]]
    assert u[1] and u[2]
    x = [Quad(p), 2*bv[0], bvv[0]]
    y = scale([Quad(1), trace(K,v),
               (trace(K,v)**2-trace(K,K.square(v)))/2, K.norm(v)], q)
    # Coefficientwise curve identity in F[s].
    u2, u3, x2 = mul(u,u), mul(mul(u,u),u), mul(x,x)
    rhs = add(add(mul(x2,x),scale(mul(u,x2),2*A)),
              add(mul(add(add([A],scale(u,3*B)),scale(u2,A*A)),x),
                  add(add([B],scale(u2,A*B)),scale(u3,-B*B))))
    assert sub(mul(y,y),rhs) == [0]
    branch = -u[1]**2/(4*u[2])
    translated_x = (q+sign*h)**2/(p*p)-p  # x(P-sign*T), T=(0,sqrt(B)).
    assert branch == 1/translated_x
    assert branch.b and 1+A*branch**2+B*branch**3
    other_s = -u[1]/u[2]
    other = (evaluate(x,other_s),evaluate(y,other_s))
    assert other != (p,q)
    branch_s = -u[1]/(2*u[2])
    assert evaluate(u,branch_s) == branch
    return {
        "sign": sign, "line_direction_in_cubic_basis": [Quad(a).encode() for a in v],
        "u_polynomial_ascending": [Quad(a).encode() for a in u],
        "x_polynomial_ascending": [Quad(a).encode() for a in x],
        "y_polynomial_ascending": [Quad(a).encode() for a in y],
        "finite_branch_value": branch.encode(),
        "translated_anchor_x": translated_x.encode(),
        "other_preimage_of_zero": other_s.encode(),
        "other_specialized_point": [Quad(a).encode() for a in other],
        "point_identity": True, "non_invariant_section": True,
        "finite_branch_is_smooth": True
    }


def bindings():
    paths = [PROTOCOL, ANCHOR, BRIDGE, Path(__file__),
             HERE/"cubic_bridge.py", HERE/"cover_experiment.py", HERE/"retrospective.py"]
    return {str(p.relative_to(r.ROOT)): r.digest(p.read_bytes()) for p in paths}


def build(check=False):
    anchor, bridge = r.read(ANCHOR), r.read(BRIDGE)
    A,B = map(r.F, anchor["anchor"]["short_model_ainvariants"][3:])
    assert B and sqrtq(B) is None and A**3+8*B**2 != 0
    K = Cubic(A,B)
    Quad.B = B
    determinants = []
    for w in bridge["norm_witnesses"]:
        beta = tuple(map(r.F,w["beta"]))
        M = form_matrix(K,beta,1)
        det = determinant(M)
        norm = K.norm(beta)
        assert det == B*norm and sqrtq(norm) is not None
        determinants.append({"mask": w["mask"], "theta_coefficient_matrix": [[str(x) for x in row] for row in M],
                             "determinant": str(det), "norm_beta": str(norm),
                             "quadric_discriminant_squareclass": "B"})
    rows = []
    for i, point in enumerate(anchor["anchor"]["known_points_on_short_model"]):
        p,q = map(r.F, point)
        assert p and q and p*p+A and q*q == p**3+A*p+B
        row = {"basis_index": i, "anchor_point": point,
               "isomorphism_to_XY_equals_Z2_minus_BW2": quadric_isomorphism(K,p,q),
               "rulings": [section(K,p,q,sign) for sign in (-1,1)]}
        rows.append(row)
        if not check:
            CHECKPOINT.mkdir(parents=True, exist_ok=True)
            path = CHECKPOINT/f"case-{i}.json"
            if path.exists():
                assert r.read(path) == row
            else:
                r.write_new(path,row)
    branches = [tuple(z["finite_branch_value"]) for row in rows for z in row["rulings"]]
    assert len(rows) == 20 and len(set(branches)) == 40
    out = {
        "schema": "rank-jump.quadric-rulings.v1", "bindings": bindings(),
        "anchor_A": str(A), "anchor_B": str(B), "quadric_determinants": determinants,
        "ruling_constructions": rows,
        "summary": {"shared_ruling_field": "Q(sqrt(B))", "distinct_finite_branch_values": 40,
                    "common_branch_value": "infinity", "each_base_change_geometric_rank": 2,
                    "each_base_change_rank_over_Q_sqrt_B_s": 2,
                    "same_quadratic_base_change_collisions": 0,
                    "simultaneous_twenty_selected_rulings_genus": 1+2**18*17},
        "boundary": "Retrospective lines use anchor points. Rank two is over Q(sqrt(B))(s), not over Q(u) or a new specialized Q-curve. Shared quadrics do not imply solubility of their individual fibres."
    }
    if check:
        assert r.read(OUTPUT) == out
        print("PASS exact quadric and ruling replay")
    else:
        r.write_new(OUTPUT,out)
    print(out["summary"])


def verify():
    """Independent symbolic Sage identities and all stored section identities."""
    from sage.all import QQ, PolynomialRing, matrix, QuadraticField
    ring = PolynomialRing(QQ, names=("A","B","p","q","h"))
    A,B,p,q,h = ring.gens()
    M = matrix(ring,[[-1,p,A],[p,A,B-p*A],[A,B-p*A,-p*B-A*A]])
    assert M.det() == B*(p**3+A*p+B)
    # Verify branch identity universally modulo the two defining equations.
    fraction = ring.fraction_field()
    b = fraction((-B+q*h)/(p*p+A))
    Q2 = -p*b*b+2*p*p*b+p*A+B
    identity = (b-p)**2 * ((q+h)**2-p**3) - p*p*Q2
    ideal = ring.ideal([q*q-p**3-A*p-B,h*h-B])
    assert ideal.reduce(ring(identity.numerator())) == 0
    data = r.read(OUTPUT)
    assert data["bindings"] == bindings()
    field = QuadraticField(QQ(data["anchor_B"]), "h")
    h = field.gen()
    polys = PolynomialRing(field,"s")
    A,B = QQ(data["anchor_A"]),QQ(data["anchor_B"])
    def decode(pair):
        return QQ(pair[0])+QQ(pair[1])*h
    for row in data["ruling_constructions"]:
        for z in row["rulings"]:
            u,x,y = [polys([decode(v) for v in z[name+"_polynomial_ascending"]])
                     for name in ("u","x","y")]
            assert y*y == x**3+2*A*u*x*x+(A+3*B*u+A*A*u*u)*x+B+A*B*u*u-B*B*u**3
            assert -u[1]**2/(4*u[2]) == decode(z["finite_branch_value"])
    print("PASS Sage symbolic determinant, branch identity and forty section identities")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("mode", choices=("build","check","verify"))
    mode = parser.parse_args().mode
    verify() if mode == "verify" else build(mode == "check")
