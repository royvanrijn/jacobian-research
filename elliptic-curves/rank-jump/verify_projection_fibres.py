#!/usr/bin/env sage -python
"""Exact replay of projection fibres; no conic solving, factoring or witness search."""
import json
import sys
from pathlib import Path
import projection_fibres as pf
import retrospective as r
import local_collision as lc
from sage.all import QQ, PolynomialRing, EllipticCurve, matrix, vector, pari, prod

sys.path.insert(0, str(pf.old.CAS))
from research_runtime.arithmetic import ArithmeticContext
from research_runtime.sage_subspace import SageSubspaceBackend
from research_runtime.fisher import invariants


class RetainedBackend(SageSubspaceBackend):
    """Only exact polynomial and retained-pair operations are available."""
    def __init__(self, context):
        self.context = context
        self.R = PolynomialRing(QQ, "x")
        self.K = self.R.quotient(self.R(list(map(QQ, context.two_torsion.polynomial))), "th")
        self.theta = self.K.gen()
        self.alpha = self.K(list(map(QQ, context.curve_generator_in_algebra)))
        self.cubic = self.R(list(map(QQ, context.minimal_model.two_division_polynomial)))
        self.a2, self.a4, self.a6 = self.cubic[2], self.cubic[1], self.cubic[0]
        self.I = self.a2**2-3*self.a4
        self.J = -2*self.a2**3+9*self.a2*self.a4-27*self.a6
        self.arithmetic = None  # Any accidental discovery call must fail.


def verify_cover(row, backend, protocol, E, basis):
    assert row["bindings"] == pf.bindings()
    index, mask = row["index"], row["mask"]
    anchor_mask = protocol["cover_anchor_masks"][index]
    assert row["anchor_mask"] == anchor_mask
    assert lc.lift(mask, protocol["anchor_basis_masks"]) == anchor_mask
    point = sum((P for i,P in enumerate(basis) if anchor_mask>>i&1), E(0))
    p, q = point[0:2]
    assert row["anchor_point"] == [str(p), str(q)]
    A, B = E.a4(), E.a6()
    assert pf.verify_geometry(A,B,p,q) == row["geometry"]
    assert row["context"] == json.loads(json.dumps(backend.context.record()))
    representatives = []
    for m in protocol["anchor_basis_masks"]:
        P = sum((P for i,P in enumerate(basis) if m>>i&1), E(0))
        representatives.append(P[0]-backend.theta)
    assert row["classes"] == [list(backend.coordinates(b)) for b in representatives]
    beta = prod(b for i,b in enumerate(representatives) if mask>>i&1)
    G = matrix(QQ, row["conic_matrix"])
    assert G == matrix(QQ, [[1,0,0],[0,-p-1,-q],[0,-q,-(p*p-p+A+1)]])
    assert G.det() == 1+A-B != 0
    point_on_conic = vector(QQ, row["conic_point"])
    assert point_on_conic != 0 and point_on_conic*G*point_on_conic == 0
    parameter = matrix(QQ, row["conic_parameter_matrix"])
    assert parameter.det() != 0
    x = backend.R.gen()
    t,v,w = parameter*vector([x*x,x,1])
    assert t*t-(p+1)*v*v-2*q*v*w-(p*p-p+A+1)*w*w == 0
    raw = backend.R(list(map(QQ,row["raw_quartic"])))
    assert raw == v*v+2*t*w+(2-p)*w*w
    change = matrix(QQ, row["parameter_change"])
    assert change.det() != 0
    a,b = change*vector([x,1])
    final = backend.R(list(map(QQ,row["quartic"])))
    assert final.degree() == 4 and final.discriminant() != 0
    assert sum(raw[i]*a**i*b**(4-i) for i in range(5)) == QQ(row["raw_y_over_final_y"])**2*final
    assert invariants(final) == (backend.I,backend.J)
    phi = -3*backend.alpha-backend.a2
    cubic = (4*final[4]*phi+3*final[3]**2-8*final[4]*final[2])/3
    root = backend.K(list(map(QQ,row["cubic_invariant_over_beta_square_root"])))
    assert cubic == beta*root*root
    return final


def verify_other_quartic(A,B,points,backend):
    """A different degree-two map translates every 2-cover label by eta."""
    th = backend.theta
    u = QQ(-1)
    a2,a4,a6 = 2*A*u,A+3*B*u+A*A*u*u,B+A*B*u*u-B*B*u**3
    I,J = a2*a2-3*a4,-2*a2**3+9*a2*a4-27*a6
    alpha = th-th*th
    phi = -3*alpha-a2
    kappa = 1-th+A+th*th
    D = 1+A-B
    assert kappa*(1+th)**2 == D*(1+th)
    Eu = EllipticCurve([0,a2,0,a4,a6])
    eta_point = Eu([A+1,D])
    assert eta_point[0]-alpha == kappa
    R = backend.R
    x = R.gen()
    for p,q in (P[0:2] for P in points):
        quartic = R([p**3-3*p*p+4*B-4*A,4*q*(2-p),
                     6*p*p-6*p+4*A,-4*q,p+1])/4
        assert invariants(quartic) == (I,J)
        cubic = (4*quartic[4]*phi+3*quartic[3]**2-8*quartic[4]*quartic[2])/3
        assert cubic == (p-th)*kappa
        F = x*x-p-1
        G = 2*x*(2*p-1)-4*q
        H = -(p+1)*x*x+4*q*x-3*p*p-4*A
        assert G*G-4*F*H == 16*quartic
    return len(points)


def verify():
    data = r.read(pf.OUTPUT)
    assert data["schema"] == "rank-jump.projection-fibres.v1"
    assert data["bindings"] == pf.bindings()
    protocol = r.read(pf.PROTOCOL)
    assert protocol["parameter_u"] == -1
    source = r.read(lc.INPUT)
    A,B = map(QQ,source["anchor"]["short_model_ainvariants"][3:])
    E = EllipticCurve([A,B])
    basis = [E(list(map(QQ,P))) for P in source["anchor"]["known_points_on_short_model"]]
    assert len(data["covers"]) == 6 and len(data["pairs"]) == 3
    context = ArithmeticContext.from_record(data["covers"][0]["context"])
    assert context.model.coefficients == tuple(map(str,[0,-2*A,0,A-3*B+A*A,B+A*B+B*B]))
    backend = RetainedBackend(context)
    scale,shift,s,t = map(QQ,context.minimal_to_input)
    assert s == t == 0
    assert backend.alpha == 4*(backend.theta-backend.theta**2-shift)/scale**2
    for i,row in enumerate(data["covers"]):
        assert row["index"] == i
        verify_cover(row,backend,protocol,E,basis)
    count = verify_other_quartic(A,B,basis,backend)
    print("PASS six projection fibres, Mumford/norm maps, exact Kummer labels, and",count,"alternate labels")
    old_matrix = next(x["matrix"] for x in source["ct"] if x["u"] == -1)
    values = []
    for i,row in enumerate(data["pairs"]):
        assert row["bindings"] == pf.bindings() and row["index"] == i
        if row.get("status") == "UNKNOWN":
            print("UNKNOWN bounded pairing",i,row["reason"])
            values.append(None)
            continue
        masks = protocol["pair_class_masks"][i]
        covers = [next(c for c in data["covers"] if c["mask"] == m) for m in masks]
        result = backend._pair(None,masks,covers,retained=row["pair"])
        for term in result["local_terms"]:
            place = term["place"]
            a,b = QQ(covers[1]["quartic"][4]),QQ(term["gamma_value"])
            symbol = (-1 if a < 0 and b < 0 else 1) if place == "infinity" else int(pari.hilbert(a,b,place))
            assert symbol == term["hilbert_symbol"]
        expected = old_matrix[{1:0,2:1,4:2}[masks[0]]][{1:0,2:1,4:2}[masks[1]]]
        assert expected == row["expected_retained_CT_value"] == result["value"]
        values.append(result["value"])
        print("PASS pairing",masks,result["value"],"and",len(result["local_terms"]),"Hilbert symbols")
    return values


if __name__ == "__main__":
    verify()

