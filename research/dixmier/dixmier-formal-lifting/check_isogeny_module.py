#!/usr/bin/env python3
"""Exact polynomial matrix controls for QIM1; no search for Weyl lifts."""
from pathlib import Path
import hashlib
import json
import sympy as s

HERE = Path(__file__).resolve().parent
u, v, b, h, w = s.symbols("u v b h w")
I3 = s.eye(3)
relation = s.Poly(v**2-u**3+27*b, v)


def red(expr):
    return s.Poly(s.expand(expr), v).rem(relation).as_expr()


def zero(expr, label):
    if s.cancel(red(expr)) != 0:
        raise AssertionError(label + ": " + str(expr))


def matrix_zero(M, label):
    for i in range(M.rows):
        for j in range(M.cols):
            zero(M[i,j], f"{label}[{i},{j}]")


def deriv(expr):
    return red(v*s.diff(expr, u)+s.Rational(3,2)*u**2*s.diff(expr, v))


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def degree(M):
    return int(max([s.degree(s.expand(entry), u) for entry in M if entry != 0] or [0]))


def main():
    W = s.Matrix([[0,0,-v/4],[1,0,3*u/4],[0,1,0]])
    C = s.Matrix([[0,u/2,v/2],[0,0,-u/2],[0,-1,0]])
    X = (4*W**2-u*I3)/3
    Y = (2*u*W-v*I3)/3
    Xi = (u*X-X**2)/(4*b)
    matrix_zero(W**3-3*u*W/4+v*I3/4, "monic cubic")
    matrix_zero(Y**2-X**3-b*I3, "source elliptic equation")
    matrix_zero(X*Xi-I3, "source position unit")
    matrix_zero(u*X**2-X**3-4*b*I3, "isogeny first coordinate")
    matrix_zero(v*X**3+Y*(X**3-8*b*I3), "isogeny second coordinate")
    zero(X.det()+4*b, "norm of source position")
    disc = red(s.discriminant(w**3-3*u*w/4+v/4,w))
    zero(disc-s.Rational(729,16)*b, "constant discriminant")

    def nabla(M):
        return M.applyfunc(deriv)+C*M-M*C
    matrix_zero(nabla(W)-(u*I3/2-W**2), "lifted derivation on primitive")
    matrix_zero(nabla(X)+Y, "lifted derivation on source x")
    matrix_zero(nabla(Y)+3*X**2/2, "lifted derivation on source y")
    matrix_zero(nabla(Xi)-Xi*Y*Xi, "lifted derivation on inverse x")

    I2=s.eye(2)
    V=s.Matrix([[0,u**3-27*b],[1,0]])
    C2=s.Matrix([[0,s.Rational(3,2)*u**2],[0,0]])

    def embed(expr):
        expr=red(expr)
        pp=s.Poly(expr,v)
        return pp.nth(0)*I2+pp.nth(1)*V

    def blocks(M):
        return s.BlockMatrix([[embed(M[i,j]) for j in range(M.cols)]
                              for i in range(M.rows)]).as_explicit()

    W6,X6,Y6,Xi6 = map(blocks,[W,X,Y,Xi])
    V6=s.diag(V,V,V)
    leading=V6
    lower=h*(s.diag(C2,C2,C2)+blocks(C))

    def actual_comm(M):
        matrix_zero(leading*M-M*leading, "first-order symbol compatibility")
        return h*leading*M.diff(u)+lower*M-M*lower

    matrix_zero(V**2-(u**3-27*b)*I2,"rank-two target curve")
    matrix_zero(actual_comm(u*s.eye(6))-h*V6,"matrix derivative on target u")
    matrix_zero(actual_comm(V6)-3*h*u**2*s.eye(6)/2,"matrix derivative on target v")
    matrix_zero(actual_comm(W6)-h*(u*s.eye(6)/2-W6**2),"matrix derivative on w")
    matrix_zero(actual_comm(X6)+h*Y6,"matrix derivative on x")
    matrix_zero(actual_comm(Y6)+3*h*X6**2/2,"matrix derivative on y")
    matrix_zero(actual_comm(Xi6)-h*Xi6*Y6*Xi6,"matrix derivative on inverse x")
    matrix_zero(X6*Xi6-s.eye(6),"polynomial matrix inverse x")
    matrix_zero(Y6**2-X6**3-b*s.eye(6),"rank-six source curve")
    # An actual Weyl pair inside the full curve operator algebra.
    pair_leading = -s.Rational(2,3)*Xi6**2*leading
    pair_lower = -s.Rational(2,3)*Xi6**2*lower
    matrix_zero(pair_leading*Y6-Y6*pair_leading, "canonical-pair top order")
    pair_comm = h*pair_leading*Y6.diff(u)+pair_lower*Y6-Y6*pair_lower
    matrix_zero(pair_comm-h*s.eye(6), "actual finite matrix Weyl pair")
    pair_degree = max(degree(pair_leading)+1, degree(pair_lower))
    degrees={
        "w":degree(W6),"x":degree(X6),"y":degree(Y6),
        "inverse_x":degree(Xi6),
        "h_delta":max(degree(leading)+1,degree(lower))}
    # b is an invertible scalar parameter; h carries no spatial degree.
    if max(degrees.values()) > 4:
        raise AssertionError("declared degree bound four")

    # The exact source chart also gives the boundary norm obstruction at b=0.
    xx,yy=s.symbols("xx yy")
    U=4*yy**2/xx**2-3*xx
    VV=8*yy**3/xx**3-9*yy
    zero((U**3-VV**2)-27*(yy**2-xx**3),"total-family energy identity")
    jac=s.det(s.Matrix([s.diff(U,z) for z in (xx,yy)]).row_join(
        s.Matrix([s.diff(VV,z) for z in (xx,yy)])))
    zero(jac-27,"total-family constant Jacobian")
    # For fixed u!=0 at b=0, the cubic for x has a double omitted root x=0.
    zero((xx**3-u*xx**2+4*b).subs(b,0)-xx**2*(xx-u),"boundary cubic split")

    # Bind the actual Long three-isogeny identification, without replaying it.
    old=json.loads((HERE/"quantum_elliptic_energy_controls.json").read_text())
    if old["status"] != "PASS":
        raise AssertionError("retained QEI1 status")
    for name,expected in old["hashes"].items():
        path=HERE/name
        if not path.exists():
            path=HERE.parent/"dixmier-audit"/name
        if sha(path)!=expected:
            raise AssertionError("retained QEI1 hash: "+name)

    out={"claim":"QIM1","status":"PASS","coefficient_scope":"Q[b,b^-1,h], every b!=0",
         "monic_etale_cubic":"w^3-3*u*w/4+v/4",
         "discriminant":str(disc),"norm_x":str(red(X.det())),
         "rank_three_curve_module":"explicit basis 1,w,w^2; full differential action checked",
         "rank_six_polynomial_matrix_degrees":degrees,
         "uniform_matrix_generator_degree_bound":4,
         "matrix_generator_count":5,
         "actual_finite_matrix_Weyl_pair":{"relation":"[Q,Y]=h*I6","Y_degree":degree(Y6),"Q_degree":pair_degree,"Q_formula":"-(2/3)*X_inverse^2*Z"},
         "curve_and_isogeny_relations":"PASS",
         "full_first_order_matrix_derivation_relations":"PASS",
         "total_family_Jacobian":27,"singular_energy_boundary_split":"x^2*(x-u)",
         "retained_Long_isogeny_hash_chain":"PASS",
         "hashes":{name:sha(HERE/name) for name in [
             "ELLIPTIC_ISOGENY_MATRIX_MODULE.md","check_isogeny_module.py",
             "ELLIPTIC_ISOGENY_AND_QUANTUM_ENERGY.md",
             "quantum_elliptic_energy_controls.json"]},
         "coverage_boundary":"Exact identities for the parameterized curve-module matrix construction. Freeness, faithfulness in every operator degree, the scalar genus obstruction intrinsic boundary nonfiniteness and the all-matrix-size central-energy trace obstruction use written proofs. This is not a full Weyl lift of Long; the central energy parameter and elliptic curve differential operators cannot be identified with its noncentral quantum energy.",
         "independent_replay":False,"formal_verification":False,"external_review":False}
    (HERE/"isogeny_module_controls.json").write_text(json.dumps(out,indent=2)+"\n")
    print(json.dumps({key:out[key] for key in ["claim","status","discriminant","norm_x",
        "rank_six_polynomial_matrix_degrees","uniform_matrix_generator_degree_bound"]},indent=2))


if __name__=="__main__":
    main()
