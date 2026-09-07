#!/usr/bin/env python3
"""Explicit soluble blocks in the fixed-cubic pencil; bounded retrospective replay."""
import argparse
from pathlib import Path
import retrospective as r
from cubic_bridge import Cubic
from cover_experiment import mul, sub, evaluate, sqrtq

HERE = Path(__file__).resolve().parent
PROTOCOL = HERE/"LINEAR_TWIST_BLOCK_PROTOCOL.json"
INPUT = r.OUT/"rank_jump_local_collision_inputs_v1.json"
OUTPUT = r.OUT/"rank_jump_linear_twist_blocks_v1.json"


def add(a,b):
    return sub(a,[-x for x in b])


def scale(a,k):
    return [k*x for x in a]


def bindings():
    paths = [PROTOCOL, INPUT, Path(__file__), HERE/"cubic_bridge.py",
             HERE/"cover_experiment.py", HERE/"retrospective.py"]
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def twist_identity(A,B,a):
    g=[r.F(1),-a];D=[r.F(1),r.F(0),A,B]
    X=[a,r.F(0),A*a+B]
    u=[r.F(0),r.F(1)]
    a2=[r.F(0),2*A];a4=[A,3*B,A*A];a6=[B,r.F(0),A*B,-B*B]
    rhs=add(add(mul(mul(X,X),X),mul(mul(a2,g),mul(X,X))),
            add(mul(mul(a4,mul(g,g)),X),mul(a6,mul(mul(g,g),g))))
    assert rhs == scale(mul(D,D),a**3+A*a+B)
    return {"a":str(a),"g_coefficients":list(map(str,g)),
            "twist_x_coefficients":list(map(str,X)),
            "twist_y_squared_coefficients":list(map(str,rhs))}


def build(check=False):
    inp=r.read(INPUT)
    A,B=map(r.F,inp["anchor"]["short_model_ainvariants"][3:])
    points=[tuple(map(r.F,P)) for P in inp["anchor"]["known_points_on_short_model"]]
    assert len(points)==20 and len({p for p,q in points})==20
    assert inp["anchor"]["known_kummer_dimension"]==20
    K=Cubic(A,B);th=K.theta;th2=K.square(th)
    identities=[twist_identity(A,B,a) for a,b in points]
    rows=[]
    for u in range(-3,4):
        D=1+A*u*u+B*u**3
        assert D
        alpha=K.add(th,K.scale(th2,u))
        gamma=K.sub(K.one,K.scale(th,u))
        eta=K.scale(gamma,D)
        kappa=K.add(K.add(K.one,K.scale(th,u)),K.scale(K.add(K.scalar(A),th2),u*u))
        assert K.mul(kappa,gamma)==K.scalar(D)
        assert K.norm(eta)==D**4
        tests=[]
        for i,(a,b) in enumerate(points):
            assert b*b==a**3+A*a+B
            g=1-a*u
            root=sqrtq(g)
            beta=K.sub(K.scalar(a),th)
            # Universal equality before the square test; it tracks the class shift.
            assert g
            X=(a+(A*a+B)*u*u)/g
            lhs=K.sub(K.scalar(X),alpha)
            rhs=K.scale(K.mul(beta,kappa),1/g)
            assert lhs==rhs
            item={"basis_index":i,"one_minus_a_u":str(g),"rational_square":root is not None}
            if root is not None:
                assert root
                Y=b*D/(root**3)
                assert Y*Y==X**3+2*A*u*X*X+(A+3*B*u+A*A*u*u)*X+B+A*B*u*u-B*B*u**3
                # Exact square correction: x-alpha = beta*eta/(z*gamma)^2.
                correction=K.scale(gamma,root)
                assert K.mul(lhs,K.square(correction))==K.mul(beta,eta)
                item.update({"square_root":str(root),"transported_point":[str(X),str(Y)],
                             "kummer_square_correction":list(map(str,correction))})
            tests.append(item)
        k=sum(x["rational_square"] for x in tests)
        local=next(x for x in inp["rows"] if int(x["parameter_u"])==u)
        ct=next(x["matrix"] for x in inp["ct"] if x["u"]==u)
        assert local["all_local_kummer_images_complete"]
        n=local["W_u_dimension"]
        assert n==len(ct)
        ct_rank=r.rank([r.pack(row) for row in ct])
        radical=n-ct_rank
        transport_cap=min(20,radical+1)
        assert k<=transport_cap
        rows.append({"u":u,"D":str(D),"common_kummer_shift":list(map(str,eta)),
                     "tests":tests,"soluble_transports":k,
                     "inherited_local_dimension":n,"inherited_CT_rank":ct_rank,
                     "inherited_CT_radical_dimension":radical,
                     "independent_anchor_class_transport_cap":transport_cap,
                     "cap_scope":"Necessary bound on the number of independent inherited anchor Kummer classes simultaneously transported, valid for every choice of representatives. Not an upper rank of E_u or of an unsaturated point subgroup.",
                     "rank_lower_bound_from_transport":k if u==0 else max(0,k-1),
                     "bound_scope":"Only the displayed transporter; a zero lower bound is not a rank-zero assertion."})
    out={"schema":"rank-jump.linear-twist-blocks.v1","bindings":bindings(),
         "anchor_A":str(A),"anchor_B":str(B),"twist_identities":identities,
         "retrospective_specializations":rows,
         "summary":{"rational_square_tests":140,"universal_twist_identities":20,
                    "transport_counts":[x["soluble_transports"] for x in rows],
                    "nonzero_u_transports":sum(x["soluble_transports"] for x in rows if x["u"]),
                    "independence_theorem":"k independent anchor Kummer classes with nonzero square 1-a_i*u imply rank E_u(Q) >= k-1, provided the cubic is irreducible and D(u) != 0.",
                    "splitting_event":"At u=0 the shared correspondence b^2=f(a), z^2=1-au splits into two copies of E_0."},
         "boundary":"The retained anchor points are oracle inputs; this is not an Agent 1 selector, a new curve search, or a whole-curve upper rank."}
    if check:
        assert r.read(OUTPUT)==out
        print("PASS exact linear-twist block replay")
    else:
        r.write_new(OUTPUT,out)
    print(out["summary"])


def verify():
    from sage.all import QQ, PolynomialRing
    ring=PolynomialRing(QQ,names=("A","B","a","u","x"))
    A,B,a,u,x=ring.gens()
    g=1-a*u;D=1+A*u*u+B*u**3;f=a**3+A*a+B
    X=a+(A*a+B)*u*u
    rhs=X**3+2*A*u*g*X*X+(A+3*B*u+A*A*u*u)*g*g*X+(B+A*B*u*u-B*B*u**3)*g**3
    assert rhs==f*D*D
    # Quartic r^2=(1-au)f(a) <-> E_u, via a=(x-Bu^2)/(1+Au^2+xu).
    F=ring.fraction_field()
    den=1+A*u*u+x*u
    inva=F((x-B*u*u)/den)
    Fu=x**3+2*A*u*x*x+(A+3*B*u+A*A*u*u)*x+B+A*B*u*u-B*B*u**3
    assert (1-inva*u)*(inva**3+A*inva+B)==Fu*D*D/(den**4)
    data=r.read(OUTPUT)
    assert data["bindings"]==bindings()
    print("PASS Sage universal twist identity and inverse quartic correspondence")


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("mode",choices=("build","check","verify"))
    mode=parser.parse_args().mode
    verify() if mode=="verify" else build(mode=="check")
