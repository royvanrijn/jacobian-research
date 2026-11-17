#!/usr/bin/env python3
"""Retain a nonzero-target collision family and its affine endpoint obstructions."""
from pathlib import Path
import hashlib,json
import sympy as sp
from itertools import combinations
HERE=Path(__file__).resolve().parent
def main():
    x,y,z,k=sp.symbols("x y z k")
    U=1+x*y
    F=sp.Matrix([U**3*z+y**2*U*(4+3*x*y),
                 y+3*x*U**2*z+3*x*y**2*(4+3*x*y),
                 2*x-3*x*x*y-x**3*z])
    a=sp.Matrix([(1-k)/2,-2,-4*(1+3*k)/k**2])
    b=sp.Matrix([(1+k)/2,-2,-4*(1-3*k)/k**2])
    fa=F.subs(dict(zip((x,y,z),a))).applyfunc(sp.cancel)
    fb=F.subs(dict(zip((x,y,z),b))).applyfunc(sp.cancel)
    assert fa==fb==sp.Matrix([0,-2,(1-k**2)/(2*k**2)]).applyfunc(sp.cancel)
    A=F.jacobian((x,y,z)).subs(dict(zip((x,y,z),a)))
    L=sp.Matrix.hstack(b-a,sp.Matrix([0,1,0]),sp.Matrix([0,0,1]))
    assert sp.cancel(A.det())==-2 and L.det()==k
    sub=dict(zip((x,y,z),a+L*sp.Matrix([x,y,z])))
    H=((A*L).inv()*(F.subs(sub,simultaneous=True)-fa)).applyfunc(sp.expand)
    assert H.subs({x:0,y:0,z:0})==sp.zeros(3,1)
    assert H.subs({x:1,y:0,z:0}).applyfunc(sp.cancel)==sp.zeros(3,1)
    assert H.jacobian((x,y,z)).subs({x:0,y:0,z:0}).applyfunc(sp.cancel)==sp.eye(3)
    axis=sp.Poly(H[0].subs({y:0,z:0}),x)
    c2=sp.cancel(axis.nth(2));c4=sp.cancel(axis.nth(4))
    assert sp.limit(k*c2,k,0)==12 and sp.limit(c4/k**2,k,sp.oo)==36
    coefficients=[]
    for i,h in enumerate(H):
        for exponent,c in sp.Poly(h,x,y,z).terms():
            if sum(exponent)>=2:
                coefficients.append(dict(row=i,exponent=list(exponent),coefficient=str(sp.cancel(c))))
    receipt=dict(status="verified-exploratory-identity",
        source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        parameter_domain="k nonzero",source_a=list(map(str,a)),source_b=list(map(str,b)),
        common_target=list(map(str,fa)),generic_support=len(coefficients),
        coefficient_rows=coefficients,first_axis=str(axis.as_expr()),
        x2_coefficient=str(c2),x4_coefficient=str(c4),
        k0_obstruction="x^2 coefficient in H1 has the invariant pole 12/k",
        kinfinity_obstruction="x^4 coefficient in H1 has the invariant leading term 36*k^2",
        scope="No finite coefficientwise limit at k=0 or infinity after any meromorphic frame fixing e1 for this specified collision pair family. Finite-parameter sparsity and other collision pair families remain open.")
    # Canonicalize the full quartic collision axis by its coefficient matrix.
    C=sp.Matrix([[sp.Poly(entry,x).nth(j) for j in (2,3,4)]
                 for entry in H.subs({y:0,z:0})])
    D=sp.Matrix([[-1,0,0],[1,-1,0],[0,1,-1]])
    assert sp.factor(C.det())==248832*k**2
    M=(C*D.inv()).applyfunc(sp.cancel)
    assert M[:,0]==sp.Matrix([1,0,0])
    G=(M.inv()*H.subs(dict(zip((x,y,z),M*sp.Matrix([x,y,z]))),
                        simultaneous=True)).applyfunc(sp.expand)
    assert G.subs({y:0,z:0}).applyfunc(sp.cancel)==sp.Matrix([x-x*x,x*x-x**3,x**3-x**4])
    valuation_rows=[]
    for i,entry in enumerate(G):
        for e,c in sp.Poly(entry,x,y,z).terms():
            if sum(e)<2:continue
            c=sp.cancel(c);num,den=c.as_numer_denom()
            P=sp.Poly(num,k);Q=sp.Poly(den,k)
            v0=min(d[0] for d,c in P.terms())-min(d[0] for d,c in Q.terms())
            vi=Q.degree()-P.degree()
            valuation_rows.append(dict(row=i,exponent=list(e),coefficient=str(c),
                charges=[e[1]-int(i==1),e[2]-int(i==2)],zero_valuation=int(v0),infinity_valuation=int(vi)))
    def farkas(endpoint):
        strongest={}
        for row in valuation_rows:
            charge=tuple(row["charges"])
            if charge not in strongest or row[endpoint]<strongest[charge][endpoint]:
                strongest[charge]=row
        rows=list(strongest.values())
        def verify(selected,multipliers):
            assert all(m>=0 for m in multipliers) and any(multipliers)
            assert all(sum(m*r["charges"][j] for r,m in zip(selected,multipliers))==0 for j in range(2))
            value=sum(m*r[endpoint] for r,m in zip(selected,multipliers))
            assert value<0
            return dict(rows=selected,multipliers=list(map(int,multipliers)),negative_constant=int(value))
        for r in rows:
            if r["charges"]==[0,0] and r[endpoint]<0:return verify([r],[1])
        for r,s in combinations(rows,2):
            a,b=r["charges"];c,d=s["charges"]
            if a*d-b*c==0 and a*c+b*d<0:
                weights=[abs(c),abs(a)] if a and c else [abs(d),abs(b)]
                if sum(m*q[endpoint] for q,m in zip((r,s),weights))<0:
                    return verify([r,s],weights)
        for r,s,trow in combinations(rows,3):
            a,b=r["charges"];c,d=s["charges"];e,f=trow["charges"]
            weights=[c*f-d*e,e*b-f*a,a*d-b*c]
            if all(q<=0 for q in weights):weights=[-q for q in weights]
            if all(q>=0 for q in weights) and any(weights):
                if sum(m*q[endpoint] for q,m in zip((r,s,trow),weights))<0:
                    return verify([r,s,trow],weights)
        raise AssertionError("No exact infeasibility certificate found")
    pole_order=min(row["zero_valuation"] for row in valuation_rows)
    principal={}
    for row in valuation_rows:
        if row["zero_valuation"]!=pole_order:continue
        num,den=sp.cancel(sp.sympify(row["coefficient"])).as_numer_denom()
        P=sp.Poly(num,k);Q=sp.Poly(den,k)
        ne=min(e[0] for e,c in P.terms());de=min(e[0] for e,c in Q.terms())
        coefficient=P.nth(ne)/Q.nth(de)
        principal.setdefault(tuple(row["exponent"]),[0,0,0])[row["row"]]=coefficient
    principal_matrix=sp.Matrix.hstack(*[sp.Matrix(principal[e]) for e in sorted(principal)])
    trace=sp.cancel(sp.trace(H.jacobian((x,y,z)).subs({y:0,z:0})))
    axis2=H.subs({x:2,y:0,z:0}).applyfunc(sp.cancel)
    image_trace=sp.cancel(sp.trace(H.jacobian((x,y,z))).subs(dict(zip((x,y,z),axis2)),simultaneous=True))
    axis_trace=sp.cancel(trace.subs(x,2))
    assert axis_trace==864*k**2-2736*k+723
    assert sp.limit(k**3*image_trace,k,0)==-49460507836416
    receipt["frame_invariants"]=dict(
        axis_jacobian_trace=str(trace),
        trace_at_axis2=str(axis_trace),trace_at_image_axis2=str(image_trace),
        zero_image_trace_valuation=-3,zero_image_trace_leading_coefficient="-49460507836416",
        scope="Trace JH(2e1) and trace JH(H(2e1)) are invariant under every conjugation fixing e1. Their endpoint poles prohibit every finite meromorphic affine-frame limit for this pair family.",
        trace_x4_coefficient=str(sp.factor(sp.Poly(trace,x).nth(4))),
        zero_principal_valuation=int(pole_order),
        zero_principal_output_rank=int(principal_matrix.rank()),
        zero_principal_columns=[dict(exponent=list(e),vector=list(map(str,principal[e]))) for e in sorted(principal)])
    receipt["axis_normalization"]=dict(
        matrix=[[str(v) for v in row] for row in M.tolist()],
        generic_support=len(valuation_rows),axis=["x-x^2","x^2-x^3","x^3-x^4"],
        zero_farkas=farkas("zero_valuation"),infinity_farkas=farkas("infinity_valuation"),
        boundary="The valuation certificates exclude diagonal limits in this axis frame. The separate invariant trace-at-image pole now excludes all meromorphic frame limits at zero, and trace-at-axis excludes infinity.")
    (HERE/"nonzero-target-probe.json").write_text(json.dumps(receipt,indent=2,sort_keys=True)+chr(10))
    print(json.dumps({k:v for k,v in receipt.items() if k not in ("coefficient_rows","axis_normalization")}))
if __name__=="__main__":main()

