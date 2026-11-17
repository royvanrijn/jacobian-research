"""Exact controls for the written pole-order <=3 exclusion.

The all-family proof is POLE_REDUCTION.md. This checks algebra and exact
linear certificates; it neither formally certifies the proof nor supplies
a compatible completion or a finite map.
"""
from pathlib import Path
from math import comb
from functools import reduce
from time import monotonic
import hashlib,json,platform
import sympy as S

HERE=Path(__file__).resolve().parent
CHECKS=[]

def check(name,condition):
    if not condition:
        raise AssertionError(name)
    CHECKS.append(name)

def zero(name,expr):
    check(name,S.cancel(expr)==0)

def wedge(F,G,variables):
    u,v=variables
    return S.cancel(S.diff(F,u)*S.diff(G,v)-S.diff(F,v)*S.diff(G,u))

def independent_matrix(record):
    """Differentiate rational quotients via polynomial pairs, then binomial shift."""
    x,y=S.symbols("x y")
    m,n=record["pole_orders"]
    B=S.Poly(x-y*y,x,y,domain=S.QQ)
    Bx,By=B.diff(x),B.diff(y)
    am=[S.Poly(S.sympify(a),x,y,domain=S.QQ) for a in record["basis_m"]]
    cn=[S.Poly(S.sympify(c),x,y,domain=S.QQ) for c in record["basis_n"]]
    columns=[]
    for i,j in record["columns"]:
        A,C=am[i],cn[j]
        Ux,Uy=B*A.diff(x)-m*A*Bx,B*A.diff(y)-m*A*By
        Vx,Vy=B*C.diff(x)-n*C*Bx,B*C.diff(y)-n*C*By
        W=(Ux*Vy-Uy*Vx).exquo(B)
        shifted={}
        for (ix,iy),value in W.terms():
            for iz in range(min(ix,m+n-3)+1):
                key=(iz,2*(ix-iz)+iy)
                shifted[key]=shifted.get(key,0)+value*comb(ix,iz)
        columns.append({k:v for k,v in shifted.items() if v})
    effective_rows=sorted({power for col in columns for power in col})
    stored_effective_rows=sorted({tuple(record["row_monomials"][i])
                                  for i,j,value in record["sparse_matrix"] if S.Rational(value)})
    check(f"independent nonzero rows {m},{n}",effective_rows==stored_effective_rows)
    # The producer retains a harmless zero row from Poly(0).terms(); verify padding too.
    rows=list(map(tuple,record["row_monomials"]))
    M=S.MutableSparseMatrix(len(rows),len(columns),{
        (i,j):columns[j][power] for j in range(len(columns))
        for i,power in enumerate(rows) if power in columns[j]})
    stored=S.MutableSparseMatrix(len(rows),len(columns),{
        (i,j):S.Rational(c) for i,j,c in record["sparse_matrix"]})
    check(f"all independently reconstructed columns {m},{n}",M==stored)
    count=record["kernel_dimension"]
    if "kernel_basis" in record:
        K=S.Matrix([[S.Rational(v) for v in row] for row in record["kernel_basis"]])
    else:
        K=S.MutableSparseMatrix(count,len(columns),{
            (i,j):S.Rational(v) for i,row in enumerate(record["sparse_kernel_basis"])
            for j,v in row})
    check(f"exact nullspace {m},{n}",M*K.T==S.zeros(M.rows,K.rows))
    # Free-column entries give an explicit diagonal minor for kernel row rank.
    # This avoids a redundant dense elimination on hundreds of sparse rows.
    pivots=record["pivot_columns"]
    free=[j for j in range(M.cols) if j not in pivots]
    free_positions={j:i for i,j in enumerate(free)}
    check(f"kernel free-column diagonal is nonzero {m},{n}",
          len(free)==count and all(K[i,j]!=0 for i,j in enumerate(free)))
    check(f"kernel free-column off-diagonal vanishes {m},{n}",
          all(i==free_positions[j] for (i,j),value in K.todok().items()
              if j in free_positions and value))
    check(f"pivot-column lower bound {m},{n}",
          M[:,pivots].rank()==record["constraint_rank"])
    check(f"complete rank-nullity {m},{n}",
          record["constraint_rank"]+count==M.cols)
    return dict(pole_orders=[m,n],columns=M.cols,
                rank=record["constraint_rank"],kernel=count)

def main():
    start=monotonic()
    z,w=S.symbols("z w")
    p,q=S.Function("p")(w),S.Function("q")(w)
    m,k=S.symbols("m k",integer=True,positive=True)
    zero("general divisible-order derivative identity",
         k*m*S.diff(p,w)*q-m*p*S.diff(q,w)
         +m*p**(k+1)*S.diff(q/p**k,w))
    X,Y,c=S.symbols("X Y c")
    zero("general target shear has Jacobian one",wedge(X,Y-c*X**k,(X,Y))-1)

    # Enumerate all terminal pole pairs at height <=4, including strictly lowering
    # divisible steps. Exhaustion here is finite; the general reduction is written.
    terminal={}
    def reduce_orders(a,b):
        a,b=sorted((a,b))
        if (a,b)==(1,1):
            return {(1,1)}
        if b%a:
            return {(a,b)}
        return set().union(*(reduce_orders(a,j) for j in range(1,b)))
    for a in range(1,5):
        for b in range(a,5):
            terminal[f"{a},{b}"]=sorted(reduce_orders(a,b))
    check("terminal pairs at pole bound two",
          set().union(*(set(map(tuple,terminal[f"{a},{b}"]))
                        for a in range(1,3) for b in range(a,3)))=={(1,1)})
    check("terminal pairs at pole bound three",
          set().union(*(set(map(tuple,terminal[f"{a},{b}"]))
                        for a in range(1,4) for b in range(a,4)))=={(1,1),(2,3)})
    check("next nondivisible pair at height four",
          reduce_orders(3,4)=={(3,4)})

    T,A,C,D=S.symbols("T A C D")
    P=T*T+A
    Q=T**3+S.Rational(3,2)*T*A+C
    for (u,v),expected in [
        ((T,A),-S.Rational(3,2)*A),((T,C),2*T),((A,C),1)]:
        zero(f"global normal-form wedge coefficient {u},{v}",
             wedge(P,Q,(u,v))-expected)
    zero("normalization residual has no pole-five contribution from the core",
         wedge(P,T**3+S.Rational(3,2)*T*A,(T,A))+S.Rational(3,2)*A)

    tau,aa,cc=[S.Function(name)(w) for name in ("tau","a","c")]
    leading=(-S.Rational(3,2)*aa**2*S.diff(tau,w)
             +S.Rational(3,2)*tau*aa*S.diff(aa,w)
             +2*tau*S.diff(tau,w)*cc-2*tau**2*S.diff(cc,w))
    zero("conic relation is the pole-four differential equation",
         leading+2*tau**3*S.diff(cc/tau-S.Rational(3,8)*(aa/tau)**2,w))

    ss,tt,uu,vv,kap=S.symbols("s t u v kappa")
    root_controls=[]
    cubic_families=[
        ("three_simple",ss*tt*(ss-tt),uu*ss*tt*(ss-tt),
         (S.Rational(3,8)*uu**2+kap)*ss*tt*(ss-tt),tt),
        ("double_and_simple",ss**2*tt,ss*tt*(uu*ss+vv*tt),
         S.Rational(3,8)*tt*(uu*ss+vv*tt)**2+kap*ss**2*tt,tt),
        ("triple",ss**3,ss**2*(uu*ss+vv*tt),
         S.Rational(3,8)*ss*(uu*ss+vv*tt)**2+kap*ss**3,ss)]
    for name,ta,al,cl,factor in cubic_families:
        zero("cubic conic relation "+name,8*ta*cl-3*al**2-8*kap*ta**2)
        check("common homogeneous factor "+name,
              all(S.rem(S.Poly(f,ss,tt),S.Poly(factor,ss,tt)).is_zero
                  for f in (ta,al,cl)))
        root_controls.append(dict(partition=name,common_factor=str(factor)))
    even=(ss**2,ss*tt,S.Rational(3,8)*tt**2)
    zero("even-degree conic control",8*even[0]*even[2]-3*even[1]**2)
    check("even-degree control has no common zero",
          reduce(S.gcd,[S.Poly(f,ss,tt) for f in even]).total_degree()==0)

    # Generic local linear terms confirm the maximal-base-point order exactly.
    coeff=S.symbols("a0:6")
    t0=coeff[0]*z+coeff[1]*w
    a0=coeff[2]*z+coeff[3]*w
    c0=coeff[4]*z+coeff[5]*w
    sections=[z*t0*t0+z*z*a0,
              t0**3+S.Rational(3,2)*z*t0*a0+z*z*c0,z**3]
    orders=[min(sum(power) for power,value in S.Poly(f,z,w).terms()) for f in sections]
    check("proper base multiplicity three from common symbol zero",min(orders)==3)
    v=S.symbols("v")
    for name,numerator,expected in [
        ("outside_E",S.Integer(1),-2),
        ("smooth_E_contact_four",w-z**4,-1)]:
        pulled=S.cancel((numerator/z**3).subs(w,z*v)*z)
        num,den=S.fraction(pulled)
        order=min(power[0] for power,value in S.Poly(num,z).terms())-S.degree(den,z)
        check("forbidden proper-exceptional differential pole "+name,order==expected)

    # A nonzero pole-three example where the zero divisor has multiplicity two.
    # It explains why smoothness of E is essential; it is not a counterexample.
    t0=w*(z+w*w)/z
    a0=w*w/z
    c0=S.Rational(3,8)*w/z
    p0=t0*t0+a0
    q0=t0**3+S.Rational(3,2)*t0*a0+c0
    volume=S.cancel(z**3*wedge(p0,q0,(z,w)))
    zero("degenerate-volume control numerator",
         volume-S.Rational(3,8)*w*w*(2*w*w+2*z+1))
    check("degenerate control has double zero, not smooth ramification",
          min(sum(power) for power,value in S.Poly(volume,z,w).terms())==2)

    # Exact multiple-fibre control in the two-chart affine-line torsor.
    x,y,u,v=S.symbols("x y u v")
    source_a=x*x+y
    source_d=x+source_a**2
    zero("source coordinate Jacobian",wedge(source_a,source_d,(x,y))+1)
    zero("source coordinate inverse x",source_d-source_a**2-x)
    zero("source coordinate inverse y",source_a-(source_d-source_a**2)**2-y)
    change={x:1/u,y:u**3*v-u*u}
    for expr,expected in [
        (y,u**3*v-u*u),(x*y,u*u*v-u),(x*x*y,u*v-1),(x**3*y+x,v)]:
        zero("global torsor function "+str(expr),expr.subs(change)-expected)
    a_unit=S.cancel(u**2*source_a.subs(change))
    d_unit=S.cancel(u**4*source_d.subs(change))
    zero("first multiple-fibre unit",a_unit-(1-u**4+u**5*v))
    zero("second multiple-fibre unit",d_unit-((1-u**4+u**5*v)**2+u**3))
    check("infinity fibre multiplicities are exactly two and four",
          a_unit.subs(u,0)==1 and d_unit.subs(u,0)==1)

    reports=[]
    for filename in ("pole_layers.json","pole_frontier.json"):
        payload=json.loads((HERE/filename).read_text())
        if filename=="pole_layers.json":
            check("initial matrix producer hash",
                  payload["checker_sha256"]==hashlib.sha256(
                      (HERE/"compute_pole_layers.py").read_bytes()).hexdigest())
            records=payload["layers"]
        else:
            for path,digest in payload["source_sha256"].items():
                check("frontier producer hash "+path,
                      hashlib.sha256((HERE/path).read_bytes()).hexdigest()==digest)
            records=[payload["layer"]]
        reports.extend(independent_matrix(record) for record in records)
    check("expected exact constraint ranks",
          [r["rank"] for r in reports]==[19,33,70])

    manifest=json.loads((HERE/"sources/wright-connection-manifest.json").read_text())
    for source in manifest["sources"]:
        raw=(HERE/"sources"/source["file"]).read_bytes()
        check("literature-lead source hash "+source["file"],
              hashlib.sha256(raw).hexdigest()==source["sha256"] and len(raw)==source["bytes"])
    check("literature lead is not adopted as a theorem",
          manifest["adopted_as_theorem"] is False)

    report=dict(as_of="2026-11-14",outcome="PASS",check_count=len(CHECKS),checks=CHECKS,
                root_partition_controls=root_controls,terminal_pole_controls=terminal,
                independently_verified_linear_layers=reports,
                written_result="Compatible quartic finite maps have no target-coordinate pole presentation of height <=3",
                subsequently_closed_pair=[3,4], subsequent_full_exclusion="PF2RULED4",
                multiple_fibre_control=[2,4],
                assurance=dict(all_family_geometry="written proof, not formally certified",
                               finite_map_constructed=False,actual_ruled_compatibility_proved=False,
                               higher_poles_excluded=False,independent_full_reproof=False,
                               formal_verification=False,external_review=False),
                source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest()
                               for p in [Path(__file__),HERE/"POLE_REDUCTION.md",
                                         HERE/"WRIGHT_CONNECTION.md"]},
                python=platform.python_version(),sympy=S.__version__,
                elapsed_seconds=round(monotonic()-start,6))
    (HERE/"pole_reduction_receipt.json").write_text(json.dumps(report,indent=2)+"\n")
    print(f"PASS {len(CHECKS)} exact controls; m<=3 reduction verified; full pole-four exclusion is PF2RULED4")

if __name__=="__main__":
    main()

