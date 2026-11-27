#!/usr/bin/env python3
"""Exact rational Hermite charts and critical-value norm certificates."""
from fractions import Fraction as F
from math import factorial
from pathlib import Path
import hashlib,json,platform
import verify_polynomial_involution_components as b

if not __debug__:
    raise RuntimeError("Hermite phase checker requires active assertions")
HERE=Path(__file__).resolve().parent
c,m,l=b.c,b.m,b.l


def trim(p):
    p=list(map(F,p))
    while len(p)>1 and not p[-1]:p.pop()
    return p or [F(0)]
def add(p,q):
    out=[F(0)]*max(len(p),len(q))
    for i,x in enumerate(p):out[i]+=x
    for i,x in enumerate(q):out[i]+=x
    return trim(out)
def scale(p,k):return trim([x*k for x in p])
def sub(p,q):return add(p,scale(q,-1))
def mul(p,q):
    out=[F(0)]*(len(p)+len(q)-1)
    for i,x in enumerate(p):
        for j,y in enumerate(q):out[i+j]+=x*y
    return trim(out)
def power(p,k):
    out=[F(1)]
    for _ in range(k):out=mul(out,p)
    return out
def divmod_poly(p,q):
    p,q=trim(p),trim(q);assert q!=[0]
    quotient=[F(0)]*max(1,len(p)-len(q)+1)
    while p!=[0] and len(p)>=len(q):
        j=len(p)-len(q);value=p[-1]/q[-1];quotient[j]+=value
        p=sub(p,[0]*j+scale(q,value))
    return trim(quotient),p
def mod(p,q):return divmod_poly(p,q)[1]
def gcd_poly(p,q):
    while q!=[0]:p,q=q,mod(p,q)
    return scale(p,1/p[-1])
def inverse_mod(p,q):
    r0,r1=p,q;s0,s1=[F(1)],[F(0)]
    while r1!=[0]:
        quotient,remainder=divmod_poly(r0,r1)
        r0,r1=r1,remainder;s0,s1=s1,sub(s0,mul(quotient,s1))
    assert len(r0)==1
    answer=mod(scale(s0,1/r0[0]),q)
    assert mod(mul(p,answer),q)==[1]
    return answer
def derivative(p):return trim([i*p[i] for i in range(1,len(p))])
def evaluate(p,x):
    out=F(0)
    for y in p[::-1]:out=out*x+y
    return out
def compose(p,q):
    out=[F(0)]
    for y in p[::-1]:out=add(mul(out,q),[y])
    return out
def from_roots(roots):
    out=[F(1)]
    for x in roots:out=mul(out,[-x,1])
    return out


def hermite_sign(plus,minus):
    U,V=power(plus,2),power(minus,2)
    assert gcd_poly(U,V)==[1]
    answer=mod(sub([1],scale(mul(U,inverse_mod(U,V)),2)),mul(U,V))
    assert mod(sub(answer,[1]),U)==[0]
    assert mod(add(answer,[1]),V)==[0]
    assert mod(derivative(answer),mul(plus,minus))==[0]
    return answer


def chart(n,a,j,plus,minus,v,residual):
    assert a>=2 and n>2*a and 0<=j<=a//2 and v
    plus,minus,residual=map(trim,(plus,minus,residual))
    assert len(plus)==j+1 and len(minus)==a-j+1
    assert len(residual)==n-2*a+1
    assert plus[-1]==minus[-1]==residual[-1]==1
    P=mul(plus,minus);Q=hermite_sign(plus,minus)
    H=add(scale(Q,v),mul(power(P,2),residual))
    assert len(H)==n+1 and H[-1]==1
    assert mod(sub(H,[v]),power(plus,2))==[0]
    assert mod(add(H,[v]),power(minus,2))==[0]
    A,rest=divmod_poly(derivative(H),scale(P,n))
    assert rest==[0] and len(A)==n-a and A[-1]==1
    return P,Q,H,A


def matmul(A,B):
    return [[sum((A[i][k]*B[k][j] for k in range(len(A))),F(0))
             for j in range(len(B[0]))] for i in range(len(A))]
def determinant(A):
    A=[row[:] for row in A];value=F(1);n=len(A)
    for j in range(n):
        pivot=next((i for i in range(j,n) if A[i][j]),None)
        if pivot is None:return F(0)
        if pivot!=j:A[j],A[pivot]=A[pivot],A[j];value=-value
        x=A[j][j];value*=x
        for i in range(j+1,n):
            ratio=A[i][j]/x
            for k in range(j+1,n):A[i][k]-=ratio*A[j][k]
    return value


def critical_value_norm(A,H):
    d=len(A)-1
    columns=[mod([0]*k+power(H,2),A) for k in range(d)]
    M=[[columns[j][i] if i<len(columns[j]) else F(0) for j in range(d)]
       for i in range(d)]
    I=[[F(i==j) for j in range(d)] for i in range(d)]
    powers=I;traces=[F(0)]
    for k in range(1,d+1):
        powers=matmul(powers,M);traces.append(sum(powers[i][i] for i in range(d)))
    descending=[F(1)]
    for k in range(1,d+1):
        descending.append(-sum(descending[k-i]*traces[i] for i in range(1,k+1))/k)
    polynomial=descending[::-1]
    # Different exact method, at d+1 points: both sides have degree d.
    for y in range(d+1):
        matrix=[[F(y)*(i==j)-M[i][j] for j in range(d)] for i in range(d)]
        assert determinant(matrix)==evaluate(polynomial,y)
    return polynomial


def admissible(P,H,A,v,C=None):
    if not v:return False
    if gcd_poly(P,derivative(P))!=[1]:return False
    if gcd_poly(A,derivative(A))!=[1] or gcd_poly(P,A)!=[1]:return False
    C=critical_value_norm(A,H) if C is None else C
    return bool(evaluate(C,0) and evaluate(C,v*v) and gcd_poly(C,derivative(C))==[1])


def chart_controls():
    rows=[]
    for n,a in ((5,2),(6,2),(7,3),(8,3),(9,4),(10,4),(11,5),(12,5)):
        for j in range(a//2+1):
            plus=from_roots(range(j));minus=from_roots(range(j,a))
            residual=from_roots(a+2+3*k for k in range(n-2*a))
            for denominator in (101,103,107):
                v=F(1,denominator)
                P,Q,H,A=chart(n,a,j,plus,minus,v,residual)
                C=critical_value_norm(A,H)
                if admissible(P,H,A,v,C):break
            else:raise AssertionError(("no valid finite witness",n,a,j))
            assert gcd_poly(H,derivative(H))==[1]
            rows.append({"n":n,"a":a,"phase_smaller_bin":j,"dimension":n-a-1,
                         "plus":list(map(str,plus)),"minus":list(map(str,minus)),
                         "residual":list(map(str,residual)),"v":str(v),
                         "H":list(map(str,H)),"critical_quotient":list(map(str,A)),
                         "other_squared_critical_value_polynomial":list(map(str,C)),
                         "all_exact_open_conditions":True})
    return rows


def degree10_surfaces():
    out=[]
    z=[0,1];P=mul(z,[-1,1])
    for j in (0,1):
        plus=[1] if j==0 else z
        minus=P if j==0 else [-1,1]
        for t,v in ((F(1),F(1,101)),(F(-2),F(2,103)),(F(3,2),F(1,107))):
            residual=[t,1]
            _,Q,H,A=chart(5,2,j,plus,minus,v,residual)
            assert Q==([-1] if j==0 else [1,0,-6,4])
            expected=[-2*t+(12*v if j else 0),4*t-3,5]
            assert scale(A,5)==expected
            assert admissible(P,H,A,v)
            new_t=-1-t;new_v=-v if j==0 else v
            _,_,new_H,_=chart(5,2,j,plus,minus,new_v,[new_t,1])
            assert new_H==scale(compose(H,[1,-1]),-1)
            coordinates=((2*t+1)/v,v*v) if j==0 else ((2*t+1)**2,v)
            new_coordinates=((2*new_t+1)/new_v,new_v*new_v) if j==0 else ((2*new_t+1)**2,new_v)
            assert coordinates==new_coordinates
            out.append({"phase":j,"t":str(t),"v":str(v),
                        "quotient_coordinates":list(map(str,coordinates))})
    return out


def size_controls():
    checked=0;hook=0
    for n in range(5,25):
        for a in range(2,(n-1)//2+1):
            s=n-a-1;L=(s+2)*(s+1)
            sizes=[]
            descriptors=b.component_descriptors(2*n,n,a,s)
            assert len(descriptors)==a//2+1
            for j in range(a//2+1):
                C=b.ordered_count(n,j,a-j,s)
                size=L*2**s*(2-(2*j==a))*C;sizes.append(size)
                row=next(row for row in descriptors if row["primitive_degree_partition"]==
                         tuple(sorted((n-2*j,n-2*(a-j)))))
                assert row["inner_degree"]==n and row["outer_degree"]==2
                assert row["literal_component_size"]==size
                assert row["group_order"]==2*factorial(n)**2
            total=L*b.ordered_count(2*n,n,a,s)
            assert sum(sizes)==total
            if n<=10:
                passport=[b.involution_type(2*n,n),b.involution_type(2*n,a)]+[
                    b.involution_type(2*n,1)]*s
                assert c.character_total(passport)*c.type_placements(passport)==total
                hook+=1
            if a==2:
                assert sizes==[(n-1)*(n-2)*2**(n-3)*(n-3)*n**(n-3),
                               (n-1)*(n-2)*2**(n-3)*n**(n-2)]
            checked+=1
    return {"all_parameter_count_domains":checked,"max_n":24,
            "independent_hook_domains":hook,
            "scope":"Arithmetic controls only; all-parameter sizes and connectedness use J1-J3."}


def degree10_orbits():
    s=2;rows=b.component_descriptors(10,5,2,s);union=set();out=[]
    for row in rows:
        seed=b.component_seed(row,s)
        orbit,_=c.closure(seed)
        assert len(orbit)==row["literal_component_size"] and not union&orbit
        union.update(orbit)
        normal=b.l.group_from_transposition(seed)
        assert normal[:2]==(5,2) and normal[2]==20
        label=b.tuple_label(seed,s,(5,2))
        j=(5-max(label[2:]))//2
        for state in orbit:assert b.tuple_label(state,s,(5,2))==label
        assert c.block_divisors(seed)==[1,2,10] and c.stabilizer_size(seed)==1
        out.append({"phase_smaller_bin":j,"literal_seed":seed,"literal_orbit":len(orbit),
                    "absolute_orbit":len(orbit)//10,"fixed_type_assignment_degree":len(orbit)//120,
                    "group_order":28800,"all_block_count_divisors":[1,2,10],"deck":1,
                    "orbit_sha256":c.digest(orbit)})
    assert len(union)==8400 and sorted(x["literal_orbit"] for x in out)==[2400,6000]
    return {"complete_literal_domain":8400,"components":out}


def negative_controls():
    failures=0
    for args in ((4,2,0,[1],[0,-1,1],F(1),[1]),
                 (5,2,1,[0,1],[0,1],F(1),[1,1]),
                 (5,2,1,[0,1],[-1,1],F(0),[1,1])):
        try:chart(*args)
        except AssertionError:failures+=1
        else:raise AssertionError("invalid Hermite parameters accepted")
    P,Q,H,A=chart(5,2,1,[0,1],[-1,1],F(1,101),[F(-1,2),1])
    assert not admissible(P,H,A,F(1,101))
    boundary=b.component_descriptors(8,4,2,1)
    assert {row["inner_degree"] for row in boundary}=={2,4}
    return {"invalid_inputs_rejected":failures,"symmetric_other_value_collision_rejected":True,
            "excluded_boundary_n4_base_degrees":[row["inner_degree"] for row in boundary]}


def main():
    inputs=["POLYNOMIAL_HERMITE_POWER_PHASE_CHARTS.md",
            "POLYNOMIAL_PERIOD_DIAGONALS_AND_POWER_PHASES.md",
            "POLYNOMIAL_INVOLUTION_COMPONENT_CLASSIFICATION.md",
            "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md",
            "verify_polynomial_involution_components.py","verify_polynomial_lift_transport.py",
            "verify_polynomial_passport_components.py",
            "sources/bogatyrev-gendron-pell-abel-2025.pdf"]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE)) for module in (c.m,c.h,c.v,c.v.one))
    inputs=list(dict.fromkeys(inputs));sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={"schema_version":1,"kind":"Hermite_power_phase_rational_charts","result":"PASS",
        "python":platform.python_version(),"checker_sha256":sha(Path(__file__).resolve()),
        "input_sha256":{p:sha(HERE/p) for p in inputs},
        "exact_chart_witnesses":chart_controls(),"degree10_surface_controls":degree10_surfaces(),
        "size_controls":size_controls(),"degree10_orbits":degree10_orbits(),
        "negative_controls":negative_controls(),
        "scope":"Exact finite Hermite/norm/open-locus certificates and8400-tuple degree10 coverage. Uniform nonemptiness, rationality and classification are the written J1-J3 proofs; no general CMS graph enumerator is implemented.",
        "independent_replay":False,"formal_verification":False,"external_review":False}
    (HERE/"polynomial_phase_interpolation_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: "+str(len(receipt["exact_chart_witnesses"]))+
          " exact Hermite phase charts;8400 complete degree10 tuples")


if __name__=="__main__":main()
