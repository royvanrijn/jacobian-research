#!/usr/bin/env python3
"""Exact algebraic/phase controls; not an implementation of the global cell model."""
from collections import deque
from fractions import Fraction as F
from itertools import product
from math import gcd
from pathlib import Path
import hashlib,json,platform
import verify_polynomial_lifting_obstructions as d

if not __debug__:
    raise RuntimeError("period-diagonal checker requires active assertions")
HERE=Path(__file__).resolve().parent
c,m,l=d.c,d.m,d.l


# Exact Q[a,b], with monomial exponent pairs.
def clean(p):return {k:F(v) for k,v in p.items() if v}
def add(p,q):
    ans=dict(p)
    for k,v in q.items():ans[k]=ans.get(k,F(0))+v
    return clean(ans)
def scale(p,x):return clean({k:v*x for k,v in p.items()})
def mul(p,q):
    ans={}
    for (i,j),u in p.items():
        for (k,h),v in q.items():ans[i+k,j+h]=ans.get((i+k,j+h),F(0))+u*v
    return clean(ans)
def power(p,k):
    ans={(0,0):F(1)}
    for _ in range(k):ans=mul(ans,p)
    return ans
ONE={(0,0):F(1)};A={(1,0):F(1)};B={(0,1):F(1)}


def formal_remainder():
    dividend={4:ONE,3:A,2:B}
    divisor={2:scale(ONE,4),1:scale(A,3),0:scale(B,2)}
    while dividend and max(dividend)>=2:
        degree=max(dividend);coef=scale(dividend[degree],F(1,4))
        for j,p in divisor.items():
            at=j+degree-2
            dividend[at]=add(dividend.get(at,{}),scale(mul(coef,p),-1))
            if not dividend[at]:del dividend[at]
    linear=scale(add(scale(power(A,3),9),scale(mul(A,B),-32)),F(1,64))
    constant=add(scale(mul(power(A,2),B),F(3,32)),scale(power(B,2),F(-1,4)))
    assert dividend=={1:linear,0:constant}
    total=add(mul(linear,scale(A,F(-3,4))),scale(constant,2))
    expected=scale(add(add(scale(power(A,4),-27),scale(mul(power(A,2),B),144)),
                       scale(power(B,2),-128)),F(1,256))
    assert total==expected
    return linear,constant,total


# Q(sqrt(3)), no numerical root approximations.
def qa(x):return F(x),F(0)
def qadd(x,y):return x[0]+y[0],x[1]+y[1]
def qscale(x,a):return x[0]*a,x[1]*a
def qmul(x,y):return x[0]*y[0]+3*x[1]*y[1],x[0]*y[1]+x[1]*y[0]
def qpow(x,k):
    out=qa(1)
    for _ in range(k):out=qmul(out,x)
    return out
def evaluate(p,a,b):
    ans=qa(0)
    for (i,j),v in p.items():ans=qadd(ans,qscale(qmul(qpow(a,i),qpow(b,j)),v))
    return ans
def printable(x):return [str(v) for v in x]


def valid_quartic(a,b):
    linear,constant,total=formal_remainder()
    if b==qa(0) or qadd(qpow(a,2),qscale(b,-4))==qa(0):return False
    if qadd(qscale(qpow(a,2),9),qscale(b,-32))==qa(0):return False
    U,V=evaluate(linear,a,b),evaluate(constant,a,b)
    prod=qadd(qadd(qscale(qmul(qpow(U,2),b),F(1,2)),
                  qscale(qmul(qmul(a,U),V),F(-3,4))),qpow(V,2))
    return prod!=qa(0) and (U==qa(0) or evaluate(total,a,b)==qa(0))


def quartic_controls():
    linear,constant,total=formal_remainder()
    rows=[]
    for sign in (-1,1):
        a=qa(1);b=F(9,16),sign*F(3,16)
        assert valid_quartic(a,b)
        assert qadd(qadd(qscale(qpow(b,2),128),qscale(b,-144)),qa(27))==qa(0)
        assert b!=qa(0) and qadd(qa(1),qscale(b,-4))!=qa(0)
        critical_disc=qadd(qa(9),qscale(b,-32))
        assert critical_disc!=qa(0)
        U,V=evaluate(linear,a,b),evaluate(constant,a,b)
        assert U!=qa(0) and evaluate(total,a,b)==qa(0)
        # Product of critical values, from r+t=-3a/4 and rt=b/2.
        prod=qadd(qadd(qscale(qmul(qpow(U,2),b),F(1,2)),
                      qscale(qmul(U,V),F(-3,4))),qpow(V,2))
        assert prod!=qa(0)
        rows.append({"b_over_a_squared":printable(b),
                     "critical_quadratic_discriminant":printable(critical_disc),
                     "common_f_critical_value":printable(qscale(prod,-1)),
                     "H_critical_value_ratio":-1,"phase_profile":["T","T"],
                     "deck":1})
    assert rows[0]["b_over_a_squared"]!=rows[1]["b_over_a_squared"]
    assert evaluate(linear,qa(0),qa(1))==qa(0)
    assert evaluate(constant,qa(0),qa(1))==qa(F(-1,4))
    assert evaluate(total,qa(0),qa(1))==qa(F(-1,2))
    assert valid_quartic(qa(0),qa(1))
    bad=((qa(1),qa(0)),(qa(1),qa(F(1,4))),
         (qa(1),qa(F(9,32))),(qa(1),qa(1)))
    assert all(not valid_quartic(a,b) for a,b in bad)
    return {"formal_remainder_identity":True,"formal_critical_sum_identity":True,
            "quadratic_pair":rows,
            "invalid_ramification_or_wrong_grouping_rejected":len(bad),
            "even_quartic":{"a":0,"b":1,"H_critical_value":str(F(-1,4)),
                            "f_critical_value":str(F(1,16)),
                            "phase_profile":["DT","I"],"deck":2}}


def phase_profile(state):
    root=next(i for i,g in enumerate(state) if gcd(*m.cycle_lengths(g))>1)
    D=gcd(*m.cycle_lengths(state[root]))
    wreath=l.to_wreath(state,D)
    assert sorted(c.cycles(wreath[root][0]))==[tuple(range(D))]
    profiles=[]
    for i,(q,maps) in enumerate(wreath):
        if i!=root:
            assert q==tuple(range(D))
            lengths=tuple(tuple(m.cycle_lengths(g)) for g in maps)
            profiles.append(min(lengths[k:]+lengths[:k] for k in range(D)))
    return D,tuple(sorted(profiles))


def degree8_controls():
    old=json.loads((HERE/"polynomial_lift_transport_verification.json").read_text())
    rows=old["block_size_four_obstruction"]["complete_full_passport_components"]
    seen=set();out=[];groups=[]
    for row in rows:
        seed=tuple(tuple(g) for g in row["literal_seed"])
        orbit,_=c.closure(seed)
        assert not seen&orbit and len(orbit)==row["literal_orbit"]
        seen.update(orbit)
        profile=phase_profile(seed)
        for state in orbit:assert phase_profile(state)==profile
        G=d.r.perm_group(seed,8)
        assert len(G)==row["group_order"]
        deck=c.stabilizer_size(seed)
        assert deck==row["cyclic_stabilizer_order"]
        value=None
        if deck==1:
            assert profile==(2,(((1,1,2),(1,1,2)),))
            value=d.cubic_h8(seed)
            for state in orbit:assert d.cubic_h8(state)==value
            groups.append(G)
        else:
            assert profile==(2,(((1,1,1,1),(2,2)),))
        out.append({"literal_seed":seed,"literal_orbit":len(orbit),"deck":deck,
                    "group_order":len(G),"power_phase_return_types":profile,
                    "cubic_value":value,"orbit_sha256":c.digest(orbit)})
    assert groups[0]==groups[1]
    assert sorted(x["cubic_value"] for x in out if x["cubic_value"] is not None)==[1,2]
    assert len(seen)==c.character_total([[2,2,4],[1,1,1,1,2,2]])*2==40
    return {"complete_literal_domain":40,"same_actual_H576":True,"components":out}


def dihedral_polynomial_control():
    # (z^3-3z)^2 = 2 T_6(z/2)+2; critical values H(-1)=2,H(1)=-2.
    H=[0,-3,0,1]
    square=[sum(H[j]*H[k] for j in range(4) for k in range(4) if j+k==i)
            for i in range(7)]
    assert square==[0,0,9,0,-6,0,1]
    assert sum(H[j]*(-1)**j for j in range(4))==2
    assert sum(H[j] for j in range(4))==-2
    R=tuple((1-x)%6 for x in range(6))
    S=tuple((-x)%6 for x in range(6))
    seed=(R,S);c.validate(seed)
    orbit,_=c.closure(seed);G=d.r.perm_group(seed,6)
    assert len(G)==12 and len(orbit)==6 and c.stabilizer_size(seed)==2
    assert phase_profile(seed)==(2,(((1,2),(1,2)),))
    return {"polynomial_coefficients":square,"H_critical_values":[2,-2],
            "f_critical_value":4,"relative_period_mod2":1,
            "literal_seed":seed,"group_order":12,"literal_orbit":6,"deck":2}


def slice_controls(L):
    # Abstract single period chart: 0<x<L, y arbitrary. This does not generate CMS cells.
    xs=[("interval",i,i+1) for i in range(L)]+[("point",i,i) for i in range(1,L)]
    ys=[("below",),("point",),("above",)]
    cells=set(product(range(len(xs)),range(3)))
    def is_face(low,high):
        if low==high:return True
        if xs[low[0]][0]!="point" or xs[high[0]][0]!="interval":return False
        return xs[low[0]][1] in xs[high[0]][1:]
    def edge(a,b):
        xface=(a[0]==b[0] or is_face(a,b) or is_face(b,a))
        yface=(a[1]==b[1] or a[1]==1 or b[1]==1)
        # To be incident the dimension comparisons must point the same way.
        dimx=lambda k:int(xs[k][0]=="interval")
        dimy=lambda k:int(k!=1)
        compatible=((dimx(a[0])<=dimx(b[0]) and dimy(a[1])<=dimy(b[1])) or
                    (dimx(a[0])>=dimx(b[0]) and dimy(a[1])>=dimy(b[1])))
        return a!=b and xface and yface and compatible
    counts=[]
    for equal in (True,False):
        keep={x for x in cells if (xs[x[0]][0]=="point" and x[1]==1)==equal}
        remaining=set(keep);components=0
        while remaining:
            todo=[remaining.pop()];components+=1
            while todo:
                x=todo.pop()
                adjacent={y for y in remaining if edge(x,y)}
                remaining-=adjacent;todo.extend(adjacent)
        counts.append(components)
    assert counts==[L-1,1]
    return {"strip_width":L,"integer_diagonal_components":counts[0],
            "complement_components":counts[1],
            "excluded_endpoints_and_deleted_points_create_no_false_edges":True}


def cubic_ambient_controls():
    old=json.loads((HERE/"polynomial_lifting_obstructions_verification.json").read_text())
    rows=[]
    for row in old["cyclic_wreath_cubic_family"]["complete_lift_orbits"]:
        if row["quotient_group_order"]!=2520:continue
        b=row["cyclic_fibre"];seed=tuple(tuple(g) for g in row["literal_seed"])
        N=7*b;parts=[1]*(3*b)+[2]*(2*b)
        root=next(i for i,g in enumerate(seed) if m.cycle_lengths(g)==parts)
        assert gcd(*parts)==1
        other=[g for i,g in enumerate(seed) if i!=root]
        zero_orders=sorted(len(C)-1 for g in other for C in c.cycles(g,False))
        assert zero_orders==sorted([1]*(4*b)+[b-1])
        codimension=(2*b-1)*2
        assert len(zero_orders)-1-codimension==2
        quotient=tuple(q for q,_ in l.to_wreath(seed,7))
        assert d.cubic_a7(quotient)==row["cubic_value"]
        rows.append({"cyclic_fibre":b,"ambient_zero_orders":zero_orders,
                     "ambient_component_count_from_primary_theorem":1,
                     "diagonal_codimension":codimension,"diagonal_dimension":2,
                     "literal_component_size_from_pinned_D6":row["literal_orbit"],
                     "cubic_value_from_pinned_D6":row["cubic_value"]})
    assert len(rows)==4
    return {"finite_inherited_records":rows,
            "scope":"D6's pinned written theorem gives every b>=2; these records check b2,b3 dimensions and seed quotient cubic values, not new full lift replays."}


def main():
    # Negative field/ramification controls: b=0 has the wrong root fiber, and
    # b=9/32,a=1 collides the two remaining critical points.
    assert evaluate(formal_remainder()[0],qa(1),qa(F(9,32)))==qa(0)
    old=json.loads((HERE/"polynomial_lifting_obstructions_verification.json").read_text())
    inputs=[
      "POLYNOMIAL_PERIOD_DIAGONALS_AND_POWER_PHASES.md",
      "POLYNOMIAL_SINGLE_CYCLE_TAIL_COMPONENTS.md","POLYNOMIAL_PASSPORT_COMPONENTS.md",
      "POLYNOMIAL_PASSPORT_LIFT_TRANSPORT.md","POLYNOMIAL_INTRINSIC_CLASSES_AND_LIFTING.md",
      "verify_polynomial_lifting_obstructions.py","polynomial_lift_transport_verification.json",
      "polynomial_lifting_obstructions_verification.json",
      "sources/chen-mullane-salter-tahar-residual-stratified-2026.pdf"]
    inputs.extend(old["input_sha256"])
    inputs=list(dict.fromkeys(inputs))
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={"schema_version":1,"kind":"period_diagonals_power_phases_quartic_cubic_split",
      "result":"PASS","python":platform.python_version(),
      "checker_sha256":sha(Path(__file__).resolve()),
      "input_sha256":{p:sha(HERE/p) for p in inputs},
      "quartic_algebraic_classification":quartic_controls(),
      "complete_degree8_control":degree8_controls(),
      "dihedral_power_phase_example":dihedral_polynomial_control(),
      "abstract_rational_slice_controls":[slice_controls(L) for L in (2,3,5)],
      "connected_ambient_cubic_diagonals":cubic_ambient_controls(),
      "global_cell_graph_implemented":False,
      "global_scope":"Written G1-G3 reduction/width/power arguments with explicitly imported CMS2026 period model; no full classification follows from these finite controls.",
      "independent_replay":False,"formal_verification":False,"external_review":False}
    (HERE/"polynomial_period_diagonals_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: exact quartic Q(sqrt3) classification;40 degree8 tuples; power-phase and rational slice controls")


if __name__=="__main__":
    main()
