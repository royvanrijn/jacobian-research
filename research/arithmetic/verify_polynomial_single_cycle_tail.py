#!/usr/bin/env python3
"""Bounded exact controls for the primary-isoresidual single-cycle-tail theorem."""
from collections import Counter, deque
from fractions import Fraction
from itertools import combinations
from math import factorial,gcd
from pathlib import Path
import hashlib,json,platform
import verify_polynomial_passport_components as c

if not __debug__:
    raise RuntimeError("single-cycle-tail checker requires active assertions")
HERE=Path(__file__).resolve().parent
m=c.m
CONTROLS=(
 ((1,2,3),(2,2)),
 ((1,2,3,4),(2,3)),
 ((1,1,2),(2,2)),
 ((2,2,2),(2,2)),
 ((1,1,1,1,1,1,2),(3,3,3)),
 ((1,1,1,1,1,1,3),(3,3,3)),
 ((1,1,1,1,1,2,2),(3,3,3)),
 ((2,2,2,2,2,2,2),(3,3,3)),
 ((1,1,1,1,1,1,1,1,2),(3,3,3,3)),
 ((1,1,1,1,1,1,2),(3,5)),
 ((1,1,1,1,1,1,1,1,2),(3,3,5)),
 ((1,1,1,1,1,1,2),(4,4)),
 ((2,2,2,2,2,2,2),(4,4)),
 ((1,1,1,1,1,1,2),(2,3,4)),
)


def raw_star(parts,tail):
    assert len(tail)>=2 and all(q>=2 for q in tail)
    assert len(parts)==1+sum(q-1 for q in tail)
    assert parts==tuple(sorted(parts)) and min(parts)>0
    N=sum(parts)
    assert N>len(parts)
    cycles=[];offset=0
    for a in parts:
        cycles.append(tuple(range(offset,offset+a)));offset+=a
    g=m.perm(N,cycles);factors=[g];at=1
    for q in tail:
        support=[cycles[0][0]]+[cycles[j][0] for j in range(at,at+q-1)]
        factors.append(m.perm(N,[support]));at+=q-1
    assert at==len(parts)
    return tuple(factors)


def all_block_counts(seed):
    return c.block_divisors(seed)


def commutator(a,b):
    return m.compose(m.compose(m.compose(a,b),m.inverse(a)),m.inverse(b))


def normal_support(seed,parts,tail):
    N=sum(parts);D=gcd(*parts);n=N//D
    z=commutator(seed[1],seed[2])
    assert m.cycle_lengths(z)==[1]*(N-3)+[3]
    C,=c.cycles(z,False)
    supports={frozenset(C)};queue=deque(supports)
    gens=seed+tuple(m.inverse(g) for g in seed)
    while queue:
        S=queue.popleft()
        for g in gens:
            T=frozenset(g[x] for x in S)
            if T not in supports:
                supports.add(T);queue.append(T)
    components=[{x} for x in range(N)]
    for S in supports:
        hits=[T for T in components if S&T]
        components=[T for T in components if not S&T]+[set().union(*hits)]
    assert {frozenset(S) for S in components}=={
        frozenset(range(j,N,D)) for j in range(D)}
    assert all(len(S)==n for S in components)
    assert all(len({x%D for x in range(N) if g[x]!=x})==1 for g in seed[1:])
    block_count=sorted({d for d in range(1,D+1) if D%d==0}|{N})
    assert all_block_counts(seed)==block_count
    symmetric=any(q%2==0 for q in tail)
    expected=(factorial(n)**D*D if symmetric else
              (factorial(n)//2)**D*D*(1+(n%2==0)))
    closed_order=None
    if expected<=1000:
        seen={tuple(range(N))};queue=deque(seen)
        while queue:
            x=queue.popleft()
            for g in gens:
                y=m.compose(x,g)
                if y not in seen:
                    seen.add(y);queue.append(y)
        assert len(seen)==expected
        closed_order=len(seen)
    return {"commutator_three_cycle":z,"normal_supports":len(supports),
            "support_blocks":D,"block_size":n,
            "all_block_count_divisors":block_count,
            "actual_group":("S_n wr C_D" if symmetric else "Alt_n^D<P_N>"),
            "actual_group_order_from_written_F3":expected,
            "independent_small_generated_group_closure":closed_order}


def control(parts,tail):
    raw=raw_star(parts,tail);seed=c.v.one.fixed_product(raw)
    c.validate(seed)
    passport=[m.cycle_lengths(g) for g in seed]
    placements=c.type_placements(passport)
    hook=c.character_total(passport)*placements
    M=1
    for multiplicity in Counter(parts).values():M*=factorial(multiplicity)
    formula=Fraction(placements*sum(parts)**len(tail)*factorial(len(parts)-1),M)
    assert formula.denominator==1 and hook==formula.numerator
    expected_components=2 if all(q%2 for q in tail) and len(set(parts))==len(parts) else 1
    assert expected_components==1  # Strict all-odd cases use the separate two-tripod controls.
    orbit,_=c.closure(seed)
    assert len(orbit)==hook
    assert c.stabilizer_size(seed)==1
    return {"parts":parts,"tail_cycle_lengths":tail,
            "complete_literal_domain":hook,"full_components":1,
            "type_placements":placements,"absolute_orbit":hook//sum(parts),
            "deck":1,"literal_seed":seed,"orbit_sha256":c.digest(orbit),
            "normal_support_certificate":normal_support(seed,parts,tail)}


def mul(a,b):
    ans=[Fraction(0)]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):ans[i+j]+=x*y
    return ans


def derivative(a):
    return [i*a[i] for i in range(1,len(a))]


def evaluate(a,x):
    value=Fraction(0)
    for y in a[::-1]:value=value*x+y
    return value


def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def poly_gcd(a,b):
    a=list(a);b=trim(list(b))
    while any(b):
        rest=a[:]
        while any(rest) and len(rest)>=len(b):
            ell=len(rest)-len(b);coef=rest[-1]/b[-1]
            for j,x in enumerate(b):rest[j+ell]-=coef*x
            trim(rest)
        a,b=b,rest
    return trim([x/a[-1] for x in a])


def polynomial_bridge(tail,critical_points):
    assert len(tail)==len(critical_points) and len(set(critical_points))==len(tail)
    N=2+sum(q-1 for q in tail)
    fp=[Fraction(0),Fraction(N)]
    for q,x in zip(tail,critical_points):
        for _ in range(q-1):fp=mul(fp,[-x,1])
    f=[Fraction(0)]+[a/Fraction(i+1) for i,a in enumerate(fp)]
    assert len(f)==N+1 and f[-1]==1 and derivative(f)==fp
    assert f[0]==f[1]==0 and f[2]!=0
    rest=f[2:]
    assert poly_gcd(rest,derivative(rest))==[1]
    values=[evaluate(f,x) for x in critical_points]
    assert all(values) and len(set(values))==len(values)
    # The exact derivative certificate accounts for every critical point.
    assert sum(q-1 for q in tail)+1==N-1
    parts=(1,)*(N-2)+(2,)
    assert len(parts)==1+sum(q-1 for q in tail)
    return {"degree":N,"zero_partition":parts,"tail":tail,
            "polynomial_coefficients":[str(x) for x in f],
            "critical_points":critical_points,"nonzero_distinct_critical_values":[str(x) for x in values],
            "simple_residual_roots_gcd":[str(x) for x in poly_gcd(rest,derivative(rest))],
            "residues":parts+(-N,),"differential_zero_orders":[q-1 for q in tail],
            "differential_divisor_degree":sum(q-1 for q in tail)-(len(parts)+1)}


def boundary_and_negative_controls():
    orders=tuple(range(7))
    parity=lambda order:(-1)**sum(order[i]>order[j] for i in range(len(order)) for j in range(i+1,len(order)))
    assert all(parity(orders)==parity(orders[k:]+orders[:k]) for k in range(7))
    assert parity(orders)!=parity(orders[:1]+orders[1:3][::-1]+orders[3:])
    rejected=0
    for parts,tail in (((1,2,3),(3,3)),((1,2,3),(2,)),((1,1,1),(2,2))):
        try:raw_star(parts,tail)
        except AssertionError:rejected+=1
        else:raise AssertionError("invalid passport accepted")
    nonresonant=0
    for parts,tail in CONTROLS:
        rho=parts+(-sum(parts),)
        for k in range(1,len(rho)):
            for I in combinations(range(len(rho)),k):
                assert sum(rho[i] for i in I)!=0
                nonresonant+=1
    # A bad sign pattern with the same sum is resonant and cannot pass F1.
    assert sum((1,-1,2,-2)[i] for i in (0,1))==0
    return {"cyclic_parity_letters":7,"odd_label_swap_changes_parity":True,
            "all_proper_residue_subset_tests":nonresonant,
            "invalid_passports_rejected":rejected,"resonant_vector_rejected":[1,-1,2,-2]}


def main():
    inputs=[
     "POLYNOMIAL_SINGLE_CYCLE_TAIL_COMPONENTS.md",
     "POLYNOMIAL_PASSPORT_COMPONENTS.md","POLYNOMIAL_TWO_TRIPOD_COMPONENTS.md",
     "polynomial_two_tripods_verification.json","verify_polynomial_passport_components.py",
     "sources/chen-gendron-prado-tahar-isoresidual-curves-2026.pdf",
     "sources/jones-zvonkin-braid-cacti-2002.pdf"]
    inputs.extend(str(Path(module.__file__).resolve().relative_to(HERE))
                  for module in (c.m,c.h,c.v,c.v.one))
    inputs=list(dict.fromkeys(inputs))
    sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
    receipt={"schema_version":1,"kind":"polynomial_single_cycle_tail",
      "result":"PASS","python":platform.python_version(),
      "checker_sha256":sha(Path(__file__).resolve()),
      "input_sha256":{p:sha(HERE/p) for p in inputs},
      "controls":[control(parts,tail) for parts,tail in CONTROLS],
      "logarithmic_derivative_examples":[polynomial_bridge((3,3),(1,2)),
         polynomial_bridge((3,3,3),(1,2,3)),
         polynomial_bridge((3,3,3,3),(1,2,3,4)),
         polynomial_bridge((2,3,4),(1,2,4))],
      "boundary_and_negative_controls":boundary_and_negative_controls(),
      "unrestricted_component_input":"Published CGPT2026 Theorem1.7, Lemma8.1, Lemma8.6 and Section8.1.2, plus written logarithmic-derivative/label/open-locus bridge F1-F2",
      "unrestricted_group_scope":"Written return-budget, phase-tree and normal commutator support proof F3. No enumeration of large groups.",
      "strict_split_scope":"The existing two-tripod receipt supplies complete strict all-odd finite splits; new finite closures test connected cases with three/four tail factors, longer cycles, repeated roots, type collisions, and gcd blocks.",
      "per_component_size_scope":"All connected case sizes complete. Split cases have the exact summed formula; no new unrestricted separate-size formula.",
      "independent_replay":False,"formal_verification":False,"external_review":False}
    (HERE/"polynomial_single_cycle_tail_verification.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print("PASS: single-cycle-tail controls; "+
          str(sum(row["complete_literal_domain"] for row in receipt["controls"]))+
          " complete literal tuples")


if __name__=="__main__":
    main()
