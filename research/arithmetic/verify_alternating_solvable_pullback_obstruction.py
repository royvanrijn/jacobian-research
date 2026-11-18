#!/usr/bin/env python3
"""Exact branch and group certificates for the solvable-pullback obstruction.

The all-rational-map and density assertions are proved in the companion
note using the named Chebotarev/HIT inputs. No bounded check replaces them.
"""
from __future__ import annotations
import argparse,hashlib,json,platform
from itertools import combinations
from math import comb
from pathlib import Path
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
INPUT=HERE/'alternating_nonzero_branch_residual_inputs.json'
A,z,T=s.symbols('A z T')

def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()
def trim(f):
    f=list(f)
    while f and f[-1]==0:f.pop()
    return f
def sub(f,g,p):
    return trim([((f[i]if i<len(f)else 0)-(g[i]if i<len(g)else 0))%p
                 for i in range(max(len(f),len(g)))])
def mul(f,g,p):
    out=[0]*(max(0,len(f)+len(g)-1))
    for i,a in enumerate(f):
        for j,b in enumerate(g):out[i+j]=(out[i+j]+a*b)%p
    return trim(out)
def rem(f,g,p):
    f=trim([int(a)%p for a in f]);g=trim([int(a)%p for a in g])
    if not g:raise ZeroDivisionError
    inverse=pow(g[-1],-1,p)
    while len(f)>=len(g):
        offset=len(f)-len(g);co=f[-1]*inverse%p
        for j,b in enumerate(g):f[offset+j]=(f[offset+j]-co*b)%p
        f=trim(f)
    return f
def powmod(f,n,g,p):
    out=[1];f=rem(f,g,p)
    while n:
        if n&1:out=rem(mul(out,f,p),g,p)
        f=rem(mul(f,f,p),g,p);n//=2
    return out
def gcd(f,g,p):
    f=trim(f);g=trim(g)
    while g:f,g=g,rem(f,g,p)
    return [(a*pow(f[-1],-1,p))%p for a in f]if f else []
def irreducible(desc,p):
    f=trim(list(reversed(desc)));n=len(f)-1
    assert n>=1 and f[-1]==1
    x=[0,1]
    if rem(sub(powmod(x,p**n,f,p),x,p),f,p):return False
    return all(gcd(f,sub(powmod(x,p**(n//q),f,p),x,p),p)==[1]
               for q in s.factorint(n))
def md(value,p):
    value=s.Rational(value)
    return int(value.p)*pow(int(value.q),-1,p)%p
def compose(a,b):return tuple(a[b[i]]for i in range(len(a)))
def perm(n,*cycles):
    out=list(range(n))
    for cycle in cycles:
        for i,x in enumerate(cycle):out[x]=cycle[(i+1)%len(cycle)]
    return tuple(out)
def closure(gens):
    seen={tuple(range(len(gens[0])))};todo=list(seen)
    while todo:
        x=todo.pop()
        for g in gens:
            y=compose(x,g)
            if y not in seen:seen.add(y);todo.append(y)
    return seen

def run():
    d=json.loads(INPUT.read_text())
    N=9*z**10-864*z**8+28161*z**6-350150*z**4+1127500*z**2+250000
    q=z*z+A*z-s.Rational(10,3)
    J=A**4-s.Rational(109,3)*A*A+s.Rational(1690,9)
    alpha=(-29800*A**3+1453400*A)/3
    G=(z**6-2*A*z**5+(3*A*A-s.Rational(268,3))*z**4
       +(-4*A**3+172*A)*z**3+(s.Rational(14251,9)-s.Rational(199,3)*A*A)*z*z
       +(s.Rational(298,3)*A**3-s.Rational(10034,3)*A)*z+2500)
    assert s.expand(s.sympify(d['N'])-N)==0
    assert s.expand(s.sympify(d['modulus'])-J)==0
    assert s.expand(s.sympify(d['s'])-alpha)==0
    assert d['q_coefficients']==list(map(str,s.Poly(q,z).all_coeffs()))
    assert len(d['G_coefficients'])==7
    assert all(s.expand(s.sympify(a)-b)==0 for a,b in zip(d['G_coefficients'],s.Poly(G,z).all_coeffs()))
    reduce=lambda expr:s.expand(s.rem(expr,J,A))
    residual=s.Poly(s.expand(N-alpha*z-9*q*q*G),z)
    assert all(reduce(co)==0 for co in residual.all_coeffs())
    # An independent division reconstructs alpha and the minimal a equation.
    quotient,remainder=s.div(N-T*z,q*q,z)
    solved=s.solve(s.Poly(remainder,z).nth(1),T)[0]
    assert reduce(solved-alpha)==0
    coeffs=s.Poly(s.expand(remainder.subs(T,solved)),z).all_coeffs()
    assert s.Poly(s.gcd_list(coeffs),A).monic().as_expr()==J
    assert all(reduce(a-b)==0 for a,b in zip(s.Poly(quotient/9,z).all_coeffs(),s.Poly(G,z).all_coeffs()))

    control=d['quartic_irreducibility_control'];p=control['p']
    assert p==11 and control['coefficients']==[md(co,p)for co in s.Poly(J,A).all_coeffs()]
    assert irreducible(control['coefficients'],p)
    jdisc=s.discriminant(J,A);assert jdisc==s.Rational(8754497440,9)
    basis=s.Matrix([[s.Poly(reduce(alpha**i),A).nth(j)for i in range(4)]for j in range(4)])
    determinant=basis.det();assert determinant!=0 and str(determinant)==d['alpha_power_basis_determinant']
    inverse=basis.inv()*s.Matrix([0,1,0,0])
    assert list(map(str,inverse))==d['A_from_alpha_coefficients']
    M=s.Matrix([[s.Poly(reduce(alpha*A**i),A).nth(j)for i in range(4)]for j in range(4)])
    W=M.charpoly(T).as_expr()
    expected=T**4-s.Rational(522379320400000,243)*T*T+s.Rational(68130645548641000000000000000,59049)
    assert s.expand(W-expected)==0 and s.discriminant(W,T)!=0
    assert reduce(s.discriminant(q,z))!=0
    assert reduce(s.resultant(q,G,z))!=0

    frobenius=[]
    for row in d['frobenius_controls']:
        p,a=row['p'],row['A'];assert s.isprime(p)
        assert md(J.subs(A,a),p)==0 and md(jdisc,p)!=0
        expected_g=[md(co.subs(A,a),p)for co in s.Poly(G,z).all_coeffs()]
        assert expected_g==row['G_coefficients']
        product=[1];degrees=[]
        for factor in row['factors']:
            coefficients=factor['coefficients'];multiplicity=factor['multiplicity']
            assert multiplicity==1 and irreducible(coefficients,p)
            product=mul(product,list(reversed(coefficients)),p);degrees.append(len(coefficients)-1)
        assert list(reversed(product))==expected_g
        derivative=[(i*product[i])%p for i in range(1,len(product))]
        assert gcd(product,derivative,p)==[1]
        assert degrees==row['factor_degrees']
        frobenius.append(dict(p=p,a=a,degrees=degrees))
    assert [row['degrees']for row in frobenius[:3]]==[[1,5],[6],[1,1,1,1,2]]
    witness=d['square_pullback_witness'];p=witness['p'];modulus=p*p
    assert p==347 and md(J.subs(A,witness['a_lift_mod_p2']),modulus)==0
    assert witness['a_lift_mod_p2']%p==witness['a']==337
    assert md(alpha.subs(A,witness['a_lift_mod_p2']),modulus)==witness['alpha_mod_p2']
    assert witness['alpha_mod_p2']%p==witness['alpha_mod_p']==143
    assert witness['beta_lift_mod_p2']%p==witness['beta_mod_p']==166
    assert (witness['beta_lift_mod_p2']**2-witness['alpha_mod_p2'])%modulus==0
    assert (witness['test_t']**2-witness['alpha_mod_p2'])%modulus==witness['contact_mod_p2']
    assert witness['contact_mod_p2']%p==0 and witness['contact_mod_p2']!=0
    assert [x for x in range(p)if md(q.subs({A:337,z:x}),p)==0]==witness['critical_roots_mod_p']==[17,340]
    # Explicit tame local pattern, checked against every possible 5+5 partition.
    inertia=perm(10,(0,1),(2,3));h=perm(10,(4,5),(6,7,8,9))
    local=closure([inertia,h]);assert len(local)==8
    orbits=[];unused=set(range(10))
    while unused:
        point=min(unused);orbit={g[point]for g in local};orbits.append(sorted(orbit));unused-=orbit
    assert sorted(map(len,orbits))==[2,2,2,4]
    partitions=[];letters=set(range(10))
    for half in combinations(range(10),5):
        left=set(half);right=letters-left
        if 0 not in left:continue
        if all({g[x]for x in left}in(left,right)for g in (inertia,h)):partitions.append(half)
    assert not partitions
    a=perm(6,(0,1),(3,4));b=perm(6,(1,2),(3,4))
    assert compose(compose(compose(a,b),a),b)==perm(6,(0,2,1))
    uniform_controls=[]
    for n in (10,14,18,26):
        b=n//2
        rank=sum(comb(n,k)for k in range(1,b,2))+comb(n,b)//2
        assert rank==2**(n-2)
        inertia_n=perm(n,(0,1),(2,3));h_n=perm(n,(4,5),tuple(range(6,n)))
        local_n=closure([inertia_n,h_n]);assert len(local_n)==2*(n-6)
        sizes=[];unused=set(range(n))
        while unused:
            point=min(unused);orbit={g[point]for g in local_n};sizes.append(len(orbit));unused-=orbit
        assert sorted(sizes)==[2,2,2,n-6]
        uniform_controls.append(dict(degree=n,element_cover_rank=rank,local_group_order=len(local_n),local_orbits=sorted(sizes)))
    report=dict(result='PASS',schema_version=1,python=platform.python_version(),sympy=s.__version__,
        branch_field_degree=4,branch_parameter_modulus=str(J),nonzero_branch_polynomial=str(W),
        field_generation_determinant=str(determinant),exact_factorization_P_alpha_9q_squared_G=True,
        relative_residual_group='S6',frobenius_certificates=frobenius,
        local_group='C2 x C4',local_group_order=len(local),local_orbits=orbits,
        fixed_unordered_5_5_partitions=partitions,perfect_subgroup_commutator_verified=True,
        square_parameter_local_witness=witness,
        conditional_uniform_transfer_controls=uniform_controls,
        uniform_scope_supplied_by_written_proof=True,all_rational_maps_tested_by_sampling=False,
        assurance=dict(independent_replay=False,formal_verification=False,external_review=False))
    paths=[INPUT,HERE/'NONZERO_BRANCH_RESIDUAL_S6_AND_SOLVABLE_PULLBACKS.md',Path(__file__).resolve(),
        HERE/'RATIONAL_DEGREE_TEN_ALTERNATING_PENCIL.md',HERE/'ALTERNATING_SUBGROUP_COVERS_AND_HASSE.md',
        HERE/'ALTERNATING_PENCIL_RAMIFICATION_DENSITY.md',HERE/'sources/lenstra-chebotarev-density-notes.pdf',
        HERE/'sources/zywina-quantitative-hit-2010.pdf']
    report['input_sha256']={path.relative_to(ROOT).as_posix():sha(path)for path in paths}
    report['checker_sha256']=sha(Path(__file__))
    return report

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--write-receipt',type=Path);args=parser.parse_args()
    payload=json.dumps(run(),indent=2,sort_keys=True)+'\n'
    if args.write_receipt:args.write_receipt.write_text(payload)
    print(payload,end='')
if __name__=='__main__':main()
