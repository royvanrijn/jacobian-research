"""Exact reproduction controls and arithmetic/block checks for KMFI69--72.

Unbounded conclusions are proved in MARKED_WREATH_ACTIONS.md.
Finite ranks and examples do not establish an action classification.
"""
import os
os.environ.setdefault("SYMPY_GROUND_TYPES","flint")
from pathlib import Path
from itertools import product
from math import factorial,gcd
import hashlib,json,platform,time
import sympy as sp
from sympy.combinatorics import Permutation,PermutationGroup
from sympy.matrices.normalforms import smith_normal_form
from marked_wreath import (polynomial_map,geometric_generators,arithmetic_generators,
 sheet_action,colored_allowed,cyclotomic_twist,generated_block,permutation_cycles)
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
polys=[]
for n in range(3,8):
    variables,T,phi,F=polynomial_map(n);x,y=variables[:2];p=1+x*y;m=n-2
    check(f"linear coefficient exactly1 n{n}",phi.nth(1)==1)
    check(f"all coordinates polynomial overQ n{n}",all(sp.Poly(f,*variables,domain=sp.QQ) is not None for f in F))
    check(f"root identity after exact denominator clearing n{n}",
      sp.expand(sum(phi.nth(j)*x**j*p**(n-j) for j in range(n+1)))==0)
    check(f"derivative marking identity cleared n{n}",
      sp.expand(sum(j*phi.nth(j)*x**(j-1)*p**(n-j) for j in range(1,n+1))-p)==0)
    high=sp.Matrix(F[:-2]).jacobian(variables[2:])
    check(f"triangular high coefficient determinant n{n}",sp.factor(high.det(method="domain-ge"))==p**(2*n-3))
    jac=sp.Matrix(F).jacobian(variables).det(method="domain-ge")
    check(f"entire expanded polynomial Jacobian n{n}",sp.cancel(jac)==sp.Rational(m,2))
    check(f"boundary coordinate slope n{n}",sp.expand(F[-2].subs(x,0)+sp.Rational(m,2)*y)==0)
    polys.append({'n':n,'degree':n*m,'jacobian':str(sp.Rational(m,2)),
      'monomials':sum(len(sp.Poly(f,*variables).terms()) for f in F),'entire_Jacobian_checked':True})
family=[]
for n in range(3,11):
    m=n-2;N=n*m;gens=geometric_generators(n)
    G=PermutationGroup([Permutation(g) for g in gens])
    expected=m**n*factorial(n)//gcd(2,m)
    check(f"exact geometric permutation group order n{n}",G.order()==expected)
    check(f"faithful transitive marked action n{n}",G.degree==N and G.is_transitive())
    H=G.stabilizer(0);orbits=sorted(H.orbits(),key=lambda a:min(a))
    check(f"stabilizer order n{n}",H.order()==expected//N)
    check(f"stabilizer fixes whole root block n{n}",all({j} in orbits for j in range(m)))
    check(f"stabilizer transitive outside root block n{n}",
      frozenset(range(m,N)) in [frozenset(o) for o in orbits])
    diagonal=sheet_action(n,[1]*n)
    check(f"cyclic geometric deck centralizes all generators n{n}",
      all(tuple(diagonal[g[i]] for i in range(N))==tuple(g[diagonal[i]] for i in range(N)) for g in gens))
    check(f"sheet stabilizer normally generates geometric group n{n}",G.normal_closure(H).order()==expected)
    if m>1:
        color=sp.Symbol("color")
        check(f"every primitive nonidentity deck has genuine x-zero pole n{n}",
          sp.rem(1-color,sp.cyclotomic_poly(m,color),domain=sp.QQ)!=0)
    fixed=sum(len(o) for o in orbits if len(o)==1)
    check(f"geometric generic deck exact order n{n}",fixed==m)
    block_sizes=set()
    for j in range(N):
        B=generated_block(gens,0,j);block_sizes.add(len(B))
        expectedB=(frozenset(range(N)) if j>=m else frozenset(range(0,m,gcd(m,j))))
        check(f"principal block exact n{n},sheet{j}",B==expectedB)
    divisors=[d for d in range(1,m+1) if m%d==0]
    check(f"complete block sizes n{n}",block_sizes==set(divisors+[N]))
    sigma=tuple(list(range(1,n))+[0]);v=[0]*n;v[0]=1
    full=sheet_action(n,v,sigma)
    check(f"full degree cycle satisfies group condition n{n}",colored_allowed(n,v,sigma))
    check(f"degree cycle actual group membership n{n}",G.contains(Permutation(full)))
    check(f"cycle-gcd invariant equals degree n{n}",[len(c) for c in permutation_cycles(full)]==[N])
    AG=PermutationGroup([Permutation(g) for g in arithmetic_generators(n)])
    phi_m=int(sp.totient(m))
    check(f"exact rational arithmetic group order n{n}",AG.order()==expected*phi_m)
    AH=AG.stabilizer(0)
    arith_fixed=sum(len(o) for o in AH.orbits() if len(o)==1)
    check(f"rational generic deck exact order n{n}",arith_fixed==gcd(2,m))
    check(f"geometric subgroup normal in arithmetic n{n}",G.is_normal(AG))
    # Incidence matrix of all colliding pairs gives the complete potential
    # relation lattice: Smith factors1,...,1,2.
    rows=[]
    for i in range(n):
        for j in range(i+1,n):rows.append([int(k==i)+int(k==j) for k in range(n)])
    D=smith_normal_form(sp.Matrix(rows),domain=sp.ZZ)
    smith=[abs(int(D[i,i])) for i in range(n)]
    check(f"pair valuation Smith factors n{n}",smith==[1]*(n-1)+[2])
    if m%2==0:
        epsilon=(-1)**(n*(n-1)//2)
        check(f"required square root lies in cyclotomic constants n{n}",epsilon==1 or m%4==0)
        for u in range(1,m):
            if gcd(u,m)==1:
                sigma0=list(range(n))
                if cyclotomic_twist(n,u):sigma0[0],sigma0[1]=sigma0[1],sigma0[0]
                check(f"cyclotomic section satisfies exact parity n{n},u{u}",colored_allowed(n,[0]*n,sigma0,u))
        if n%4==2:
            check(f"naive root-fixed conjugation rejected n{n}",
              not colored_allowed(n,[0]*n,tuple(range(n)),m-1))
    family.append({'n':n,'marking_order':m,'sheets':N,'geometric_order':int(expected),
      'arithmetic_order_over_Q':int(AG.order()),'normal_closure_constants':f'Q(zeta_{m})',
      'geometric_deck_order':m,'rational_deck_order_over_Q':gcd(2,m),
      'proper_block_sizes':divisors,'stable_atomicity':'proved in written general argument',
      'cycle_gcd_invariant':N})
# Generic derivative-product identity in independent root coordinates.
for n in range(3,8):
    roots=list(range(1,n+1));q=sp.Rational(3,7)
    ds=[q*sp.prod(roots[i]-roots[j] for j in range(n) if i!=j) for i in range(n)]
    delta=sp.prod(roots[i]-roots[j] for i in range(n) for j in range(i+1,n))
    check(f"derivative product sign normalization n{n}",
      sp.prod(ds)==(-1)**(n*(n-1)//2)*q**n*delta**2)
manifest=json.loads((ROOT/'sources/marked-wreath-manifest.json').read_text())
for s in manifest['sources']:
    f=ROOT/s['path']
    check("primary source bytes",f.stat().st_size==s['bytes'])
    check("primary source SHA256",hashlib.sha256(f.read_bytes()).hexdigest()==s['sha256'])
model={'claims':['KMFI69','KMFI70','KMFI71','KMFI72'],
 'known_geometric_family_reproduced':True,'source_attribution':'Harish construction; Fonollosa geometric monodromy',
 'polynomial_controls':polys,'group_controls':family,
 'universal_statements':{'n':'all integers n>=3','m':'n-2','dimension':'n',
 'geometric_order':'m^n*n!/gcd(2,m)','constants_over_Q':'Q(zeta_m)',
 'all_proper_blocks':'within one root block, cyclic subgroup sizes d|m',
 'all_upper_intermediate_extensions':'L/K(t,a^d) cyclic of degree d for every d|m',
 'stable_atomicity':True,'every_finite_cyclic_geometric_deck_group_realized':True,
 'every_finite_abelian_geometric_deck_group_realized_by_products':True},
 'natural_PSL2_11_degree12_absolute_occurrence':'excluded by KMFI67',
 'exceptional_PSL2_11_degree11_absolute_occurrence_decided':False,
 'full_action_classification_complete':False}
(ROOT/'marked_wreath_model.json').write_text(json.dumps(model,indent=2)+'\n')
names=['MARKED_WREATH_ACTIONS.md','marked_wreath.py','verify_marked_wreath.py','marked_wreath_model.json',
 'sources/marked-wreath-manifest.json',manifest['sources'][0]['path']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,'python':platform.python_version(),
 'sympy':sp.__version__,'runtime_seconds':round(time.monotonic()-started,3),
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in names},
 'assurance':'Written all-n proof and source-attributed geometric reproduction with exact polynomial, permutation, arithmetic and block controls. No formal verification or independent whole-proof review.',
 'goal_complete':False}
(ROOT/'marked_wreath_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'cyclic_deck_orders':'all positive integers','full_goal_complete':False}))
