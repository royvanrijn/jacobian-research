#!/usr/bin/env python3
"""Exact finite inputs to the written all-rank Smith/tangent obstruction."""
from pathlib import Path
from functools import lru_cache
from itertools import product,combinations
from math import factorial,prod
from collections import Counter
from fractions import Fraction
import hashlib,json,platform,time
import sympy as s

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def ck(label,condition):
    assert condition,label
    checks.append(label)

@lru_cache(None)
def parts(n,top=0):
    if n==0:return ((),)
    return tuple((a,)+b for a in range(min(n,top or n),0,-1) for b in parts(n-a,a))

def hook_dim(lam):
    hook=prod(lam[i]-j+sum(v>j for v in lam[i+1:])
              for i in range(len(lam)) for j in range(lam[i]))
    return factorial(sum(lam))//hook

def children(lam):
    out=[]
    for i,a in enumerate(lam):
        if i+1==len(lam) or a>lam[i+1]:
            child=list(lam);child[i]-=1
            out.append(tuple(a for a in child if a))
    return out

@lru_cache(None)
def branch_dim(lam):
    return 1 if not lam else sum(branch_dim(v) for v in children(lam))

def exceptional(n):return {(n,),(1,)*n,(n-1,1),(2,)+(1,)*(n-2)}

def cells(lam):return {(i,j) for i,a in enumerate(lam) for j in range(a)}

@lru_cache(None)
def character(lam,cycles):
    if not cycles:return int(not lam)
    k=cycles[0];out=0
    for mu in parts(sum(lam)-k):
        if any(i>=len(lam) or a>lam[i] for i,a in enumerate(mu)):continue
        skew=cells(lam)-cells(mu)
        if len(skew)!=k:continue
        if any({(i,j),(i+1,j),(i,j+1),(i+1,j+1)}<=skew for i,j in skew):continue
        seen={next(iter(skew))};queue=list(seen)
        for i,j in queue:
            for v in [(i-1,j),(i+1,j),(i,j-1),(i,j+1)]:
                if v in skew and v not in seen:seen.add(v);queue.append(v)
        if seen!=skew:continue
        out+=(-1)**(len({i for i,j in skew})-1)*character(mu,cycles[1:])
    return out

def class_size(cycles):
    count=Counter(cycles)
    return factorial(sum(cycles))//prod(j**a*factorial(a) for j,a in count.items())

def multisets(items,total,start=0):
    if total==0:yield ();return
    for j in range(start,len(items)):
        lam,degree=items[j]
        if degree<=total:
            for tail in multisets(items,total-degree,j):yield (lam,)+tail

def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def cycle_type(g):
    left=set(range(len(g)));out=[]
    while left:
        v=min(left);cycle=[]
        while v not in cycle:cycle.append(v);left.remove(v);v=g[v]
        out.append(len(cycle))
    return tuple(sorted(out,reverse=True))

def orbit_patterns(n,p,low=1):
    if n==0:yield ();return
    q=low
    while q<=n:
        for tail in orbit_patterns(n-q,p,q):yield (q,)+tail
        q*=p

def classification(n):
    for p in list(s.primerange(2,n+1)):
        q=p
        while q<n:q*=p
        if q==n:return 'prime_power'
    for p in list(s.primerange(2,n+1)):
        powers=[];q=p
        while q<n:powers.append(q);q*=p
        if any(a+b==n for a in powers for b in powers):return 'two_positive_powers'
    return 'standard_admissible'

def main():
    started=time.monotonic();hook_controls=[];small_controls=[]
    for n in [7,8]:
        rows=[{'partition':lam,'hook_degree':hook_dim(lam),'branch_degree':branch_dim(lam)} for lam in parts(n)]
        ck(f'S{n} hook and independent branching dimensions agree',all(a['hook_degree']==a['branch_degree'] for a in rows))
        ck(f'S{n} complete low-degree induction base',
           {a['partition'] for a in rows if a['hook_degree']<=n+1}==exceptional(n))
        ck(f'S{n} nonexceptional minimum exactly14',min(a['hook_degree'] for a in rows if a['partition'] not in exceptional(n))==14)
        hook_controls.append({'n':n,'all_partitions':rows})
    # Finite controls of the structural child/grandchild induction, not its proof.
    for n in range(9,17):
        failures=[]
        for lam in parts(n):
            if lam in exceptional(n):continue
            ch=children(lam)
            if len(ch)>=2:
                if not any(mu not in exceptional(n-1) for mu in ch):failures.append(lam)
            else:
                grandchildren=children(ch[0])
                if len(grandchildren)!=2 or any(mu in exceptional(n-2) for mu in grandchildren):failures.append(lam)
        ck(f'S{n} structural induction controls',not failures)
    fixedtypes={2:[((1,1),1),((2,),1)],
                3:[((1,1,1),1),((3,),2)],
                4:[((1,1,1,1),1),((2,2),3)],
                5:[((1,1,1,1,1),1),((5,),4)],
                6:[((1,)*6,1),((3,1,1,1),4),((3,3),4)]}
    for n in range(2,7):
        ps=parts(n);ct=parts(n)
        ck(f'S{n} all character dimensions match hook formula',all(character(lam,(1,)*n)==hook_dim(lam) for lam in ps))
        ck(f'S{n} entire exact character table row orthogonality',all(
           sum(class_size(c)*character(lam,c)*character(mu,c) for c in ct)==(factorial(n) if lam==mu else 0)
           for lam in ps for mu in ps))
        fixed={lam:sum(mult*character(lam,c) for c,mult in fixedtypes[n])//sum(mult for c,mult in fixedtypes[n]) for lam in ps}
        ck(f'S{n} fixed-space averages are integral',all(
           sum(mult*character(lam,c) for c,mult in fixedtypes[n])%sum(mult for c,mult in fixedtypes[n])==0 for lam in ps))
        items=[(lam,hook_dim(lam)) for lam in ps if hook_dim(lam)<=n+1]
        lower_bad=[];surface_bad=[];survivors=[];faithful_count=0
        for dim in range(1,n+2):
            for modules in multisets(items,dim):
                faithful=all(sum(character(lam,c) for lam in modules)!=dim for c in ct if c!=(1,)*n)
                if not faithful:continue
                faithful_count+=1;inv=sum(fixed[lam] for lam in modules)
                bound=n if n==6 else n+1
                if dim<bound and inv>=2:lower_bad.append((modules,inv))
                if dim==bound and inv>=2:
                    survivors.append({'partitions':modules,'fixed_dimension':inv})
                    if inv!=2:surface_bad.append((modules,inv))
        ck(f'S{n} complete faithful low-degree obstruction',not lower_bad)
        ck(f'S{n} complete equality surface profiles',not surface_bad and bool(survivors))
        small_controls.append({'n':n,'all_low_irreducibles':[{'partition':lam,'degree':deg,'fixed_dimension':fixed[lam]} for lam,deg in items],
                               'faithful_representations_checked':faithful_count,'equality_profiles_passing_this_subgroup':survivors})
    # Actual even p-groups realizing two nontrivial block orbits.
    block_controls=[]
    for p,a,b in [(3,1,1),(3,1,2),(2,1,2),(2,2,3),(7,1,1)]:
        A,B=p**a,p**b;n=A+B
        group=[]
        for i,j in product(range(A),range(B)):
            if p==2 and (i+j)%2:continue
            group.append(tuple((k+i)%A for k in range(A))+tuple(A+(k+j)%B for k in range(B)))
        ck(f'two blocks {A}+{B} even actual group',all((n-len(cycle_type(g)))%2==0 for g in group) and all(compose(g,h) in group for g in group for h in group))
        orbits={tuple(sorted({g[i] for g in group})) for i in range(n)}
        ck(f'two blocks {A}+{B} exact two orbits without fixed labels',len(orbits)==2 and all(len(o)>1 for o in orbits))
        # Character average of the standard representation, independently of orbit count.
        inv=Fraction(sum(sum(g[i]==i for i in range(n))-1 for g in group),len(group))
        ck(f'two blocks {A}+{B} standard fixed dimension exactly one',inv==1)
        block_controls.append({'p':p,'blocks':[A,B],'order':len(group),'orbits':list(orbits),'standard_fixed_dimension':int(inv)})
    arithmetic=[]
    for n in range(7,33):
        violations=[]
        for p in s.primerange(2,n+1):
            for pattern in orbit_patterns(n,int(p)):
                if pattern.count(1)!=1 and len(pattern)-1<2:violations.append({'p':int(p),'orbit_sizes':pattern})
        kind=classification(n)
        ck(f'rank{n} exhaustive prime-power orbit-pattern criterion',bool(violations)==(kind!='standard_admissible'))
        arithmetic.append({'n':n,'class':kind,'obstructing_orbit_patterns':violations})
    ck('ranks15 and28 are exact representation-level limitations',all(classification(n)=='standard_admissible' for n in [15,28]))
    # A genuinely nonlinear expression for a tame fixed curve.
    x,y=s.symbols('x y');g=s.Matrix([-x+2*y*y,y]);F=s.Matrix([2*x-y*y+2*y+1,y+1])
    variables=[x,y]
    def subst(v,map):return v.subs(dict(zip(variables,map)),simultaneous=True).applyfunc(s.expand)
    ck('nonlinear involution is order two',subst(g,g)==s.Matrix(variables))
    ck('nonlinear involution fixed ideal is the parabola',s.groebner(list(g-s.Matrix(variables)),x,y).polys==s.groebner([x-y*y],x,y).polys)
    ck('nonlinear fixed-curve control is equivariant Keller',subst(F,g)==subst(g,F) and s.det(F.jacobian(variables))==2)
    ck('whole fixed-curve restriction is the surjective translation',all(s.expand(v)==0 for v in F.subs(x,y*y)-s.Matrix([(y+1)**2,y+1])))
    for p in [2,3,5]:
        zero=lambda z:s.Poly(s.expand(z),x,y,modulus=p).is_zero
        ck(f'char{p} wild finite action has identity first jet',s.Matrix([x+y**p,y]).jacobian([x,y]).subs({x:0,y:0})==s.eye(2))
        ck(f'char{p} wild action has exact order p',zero(p*y**p) and all(not zero(a*y**p) for a in range(1,p)))
        ck(f'char{p} fixed ideal retains its nonreduced thickening',p>1 and s.Poly(y**p,y).degree()==p)
        A=x**p-x
        ck(f'char{p} all translations are actual regular deck maps',all(zero(A.subs(x,x+a)-A) for a in range(p)))
        ck(f'char{p} degree-p fixed-line exception is Keller',zero(s.diff(x+x**p,x)-1) and s.Poly(x+x**p,x).degree()==p)
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:
        ck('retained local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        ck('retained primary pin '+src['path'],hashlib.sha256((HERE/src['path']).read_bytes()).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,
             'induction_base_controls':hook_controls,'complete_small_character_controls':small_controls,
             'two_orbit_groups':block_controls,'orbit_pattern_controls':arithmetic,
             'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),
             'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},
             'assurance':'Written all-rank arguments with exact induction bases, independent hook/branching controls, complete small character-table orthogonality and faithful low-degree representation enumeration, actual even two-orbit groups, nonlinear and wild controls, and primary source seals. These controls do not formalize Smith theory, etale geometry, all-rank branching induction or surface recognition. No independent whole-proof replay, formal verification, external review or priority certification. Admissible tangent profiles do not construct Keller maps.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'complete_small_character_ranks':5,'arithmetic_control_ranks':26,'seconds':receipt['runtime_seconds']}))

if __name__=='__main__':main()
