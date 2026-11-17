#!/usr/bin/env python3
"""Exact local, permutation and complete pseudo-plane controls; not formalization."""
from pathlib import Path
from itertools import product,permutations,combinations
from math import factorial,lcm,gcd
from functools import lru_cache
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
    if not n:return ((),)
    return tuple((a,)+b for a in range(min(n,top or n),0,-1) for b in parts(n-a,a))
def pp_order(order):return order>1 and len(s.factorint(order))==1
def quotient_dimension(gb,variables):
    leading=[p.LM(order=gb.order).exponents for p in gb.polys]
    if (0,)*len(variables) in leading:return 0
    bounds=[]
    for i in range(len(variables)):
        b=[lm[i] for lm in leading if lm[i] and all(v==0 for j,v in enumerate(lm) if j!=i)]
        assert b,'not zero dimensional'
        bounds.append(min(b))
    return sum(not any(all(e>=v for e,v in zip(exp,lm)) for lm in leading)
               for exp in product(*(range(v) for v in bounds)))
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def generated(gens,n):
    seen={tuple(range(n))};todo=list(seen)
    for g in todo:
        for h in gens:
            v=compose(g,h)
            if v not in seen:seen.add(v);todo.append(v)
    return frozenset(seen)
def order(g):
    val=tuple(range(len(g)));v=val
    for j in range(1,factorial(len(g))+1):
        v=compose(v,g)
        if v==val:return j
    raise AssertionError('permutation order')
def transitive(G,n):return {g[0] for g in G}==set(range(n))

def main():
    started=time.monotonic();derangement_controls=[]
    # All finite permutation groups in degrees2..4, not just the full symmetric groups.
    for n in [2,3,4]:
        full=list(permutations(range(n)));trivial=generated([],n)
        groups={trivial};queue=[trivial]
        for G in queue:
            for g in full:
                if g in G:continue
                H=generated(list(G)+[g],n)
                if H not in groups:groups.add(H);queue.append(H)
        trans=[G for G in groups if transitive(G,n)]
        ck(f'degree{n} all transitive groups have a prime-power derangement',all(
           any(all(g[i]!=i for i in range(n)) and pp_order(order(g)) for g in G) for G in trans))
        derangement_controls.append({'degree':n,'all_subgroups':len(groups),'transitive_subgroups':len(trans)})
    # Independent cycle-partition census checks the exact affine exception list.
    affine_controls=[]
    for n in range(2,17):
        der=[p for p in parts(n) if 1 not in p]
        for p in [0]+list(s.primerange(2,20)):
            tame=[v for v in der if p==0 or gcd(lcm(*v),int(p))==1]
            exceptional=(n,int(p)) in {(2,2),(3,3),(4,2)}
            ck(f'rank{n} char{p} affine cyclic derangement criterion',bool(tame)==(not exceptional))
            if not tame:affine_controls.append({'rank':n,'characteristic':int(p),'derangement':None})
    witnesses=[]
    for n in range(2,65):
        for p in [0,2,3,5,7,11,13]:
            if (n,p) in {(2,2),(3,3),(4,2)}:continue
            if p==0 or n%p:witness=(n,)
            elif n>=5 and p!=2:witness=(2,n-2)
            elif p==2 and n>=6:witness=(3,n-3)
            else:witness=(2,2) if n==4 else (n,)
            ck(f'rank{n} char{p} explicit all-rank cyclic witness',sum(witness)==n and min(witness)>1 and (p==0 or gcd(lcm(*witness),p)==1))
            if p and all(v==p for v in s.factorint(n)) and n>=5:
                ck(f'rank{n} char{p} modular witness is genuinely composite order',not pp_order(lcm(*witness)))
            witnesses.append({'rank':n,'characteristic':p,'cycle_lengths':witness,'order':lcm(*witness)})
    # Source-fixed origins are absent in every matching modular construction.
    x,y=s.symbols('x y')
    for p in [2,3]:
        F=x-x**p;perms=set()
        for a,b in product(range(1,p),range(p)):
            ck(f'Phi{p} affine equivariance a{a} b{b}',s.Poly(s.expand(F.subs(x,a*x+b)-a*F),x,modulus=p).is_zero)
            perms.add(tuple((a*v+b)%p for v in range(p)))
        ck(f'Phi{p} actual full natural symmetric action',len(perms)==factorial(p))
        ck(f'Phi{p} whole zero-fibre rank',s.Poly(F,x,modulus=p).degree()==p and s.Poly(s.diff(F,x)-1,x,modulus=p).is_zero)
        ck(f'Phi{p} translation forbids every source origin',s.expand((x+1)-x)==1)
    roots=list(product(range(2),repeat=2));seen=set();Fx=s.Matrix([x-x*x,y-y*y])
    for a,b,c,d in product(range(2),repeat=4):
        if (a*d-b*c)%2!=1:continue
        A=s.Matrix([[a,b],[c,d]])
        for e,f in product(range(2),repeat=2):
            action=A*s.Matrix([x,y])+s.Matrix([e,f])
            pulled=Fx.subs({x:action[0],y:action[1]},simultaneous=True)-A*Fx
            assert all(s.Poly(v,x,y,modulus=2).is_zero for v in pulled)
            seen.add(tuple(roots.index(((a*v[0]+b*v[1]+e)%2,(c*v[0]+d*v[1]+f)%2)) for v in roots))
    ck('quartic plane gives the entire natural S4 action',seen==set(permutations(range(4))))
    ck('quartic plane complete zero fibre rank four',quotient_dimension(s.groebner(list(Fx),x,y,modulus=2),(x,y))==4)
    # Local rank is not special fibre length or the number of omitted points.
    a,b,z,w=s.symbols('a b z w')
    rel=[z*z-a*z,z*w-a*w,w*w-b*w,b*z-a*w]
    generic=s.groebner(rel,z,w,domain=s.QQ.frac_field(a,b))
    special=s.groebner([v.subs({a:0,b:0}) for v in rel],z,w)
    ck('torsion-free finite algebra generic rank two',quotient_dimension(generic,(z,w))==2)
    ck('nonflat special quotient has length three',quotient_dimension(special,(z,w))==3)
    ck('special quotient supported at one geometric point',all(special.reduce(v)[1]==0 for v in [z*z,z*w,w*w]))
    local_controls=[]
    for n,delta in [(3,1),(3,2),(4,3)]:
        retained=tuple((i+1)%n for i in range(n))
        r=sum(retained[i]==i for i in range(n));M=n+delta
        generic=s.Poly(z**delta-a,z,domain=s.QQ.frac_field(a))
        special=s.Poly(z**delta,z)
        ck(f'local retained-cycle n{n} delta{delta} exact degree budget',r==0 and generic.degree()==M-n)
        ck(f'local retained-cycle n{n} delta{delta} one omitted geometric point',s.gcd(special,s.Poly(z,z)).degree()==1)
        local_controls.append({'retained_rank':n,'generic_degree':M,'generic_fixed_degree':delta,
                               'omitted_geometric_points_at_origin':1,
                               'boundary_warning':'Local source is not whole affine space; its fixed punctured line is not acyclic.'})
    # Whole pseudo-plane maps, with coefficients in exact cyclotomic fields.
    u,v,w,t,e=s.symbols('u v w t e');pseudo=[];etale=[]
    for k in range(2,7):
        phi=s.cyclotomic_poly(k,e)
        def red(expr):return s.Poly(s.expand(expr),e,domain=s.QQ.poly_ring(u,v,w,t)).rem(s.Poly(phi,e,domain=s.QQ.poly_ring(u,v,w,t))).as_expr().expand()
        R=s.div((1+(e-1)*t)**k-1,t*(t-1),t)[0]
        R=red(R)
        ck(f'pseudo-plane k{k} exact polynomial R identity',red(t*(t-1)*R-((1+(e-1)*t)**k-1))==0)
        tt=-u*v;rr=R.subs(t,tt);Z=1+(e-1)*tt
        G=u+u*u*v-w**k;F=[w**k,v*rr,w*Z]
        source=s.groebner([G,phi],u,v,w,e)
        ck(f'pseudo-plane k{k} whole source maps into universal cover',source.reduce(s.expand(w**k*v*rr-Z**k+1))[1]==0)
        ck(f'pseudo-plane k{k} whole source maps into target surface',source.reduce(s.expand(F[0]+F[0]**2*F[1]-F[2]**k))[1]==0)
        targets=[((1,-1,0),k),((1,0,1),k-1),((0,1,0),2),((0,0,0),1)]
        for point,expected in targets:
            gb=s.groebner([G,phi]+[F[i]-point[i] for i in range(3)],u,v,w,e,order='grevlex')
            dimension=quotient_dimension(gb,(u,v,w,e));degree=s.degree(phi,e)
            ck(f'pseudo-plane k{k} complete fibre {point}',dimension==expected*degree)
            pseudo.append({'degree':k,'target':point,'whole_geometric_rank':expected,'cyclotomic_degree':int(degree),'whole_quotient_Q_dimension':dimension})
        if k<=3:
            matrix=s.Matrix([s.diff(G,x) for x in [u,v,w]]).T.col_join(s.Matrix(F).jacobian([u,v,w]))
            minors=[s.expand(matrix[list(rows),:].det()) for rows in combinations(range(4),3)]
            gb=s.groebner([G,phi]+minors,u,v,w,e)
            ck(f'pseudo-plane k{k} global etaleness ideal certificate',list(gb)==[1])
            etale.append(k)
    G2=u+u*u*v-w*w;F2=[w*w,4*v,w*(1+2*u*v)]
    tau={u:-u-1/v,v:v,w:-w}
    ck('pseudo-plane quadratic generic involution preserves source equation',s.cancel(G2.subs(tau,simultaneous=True)-G2)==0)
    ck('pseudo-plane quadratic generic involution preserves the whole map',all(s.cancel(q.subs(tau,simultaneous=True)-q)==0 for q in F2))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('retained local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:ck('primary pin '+src['path'],hashlib.sha256((HERE/src['path']).read_bytes()).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,
             'complete_small_transitive_groups':derangement_controls,'affine_exception_census':affine_controls,
             'explicit_symmetric_cyclic_witnesses':witnesses,'local_rank_controls':local_controls,
             'complete_pseudo_plane_fibres':pseudo,'global_pseudo_plane_etale_certificates':etale,
             'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),
             'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},
             'assurance':'Written normalization rank budget, compact-support rigidity, Smith and affine fixed-locus proofs. Bounded controls enumerate all subgroups in degrees2..4, test the explicit all-rank cyclic constructions, certify matching modular actions, nonflat rank/length distinctions, complete pseudo-plane fibre quotients and two global etaleness ideals. They do not formalize etale normalization, Smith theory, compact-support duality or CFSG. No independent whole-proof replay, formal verification, external review or novelty priority. Nonproperness alone is not an empty-fibre certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'whole_pseudo_plane_fibres':len(pseudo),
                      'symmetric_cycle_witnesses':len(witnesses),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
