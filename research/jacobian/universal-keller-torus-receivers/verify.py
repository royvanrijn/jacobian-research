#!/usr/bin/env python3
"""Exact norm-torus, group, lattice and factor controls; not formalization."""
from pathlib import Path
from itertools import product,permutations,combinations
from math import prod,gcd,lcm,factorial
import hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def compose(a,b):return tuple(a[b[i]] for i in range(len(a)))
def inv(a):return tuple(a.index(i) for i in range(len(a)))
def group(gens,n):
    one=tuple(range(n));seen={one};todo=[one]
    for g in todo:
        for h in gens:
            x=compose(g,h)
            if x not in seen:seen.add(x);todo.append(x)
    return frozenset(seen)
def abelian(G):return all(compose(g,h)==compose(h,g) for g in G for h in G)
def transitive(G,n):return {g[0] for g in G}==set(range(n))
def primitive(G,n):
    for size in range(2,n):
        if n%size:continue
        for rest in combinations(range(1,n),size-1):
            B={0,*rest}
            if all(not B&{g[b] for b in B} or B=={g[b] for b in B} for g in G):return False
    return True
def cycle_type(g):
    unseen=set(range(len(g)));sizes=[]
    while unseen:
        a=next(iter(unseen));b=a;size=0
        while b in unseen:unseen.remove(b);size+=1;b=g[b]
        sizes.append(size)
    return tuple(sorted(sizes))
def dim(gb,n):
    leading=[p.LM(order=gb.order).exponents for p in gb.polys]
    if (0,)*n in leading:return 0
    bounds=[]
    for i in range(n):
        values=[m[i] for m in leading if m[i] and all(v==0 for j,v in enumerate(m) if j!=i)]
        assert values,'positive dimensional control quotient';bounds.append(min(values))
    return sum(not any(all(a>=b for a,b in zip(v,m)) for m in leading) for v in product(*(range(b) for b in bounds)))
def main():
    started=time.monotonic();a,b=s.symbols('a b');N=a*a-a*b+b*b
    F=s.Matrix([a**3-3*a*b*b+b**3,3*a*b*(a-b)])
    tau=s.Matrix([[0,-1],[1,-1]]);iota=s.Matrix([[1,-1],[0,-1]]);v=s.Matrix([a,b])
    ck('cubic source relations over integers',tau**3==s.eye(2) and iota**2==s.eye(2) and iota*tau*iota==tau**2)
    ck('cubic source group has order six',len({tuple((tau**i*iota**j)) for i in range(3) for j in range(2)})==6)
    def pull(poly,M):return s.expand(poly.subs({a:(M*v)[0],b:(M*v)[1]},simultaneous=True))
    ck('source generators preserve norm over integers',pull(N,tau)==N and pull(N,iota)==N)
    ck('cube map preserves norm hypersurface',s.expand(N.subs({a:F[0],b:F[1]},simultaneous=True)-N**3)==0)
    ck('cube equivariance translation killed on target',all(s.expand(pull(F[j],tau)-F[j])==0 for j in range(2)))
    ck('cube equivariance inversion on target',all(s.expand(pull(F[j],iota)-(iota*F)[j])==0 for j in range(2)))
    W=s.Matrix([a-2*b,2*a-b]);WF=s.Matrix([F[0]-2*F[1],2*F[0]-F[1]])
    ck('global invariant differential is multiplication by three',all(s.expand(q)==0 for q in F.jacobian([a,b])*W-3*WF))
    ck('tangent generator has norm three on the norm-one torus',s.expand(N.subs({a:W[0],b:W[1]},simultaneous=True)-3*N)==0)
    points=[(1,0),(0,1),(-1,-1)];fibre_controls=[]
    matrix=s.Matrix([[s.diff(N,a),s.diff(N,b)]]).col_join(F.jacobian([a,b]));minors=[s.expand(matrix[list(rows),:].det()) for rows in combinations(range(3),2)]
    for p in [0,2,5,7,11]:
        opts={'domain':s.QQ} if not p else {'modulus':p}
        gb=s.groebner([N-1,F[0]-1,F[1]],a,b,**opts)
        ck(f'char{p} complete cubic marked quotient rank three',dim(gb,2)==3)
        ck(f'char{p} whole cubic curve etale ideal',list(s.groebner([N-1]+minors,a,b,**opts))==[1])
        ck(f'char{p} source three-translation has no fixed point',list(s.groebner([N-1]+list(tau*v-v),a,b,**opts))==[1])
        permutations_found=set()
        pts=[tuple(c%p for c in P) if p else P for P in points]
        for i,j in product(range(3),range(2)):
            M=tau**i*iota**j;permutation=[]
            for P in pts:
                value=tuple(M*s.Matrix(P));value=tuple(int(c)%p for c in value) if p else value
                permutation.append(pts.index(value))
            permutations_found.add(tuple(permutation))
        ck(f'char{p} actual full natural S3 marked action',permutations_found==set(permutations(range(3))))
        fibre_controls.append({'rank':3,'characteristic':p,'whole_fibre_dimension':3})
    ck('characteristic three excludes the cubic differential',all(s.Poly(s.diff(q,x),a,b,modulus=3).is_zero for q in F for x in [a,b]))
    # Quartic finite torus map and the entire natural S4 action.
    x,y,z=s.symbols('x y z');Q=x*y*z-1;quartic=[x*x,y*y,z*z];roots=[v for v in product([-1,1],repeat=3) if prod(v)==1]
    full=set();translation_group=set();linear_group={}
    for pi in permutations(range(3)):
        pg=[]
        for root in roots:pg.append(roots.index(tuple(root[pi[i]] for i in range(3))))
        linear_group[pi]=tuple(pg)
        for signs in roots:
            perm=tuple(roots.index(tuple(signs[i]*root[pi[i]] for i in range(3))) for root in roots)
            full.add(perm)
            if pi==(0,1,2):translation_group.add(perm)
            pulled=[(signs[i]*[x,y,z][pi[i]])**2 for i in range(3)]
            assert all(s.expand(pulled[i]-quartic[pi[i]])==0 for i in range(3))
    ck('quartic natural action is all twenty-four permutations',full==set(permutations(range(4))))
    ck('quartic geometric group is regular V4',len(translation_group)==4 and abelian(translation_group) and transitive(translation_group,4))
    ck('quartic V4 is imprimitive',not primitive(translation_group,4))
    c3=linear_group[(1,2,0)];swap=linear_group[(1,0,2)]
    A4=group(list(translation_group)+[c3],4);S4=group(list(translation_group)+[c3,swap],4);D8=group(list(translation_group)+[swap],4)
    ck('connected cyclic cubic resolvent gives primitive arithmetic A4',len(A4)==12 and primitive(A4,4))
    ck('full cubic resolvent gives primitive arithmetic S4',len(S4)==24 and primitive(S4,4))
    ck('quadratic resolvent arithmetic D8 remains imprimitive',len(D8)==8 and not primitive(D8,4))
    pairings=[frozenset([frozenset([0,j]),frozenset(set(range(4))-{0,j})]) for j in [1,2,3]]
    def induced(g):return tuple(pairings.index(frozenset(frozenset(g[v] for v in block) for block in P)) for P in pairings)
    ck('all V4 translations fix all three quartic resolvent labels',all(induced(g)==(0,1,2) for g in translation_group))
    ck('quartic resolvent quotient is full S3 with kernel V4',len({induced(g) for g in full})==6 and {g for g in full if induced(g)==(0,1,2)}==translation_group)
    for p in [0,3,5,7]:
        opts={'domain':s.QQ} if not p else {'modulus':p};gb=s.groebner([Q]+[q-1 for q in quartic],x,y,z,**opts)
        ck(f'char{p} whole quartic marked quotient rank four',dim(gb,3)==4)
        matrix=s.Matrix([[s.diff(Q,v) for v in [x,y,z]]]).col_join(s.Matrix(quartic).jacobian([x,y,z]));ms=[s.expand(matrix[list(rows),:].det()) for rows in combinations(range(4),3)]
        ck(f'char{p} whole quartic torus etale ideal',list(s.groebner([Q]+ms,x,y,z,**opts))==[1]);fibre_controls.append({'rank':4,'characteristic':p,'whole_fibre_dimension':4})
    # Exact index and equivariant kernel controls for the all-group construction.
    lattice_controls=[]
    for moduli in [(2,),(3,),(4,),(2,2),(2,2,2),(6,)]:
        B=list(product(*(range(n) for n in moduli)));zero=(0,)*len(moduli);B.remove(zero);d=len(B);indices=[B.index(tuple(int(i==j) for i in range(len(moduli)))) for j in range(len(moduli))];K=s.zeros(d);column=0
        for i,index in enumerate(indices):K[index,column]=moduli[i];column+=1
        for j,bv in enumerate(B):
            if j in indices:continue
            K[j,column]=1
            for i,index in enumerate(indices):K[index,column]=-bv[i]
            column+=1
        P=s.Matrix(len(moduli),d,lambda i,j:B[j][i]);size=prod(moduli)
        ck(f'abelian{moduli} exact kernel lattice index',abs(K.det())==size)
        ck(f'abelian{moduli} character quotient relations',all(int((P*K)[i,j])%moduli[i]==0 for i in range(len(moduli)) for j in range(d)))
        exponent=lcm(*moduli)
        for unit in range(1,exponent):
            if gcd(unit,exponent)!=1:continue
            perm=[B.index(tuple(unit*bv[i]%moduli[i] for i in range(len(moduli)))) for bv in B];M=s.zeros(d)
            for j,index in enumerate(perm):M[index,j]=1
            action=K.inv()*M*K
            ck(f'abelian{moduli} cyclotomic{unit} kernel invariant',all(v.q==1 for v in action) and abs(action.det())==1)
            inverse=pow(unit,-1,exponent)
            ck(f'abelian{moduli} cyclotomic{unit} constant kernel descent',all((unit*sum((exponent//n)*(inverse*bv[i]%n)*av[i] for i,n in enumerate(moduli))-sum((exponent//n)*bv[i]*av[i] for i,n in enumerate(moduli)))%exponent==0 for bv in B for av in product(*(range(n) for n in moduli))))
        lattice_controls.append({'abelian_factors':moduli,'dimension':d,'kernel_order':size,'index':int(abs(K.det()))})
    # The normalizer condition does not say the regular abelian group is inside G.
    cyc6=group([tuple((i+1)%6 for i in range(6))],6);regS3=group([tuple((i+2)%6 for i in range(6)),tuple((1-i)%6 for i in range(6))],6)
    ck('regular degree-six S3 normalizes C6 without containing it',len(regS3)==6 and transitive(regS3,6) and not abelian(regS3) and not cyc6<=regS3 and all(compose(compose(g,h),inv(g)) in cyc6 for g in regS3 for h in cyc6))
    # Small normal closures check the exact low-rank symmetric exceptions.
    normal_controls=[]
    for n in range(2,6):
        Sn=list(permutations(range(n)));representatives={cycle_type(g):g for g in Sn};regular=[]
        for typ,g in representatives.items():
            closure=group(list({compose(compose(h,g),inv(h)) for h in Sn}),n)
            if len(closure)>1 and abelian(closure):
                ck(f'S{n} abelian nontrivial normal closure is regular',len(closure)==n and transitive(closure,n));regular.append(len(closure))
        ck(f'S{n} exact abelian regular normal census',bool(regular)==(n<=4));normal_controls.append({'rank':n,'abelian_regular_orders':regular})
    for n in range(2,25):
        cyc=group([tuple((i+1)%n for i in range(n))],n);swap=list(range(n));swap[0],swap[1]=1,0;swap=tuple(swap)
        normal=all(compose(compose(swap,g),swap) in cyc for g in cyc)
        ck(f'cyclic degree{n} symmetric-normality boundary',normal==(n<=3))
    # Universal stable norm obstruction and a different intermediate torus.
    aa,bb=s.symbols('aa bb');swap=s.Matrix([[0,1],[1,0]]);endom=s.Matrix([[aa,bb],[bb,aa]])
    ck('quadratic restriction torus endomorphism norm',endom*swap==swap*endom and endom.det()==aa*aa-bb*bb)
    ck('no degree two endomorphism over quadratic restriction torus',all((a*a-b*b)%4!=2 for a,b in product(range(4),repeat=2)))
    L=s.Matrix([[1,1],[1,-1]]);twist=L.inv()*swap*L
    ck('two quadratic factors exist through a different torus lattice',L*L==2*s.eye(2) and abs(L.det())==2 and twist==s.diag(1,-1))
    ck('intermediate torus cannot be made the source by affine stabilization',twist.applyfunc(lambda x:x%2)==s.eye(2) and swap!=s.eye(2))
    # Nonnormal intermediate defeats prime finite atomicity.
    ss,t,q=s.symbols('ss t q');node=q*q-t*(t-1)**2
    ck('nonnormal prime-factor countercontrol whole equation',s.expand(node.subs({t:ss*ss,q:ss*(ss*ss-1)}))==0)
    ck('nonnormal factor is birational',s.cancel((q/(t-1)).subs({t:ss*ss,q:ss*(ss*ss-1)})-ss)==0)
    ck('nonnormal intermediate really has a node',node.subs({t:1,q:0})==0 and all(s.diff(node,v).subs({t:1,q:0})==0 for v in [t,q]))
    ck('nonnormal factor merges distinct source points',[(r*r,r*(r*r-1)) for r in [1,-1]]==[(1,0),(1,0)])
    # Frozen cubic sign and cyclic-quartic failure are exact permutation statements.
    for scalar in [1,2]:
        for translation in range(3):
            p=tuple((scalar*i+translation)%3 for i in range(3));sign=(-1)**sum(p[i]>p[j] for i in range(3) for j in range(i+1,3))
            ck(f'cubic sign frozen scalar{scalar} shift{translation}',sign==(1 if scalar==1 else -1))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('retained source pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:ck('primary source pin '+src['path'],hashlib.sha256((HERE/src['path']).read_bytes()).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'complete_marked_fibre_controls':fibre_controls,'permutation_kernel_lattice_controls':lattice_controls,'symmetric_normal_closure_controls':normal_controls,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Written all-degree finite factor closure, characteristic-zero torus-cylinder structure, all-characteristic tame torus normalizer classification and exact symmetric constructions. Exact controls check integer norm/cube identities, whole marked fibres and etaleness ideals, actual label actions, invariant kernel lattices and singular-intermediate boundaries. They do not formalize Galois descent, Hopf-Galois theory, all-rank torus duality, simple connectedness or fieldwise cohomology. Primary sources credited. No independent whole-proof replay, formal verification, external review or novelty certification.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'complete_fibres':len(fibre_controls),'lattice_cases':len(lattice_controls),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
