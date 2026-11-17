#!/usr/bin/env python3
"""Small exact representation, block and volume controls for written all-degree proofs."""
from pathlib import Path
from itertools import product, permutations
import hashlib,importlib.util,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('retained_permutation_helpers',HERE.parent/'universal-keller-homogeneous-receivers/verify.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
generated,mul,inv=old.generated,old.mul,old.inv
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def partitions(values):
    if not values:yield ();return
    first,*tail=values
    for P in partitions(tail):
        yield (frozenset([first]),)+P
        for j in range(len(P)):yield P[:j]+(P[j]|{first},)+P[j+1:]
def main():
    started=time.monotonic();R=s.Matrix([[-1,-1],[1,0]]);S=s.Matrix([[0,1],[1,0]])
    # Holonomy commutants: one copy is scalar; repeated copies genuinely enlarge it.
    for m in [1,2,3]:
        n=2+2*m;RR=s.diag(1,1,*([R]*m));SS=s.diag(-1,1,*([S]*m))
        V=s.symbols('v:'+str(n*n));X=s.Matrix(n,n,V)
        equations=list(X*RR-RR*X)+list(X*SS-SS*X)
        A,_=s.linear_eq_to_matrix(equations,V)
        ck(f'm{m} exact rational commutant dimension',n*n-A.rank()==2+m*m)
    # Actual integral S4 augmentation action on the multiplicity-three lattice.
    def aug(g):
        columns=[]
        for j in range(3):
            v=[0]*4;v[g[j]]+=1;v[g[3]]-=1;columns.append(s.Matrix(v[:3]))
        return s.Matrix.hstack(*columns)
    G4=list(permutations(range(4)));mats={g:aug(g) for g in G4}
    ck('faithful integral S4 representation on three copies',len({tuple(M) for M in mats.values()})==24 and all(abs(M.det())==1 for M in mats.values()))
    ck('S4 augmentation representation multiplication',all(mats[mul(g,h)]==mats[g]*mats[h] for g in G4 for h in G4))
    RR=s.diag(R,R,R);SS=s.diag(S,S,S)
    ck('S4 monomial mixing commutes with actual S3 holonomy',all(s.kronecker_product(M,s.eye(2))*RR==RR*s.kronecker_product(M,s.eye(2)) and s.kronecker_product(M,s.eye(2))*SS==SS*s.kronecker_product(M,s.eye(2)) for M in mats.values()))
    # General block projection: S3 wr C2 has two blocks of three, with nonfaithful block stabilizer.
    generators=[(1,0,2,3,4,5),(1,2,0,3,4,5),(3,4,5,0,1,2)]
    G=generated(generators,6);valid=[]
    for P in partitions(list(range(6))):
        blocks=set(P)
        if all({frozenset(g[i] for i in B) for B in P}==blocks for g in generators):valid.append(P)
    ck('degree6 wreath group order',len(G)==72)
    ck('all invariant partitions found exactly',sorted((len(P),len(P[0])) for P in valid)==[(1,6),(2,3),(6,1)])
    B=frozenset(range(3));J=[g for g in G if frozenset(g[i] for i in B)==B]
    image={tuple(g[i] for i in range(3)) for g in J}
    ck('block stabilizer remains transitive with kernel',len(J)==36 and len(image)==6 and len({g[0] for g in J})==3)
    # Primitive prime-square boundary: exact invariant lines, not a growing block search.
    controls=[]
    for p in [2,3,5,7,11]:
        def mv(M,v):return tuple(sum(int(M[i,j])*v[j] for j in range(2))%p for i in range(2))
        lines=[(1,a) for a in range(p)]+[(0,1)]
        def same_line(v,w):return (v[0]*w[1]-v[1]*w[0])%p==0
        stable=[v for v in lines if same_line(v,mv(R,v)) and same_line(v,mv(S,v))]
        ck(f'F{p} standard invariant-line boundary',len(stable)==(1 if p==3 else 0))
        if p in [2,3,5]:
            pts=list(product(range(p),repeat=2))
            T1=tuple(pts.index(((x+1)%p,y)) for x,y in pts);T2=tuple(pts.index((x,(y+1)%p)) for x,y in pts)
            r=tuple(pts.index(mv(R,v)) for v in pts);ss=tuple(pts.index(mv(S,v)) for v in pts)
            H=generated([T1,T2,r,ss],p*p)
            ck(f'F{p} full affine monodromy order',len(H)==6*p*p)
            if p==3:
                diagonal=frozenset(i for i,v in enumerate(pts) if v[0]==v[1])
                blockset={frozenset(g[i] for i in diagonal) for g in H}
                ck('degree9 actual three-block intermediate',len(blockset)==3 and all(len(B)==3 for B in blockset) and set().union(*blockset)==set(range(9)))
            if p==2:
                ck('quartic monodromy exactly S4',H==frozenset(permutations(range(4))))
                central=[g for g in permutations(range(4)) if g[0]==0 and mul(g,r)==mul(r,g) and mul(g,ss)==mul(ss,g)]
                ck('faithful S3 centralizer on V4 automorphisms trivial',central==[tuple(range(4))])
        controls.append({'prime':p,'stable_lines':[list(v) for v in stable],'primitive_standard_module':not stable})
    # Real global differential identities for the entire displayed dimension grid.
    x,z,k=s.symbols('x z k',nonzero=True)
    for m in [1,2,3]:
        coords=[x,z]+list(s.symbols('a:'+str(2*m),nonzero=True));aa=2 if m%2 else 1;bb=1 if m%2 else 2
        density=z**aa/s.prod(coords)
        rotation=s.Matrix([k*x,z]+[v for j in range(m) for v in [1/(coords[2+2*j]*coords[3+2*j]),coords[2+2*j]]])
        reflection=s.Matrix([1/x,-z]+[v for j in range(m) for v in [coords[3+2*j],coords[2+2*j]]])
        for name,F in [('rotation',rotation),('reflection',reflection)]:
            pulled=density.subs(dict(zip(coords,F)),simultaneous=True)
            ck(f'm{m} invariant global volume '+name,s.factor(F.jacobian(coords).det()*pulled/density)==1)
        for n in [2,3,4,5]:
            F=s.Matrix([x,z/n**bb,coords[2]**n,coords[3]**n]+coords[4:])
            pulled=density.subs(dict(zip(coords,F)),simultaneous=True)
            ck(f'm{m} Phi{n} exact volume compensation',s.factor(F.jacobian(coords).det()*pulled/density)==1)
    # Test an off-diagonal integral multiplicity map and its degree-square determinant.
    u,v,w,t=s.symbols('u v w t',nonzero=True);coords=[u,v,w,t]
    F=s.Matrix([u**2*w,v**2*t,u*w**2,v*t**2])
    ck('off-diagonal multiplicity determinant is squared',s.factor(F.jacobian(coords).det()*s.prod(coords)/s.prod(F))==9)
    manifest=json.loads((HERE/'sources.json').read_text())
    for item in manifest['local_sources']:ck('local pin '+item['path'],hashlib.sha256((ROOT/item['path']).read_bytes()).hexdigest()==item['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'standard_module_controls':controls,'block_control':{'group':'S3 wr C2','degree':6,'order':72,'block_stabilizer_order':36,'block_image_order':6},'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact matrix, finite block, small prime-square monodromy and actual differential controls. No formalization of all-endomorphism completeness, cohomological twisting, valuation obstruction, topology, all-degree factorization or canonical lifting. No independent whole-proof replay, external review or novelty certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
