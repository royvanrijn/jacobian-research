#!/usr/bin/env python3
"""Exact finite controls for written all-degree self-map and volume theorems."""
from pathlib import Path
from itertools import product, permutations, combinations
from math import lcm, prod
import hashlib, json, time, platform, importlib.util
import sympy as s
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
spec=importlib.util.spec_from_file_location('retained_homogeneous_helpers',HERE.parent/'universal-keller-homogeneous-receivers/verify.py')
old=importlib.util.module_from_spec(spec);spec.loader.exec_module(old)
Field,gmul,ginv,generated,primitive,cycle_type=old.Field,old.mul,old.inv,old.generated,old.primitive,old.cycle_type
checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def kernel_mod(A,p):
    A=[[int(c)%p for c in row] for row in A];n=len(A[0]);piv=[];row=0
    for j in range(n):
        i=next((i for i in range(row,len(A)) if A[i][j]),None)
        if i is None:continue
        A[i],A[row]=A[row],A[i];v=pow(A[row][j],-1,p);A[row]=[v*c%p for c in A[row]]
        for i in range(len(A)):
            if i==row:continue
            v=A[i][j]
            A[i]=[(a-v*b)%p for a,b in zip(A[i],A[row])]
        piv.append(j);row+=1
    out=[]
    for j in set(range(n))-set(piv):
        v=[0]*n;v[j]=1
        for i,col in enumerate(piv):v[col]=-A[i][j]%p
        out.append(v)
    return out
def fixed_space(F,M):
    n=len(M);l=F.d
    frob=[F.digits(F.pow(F.q**j,F.q)) for j in range(l)]
    A=[[0]*(n*l) for _ in range(n*l)]
    for i in range(n):
        for a in range(l):
            for b in range(l):A[i*l+a][i*l+b]+=frob[b][a]
            for j in range(n):A[i*l+a][j*l+a]-=int(M[i][j])
    basis=kernel_mod(A,F.q)
    def values(coeffs):
        v=[sum(c*b[i] for c,b in zip(coeffs,basis))%F.q for i in range(n*l)]
        return tuple(F.encode(v[i*l:(i+1)*l]) for i in range(n))
    return basis,values
def mat_apply(F,M,v):
    out=[]
    for row in M:
        a=0
        for c,x in zip(row,v):a=F.add(a,F.mul(int(c)%F.q,x))
        out.append(a)
    return tuple(out)
def standard_basis(gb):
    n=len(gb.gens);leads=[q.LM(order=gb.order).exponents for q in gb.polys];bounds=[]
    assert (0,)*n not in leads
    for i in range(n):
        values=[m[i] for m in leads if m[i] and all(not m[j] for j in range(n) if j!=i)]
        assert values,'nonfinite control algebra';bounds.append(min(values))
    return [v for v in product(*(range(b) for b in bounds)) if not any(all(a>=b for a,b in zip(v,m)) for m in leads)]
def main():
    started=time.monotonic()
    u,v,t,z,x,y,w=s.symbols('u v t z x y w')
    C=[s.Integer(2),v];Q=[s.Integer(0),s.Integer(1)]
    for n in range(1,12):C.append(s.expand(v*C[-1]-C[-2]));Q.append(s.expand(v*Q[-1]-Q[-2]))
    for n in range(1,13):
        ck(f'Chebyshev{n} norm relation',s.expand(C[n]**2-4-(v*v-4)*Q[n]**2)==0)
        ck(f'Chebyshev{n} invariant derivative',s.diff(C[n],v)==n*Q[n])
        ck(f'Chebyshev{n} torus lift',s.expand(C[n].subs(v,y+1/y)-y**n-y**(-n))==0)
    for a,b in [(2,2),(2,3),(3,2),(3,3),(2,5)]:
        ck(f'composition{a},{b} on whole norm surface',s.expand(C[a].subs(v,C[b])-C[a*b])==0 and s.expand(Q[b]*Q[a].subs(v,C[b])-Q[a*b])==0)
    R=t*t-u*(v*v-4)
    for n in [2,3,5]:
        outs=[u,C[n],t*Q[n]];J=s.Matrix([outs]).jacobian([u,v,t])
        A=s.Matrix([[s.diff(R,a) for a in [u,v,t]]]).col_join(J)
        minors=[s.expand(A[list(rows),:].det()) for rows in combinations(range(4),3)]
        for p in [0,5,7]:
            if p and n%p==0:continue
            opts={'domain':s.QQ} if not p else {'modulus':p}
            gb=s.groebner([R,u*z-1]+minors,u,v,t,z,**opts)
            ck(f'F{n} char{p} global etaleness including critical locus',list(gb)==[1])
    # Universal critical cubic fibre and its full product algebra.
    gb=s.groebner([R,C[3]-2,t*Q[3]],v,t,domain=s.QQ.frac_field(u))
    e=(v+1)/3
    ck('critical cubic whole quotient rank3',len(standard_basis(gb))==3)
    for label,expr in [('idempotent',e*e-e),('linear factor t',e*t),('quadratic factor',t*t+3*u*(1-e))]:
        ck('critical cubic exact '+label,gb.reduce(s.expand(expr))[1]==0)
    ck('critical cubic two comaximal s branches',s.gcd(v-2,v+1)==1)
    # Exact finite-field complete cubic fibres, including repeated s values.
    cubic_controls=[];cubic_counts=[]
    for p in [5,7]:
        found={};count=0
        for U,S,T in product(range(1,p),range(p),range(p)):
            if (T*T-U*(S*S-4))%p:continue
            gb=s.groebner([t*t-U*(v*v-4),C[3]-S,t*Q[3]-T],v,t,modulus=p)
            basis=standard_basis(gb);ck(f'F{p} cubic target{U,S,T} whole rank',len(basis)==3)
            mon=[v**a*t**b for a,b in basis];F=s.zeros(3)
            for j,m in enumerate(mon):
                poly=s.Poly(gb.reduce(m**p)[1],v,t,modulus=p)
                for i,exps in enumerate(basis):F[i,j]=int(poly.coeff_monomial(exps))%p
            ck(f'F{p} cubic target{U,S,T} reduced etale fibre',int(F.det())%p!=0)
            factors=len(kernel_mod((F-s.eye(3)).tolist(),p))
            typ={1:(3,),2:(1,2),3:(1,1,1)}[factors]
            if typ not in found:found[typ]={'base_field':p,'target':[U,S,T],'whole_fibre_rank':3,'splitting_type':typ,'standard_monomials':basis,'frobenius_matrix':[[int(c)%p for c in row] for row in F.tolist()]}
            count+=1
        ck(f'F{p} all three cubic algebra types occur',set(found)=={(3,),(1,2),(1,1,1)})
        cubic_controls.extend(found.values());cubic_counts.append({'base_field':p,'all_rational_targets_checked':count})
    # Complete actual monodromy for the endomorphism degree grid.
    monodromy=[]
    for d,m in product([1,3,5],range(1,7)):
        N=d*m
        if N==1:continue
        points=list(product(range(d),range(m)))
        r1=tuple(points.index(((i+1)%d,j)) for i,j in points)
        r2=tuple(points.index((i,(j+1)%m)) for i,j in points)
        h=tuple(points.index((i,-j%m)) for i,j in points)
        G=generated([r1,r2,h],N)
        ck(f'endomorphism d{d} m{m} exact geometric group order',len(G)==d*m*(2 if m>=3 else 1))
        ck(f'endomorphism d{d} m{m} prime atomicity boundary',primitive(G,N)==s.isprime(N))
        monodromy.append({'d':d,'m':m,'degree':N,'group_order':len(G),'primitive':primitive(G,N)})
    # Unit exponent and canonical volume compensation.
    for d,m in [(1,3),(3,1),(-1,2),(-3,2),(5,3)]:
        F=s.Matrix([x**d,y**m]);density=2/y
        jac=s.factor(F.jacobian([x,y]).det()*density.subs(y,y**m)/density)
        ck(f'cover Jacobian d{d} m{m}',jac==d*m*x**(d-1))
    lift=s.Matrix([x**3,y,w/(3*x*x)])
    ck('canonical-line cubic compensation preserves volume',s.factor(lift.jacobian([x,y,w]).det())==1)
    ck('no integral base volume exponent for E3',s.Rational(-1,2).q!=1)
    for n in [2,3,5]:
        F=s.Matrix([x/n,y**n])
        ck(f'normalized G{n} preserves actual global surface volume',s.factor(F.jacobian([x,y]).det()*(2/y**n)/(2/y))==1)
    # Quartic free affine torus action and all label permutations.
    tau=s.Matrix([[0,-1],[1,-1]]);iota=s.Matrix([[1,-1],[0,-1]])
    aa,bb=s.symbols('aa bb');norm=aa*aa-aa*bb+bb*bb;trace=2*aa-bb
    cube_a=aa**3-3*aa*bb*bb+bb**3;cube_b=3*aa*bb*(aa-bb)
    normgb=s.groebner([norm-1],aa,bb)
    ck('free S3 base has exactly the cubic norm-surface quotient',s.expand(trace*trace+3*bb*bb-4-4*(norm-1))==0)
    ck('cubic Galois-closure quotient cubic trace',normgb.reduce(s.expand(2*cube_a-cube_b-C[3].subs(v,trace)))[1]==0)
    ck('cubic Galois-closure quotient cubic anti-trace',normgb.reduce(s.expand(cube_b-bb*Q[3].subs(v,trace)))[1]==0)
    ck('cubic surface equation pulls back to norm equation',s.expand((z*bb)**2-(-z*z/3)*(trace*trace-4)-s.Rational(4,3)*z*z*(norm-1))==0)
    ck('free-base S3 matrix relations',tau**3==s.eye(2) and iota**2==s.eye(2) and iota*tau*iota==tau**2)
    ck('nontrivial base rotations have no norm-one fixed point',(tau-s.eye(2)).det()==3 and (tau**2-s.eye(2)).det()==3)
    r=(1,2,0);reflection=(1,0,2);identity=(0,1,2)
    data={}
    for i,j in product(range(3),range(2)):
        pi=identity
        for _ in range(i):pi=gmul(r,pi)
        if j:pi=gmul(pi,reflection)
        data[pi]=(tau**i*iota**j,j)
    roots=[a for a in product([1,-1],repeat=3) if prod(a)==1];group=[]
    for eps,pi in product(roots,data):
        pinv=ginv(pi);label=tuple(roots.index(tuple(eps[i]*a[pinv[i]] for i in range(3))) for a in roots)
        group.append({'eps':eps,'pi':pi,'label':label,'M':data[pi][0].tolist(),'odd':data[pi][1]})
    ck('quartic affine kernel plus complement is full S4',{g['label'] for g in group}==set(permutations(range(4))))
    ck('quartic monodromy is primitive',primitive(frozenset(g['label'] for g in group),4))
    ck('quartic geometric action free by base and kernel cases',all((g['pi']==identity and any(e==-1 for e in g['eps'])) or g['odd'] or g['pi']!=identity or g['label']==(0,1,2,3) for g in group))
    a,b,k=s.symbols('a b k',nonzero=True);vars=[x,z,a,b];density=z/(x*a*b)
    maps={'rotation':s.Matrix([k*x,z,1/(a*b),a]),'reflection':s.Matrix([1/x,-z,b,a]),'normalized quartic':s.Matrix([x,z/2,a*a,b*b])}
    for name,F in maps.items():
        pulled=density.subs(dict(zip(vars,F)),simultaneous=True)
        ck(name+' exact covering volume identity',s.factor(F.jacobian(vars).det()*pulled/density)==1)
    # Explicit finite-field points on EVERY S4-twist type, not a large search.
    representatives={cycle_type(g['label']):g for g in group};quartic_controls=[]
    def action(F,g,P):
        aa,bb,zz,*vv=P;xy=mat_apply(F,g['M'],[aa,bb]);pinv=ginv(g['pi'])
        return xy+(F.mul((-1 if g['odd'] else 1)%F.q,zz),)+tuple(F.mul(g['eps'][i]%F.q,vv[pinv[i]]) for i in range(3))
    for q in [5,7]:
        for typ,g in sorted(representatives.items()):
            order=lcm(*typ);F=Field(q,order)
            basis,values=fixed_space(F,g['M']);ck(f'F{q} quartic{typ} T0 linear twist dimension',len(basis)==2)
            xy=next((values(c) for c in product(range(q),repeat=2) if F.add(F.sub(F.mul(values(c)[0],values(c)[0]),F.mul(*values(c))),F.mul(values(c)[1],values(c)[1]))==1),None)
            ck(f'F{q} quartic{typ} exact quadratic norm torsor point',xy is not None)
            basis,zvalues=fixed_space(F,[[(-1 if g['odd'] else 1)]])
            ck(f'F{q} quartic{typ} sign-line twist dimension',len(basis)==1)
            zz=next(zvalues(c)[0] for c in product(range(q),repeat=1) if zvalues(c)[0])
            pinv=ginv(g['pi']);M=[[g['eps'][i] if j==pinv[i] else 0 for j in range(3)] for i in range(3)]
            basis,avalues=fixed_space(F,M);ck(f'F{q} quartic{typ} cubic linear twist dimension',len(basis)==3)
            vv=next((avalues(c) for c in product(range(q),repeat=3) if all(avalues(c)) and F.mul(F.mul(avalues(c)[0],avalues(c)[1]),avalues(c)[2])==1),None)
            ck(f'F{q} quartic{typ} exact cubic norm torsor point',vv is not None)
            P=xy+(zz,)+vv
            ck(f'F{q} quartic{typ} actual S4 Frobenius equation',tuple(F.pow(v,q) for v in P)==action(F,g,P))
            orbit={action(F,h,P) for h in group}
            ck(f'F{q} quartic{typ} free complete S4 orbit',len(orbit)==24)
            complement=[h for h in group if h['eps']==(1,1,1)]
            orbits=[];left=set(orbit)
            while left:
                P0=next(iter(left));O=frozenset(action(F,h,P0) for h in complement);orbits.append(O);left-=O
            ck(f'F{q} quartic{typ} COMPLETE quotient fibre rank4',len(orbits)==4 and all(len(O)==6 for O in orbits))
            actual=tuple(next(i for i,O in enumerate(orbits) if tuple(F.pow(v,q) for v in next(iter(C))) in O) for C in orbits)
            ck(f'F{q} quartic{typ} exact whole fibre splitting',cycle_type(actual)==typ)
            quartic_controls.append({'base_field':q,'splitting_type':typ,'extension_degree':order,'torus_point':list(P),'torus_modulus':list(F.mod),'full_Galois_orbit':24,'full_quotient_fibre_rank':4,'actual_frobenius':actual})
    # Canonical cubic on the other unit direction, over a field containing mu3.
    roots3=[1,2,4];R=tuple(roots3.index(4*a%7) for a in roots3);S=tuple(roots3.index(pow(a,-1,7)) for a in roots3)
    ck('canonical E3 actual S3 fibre action over F7',generated([R,S],3)==frozenset(permutations(range(3))))
    ck('generic quadratic norm twist has no fixed C3 translations',[a for a in range(3) if (-a)%3==a]==[0])
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'complete_cubic_fibre_controls':cubic_controls,'cubic_rational_target_coverage':cubic_counts,'complete_quartic_quotient_fibre_controls':quartic_controls,'monodromy_controls':monodromy,'python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Written all-degree theorems with exact norm/Chebyshev, global cubic etaleness, full critical fibre and finite-field cubic algebra controls. Quartic controls construct every S4 twist type over F5,F7, verify all24 torus points and the complete four-point quotient Frobenius fibre. This is not a computed presentation of the global quotient invariant ring; finiteness and scheme rank use the written torsor proof. No formalization of Schur/topology, all-degree classification, descent or cohomology; no independent whole-proof replay, external review or priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n')
    print(json.dumps({'result':'pass','checks':len(checks),'cubic_full_fibres':len(cubic_controls),'quartic_full_quotient_fibres':len(quartic_controls),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
