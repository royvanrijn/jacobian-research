#!/usr/bin/env python3
"""Full Mordell ray Cassels-Tate matrices from exact norm/local certificates.

Run with the repository's SageMath 10.9 / PARI 2.17.3 environment.
Generation searches for norms; verification never invokes a norm solver.
All off-diagonal entries are evaluated, including the coordinate at 3.
The alternating outer-product formula and compressed rank are checked too.
Diagonal zero comes from alternation; skew symmetry is checked afterwards.
Global correctness uses the written proof in MORDELL_RAY_CASSELS_TATE_PAIRING.md.
"""
import argparse
import json
from pathlib import Path
from time import monotonic
from sage.all import *
from itertools import product as cartesian_product
from copy import deepcopy as certificate_deepcopy
import subprocess as norm_subprocess
import sys as worker_sys
from sage.version import version as sage_version

def encode_k(c):
    v=(c.list()+[QQ(0)]*2)[:2]
    return [str(a) for a in v]


def wild_blocks(primes):
    """Zero-sum blocks of size at most three for the cubic character mod 9."""
    groups={0:[],1:[],2:[]}
    for q in sorted(map(ZZ,primes)):
        e={1:0,8:0,2:1,7:1,4:2,5:2}[int(q%9)]
        groups[e].append(q)
    blocks=[[q] for q in groups[0]]
    while groups[1] and groups[2]:blocks.append([groups[1].pop(),groups[2].pop()])
    for e in (1,2):
        assert len(groups[e])%3==0,'Wild norm target is not a cube modulo 9'
        blocks.extend(groups[e][j:j+3] for j in range(0,len(groups[e]),3))
    assert sorted(q for b in blocks for q in b)==sorted(map(ZZ,primes))
    assert all(prod(b)%9 in (1,8) and len(b)<=3 for b in blocks)
    return sorted([sorted(b) for b in blocks],key=lambda b:int(prod(b)))


def swap_norm(K,R,a,el):
    """VBF Theorem 3.1, with its sign and denominator checked exactly."""
    P=K['x'];x=P.gen();G=K.extension(x**3-R,'rho');rho=G.gen()
    assert el.relative_norm()==R
    A,B,C=(el.list()+[K(0)]*3)[:3]
    if not B and not C:
        # R>1 squarefree makes this impossible here; keep the general identity.
        raise ValueError('A scalar reverse witness requires a separate cube case')
    xi=((A*B-C*C*a)-B*rho)/(-(B*B-A*C)-C*rho)
    assert xi.relative_norm()==a
    return xi


def complete_norm(K,R,a,report=lambda *args,**kwargs:None):
    """Unbounded exact fallback; existence proves termination, not a useful bound."""
    P=K['x'];x=P.gen();G=K.extension(x**3-R,'rho');rho=G.gen();s=K.gen()
    H=1
    while True:
        report('exhaustive norm height',H,flush=True)
        for den in range(1,H+1):
            for nums in cartesian_product(range(-H,H+1),repeat=6):
                if max([den]+list(map(abs,nums)))!=H:continue
                if gcd([den]+list(nums))!=1:continue
                A,B,C=[K(nums[j]+nums[j+1]*s) for j in (0,2,4)]
                if A**3+R*B**3+R**2*C**3-3*R*A*B*C==a*den**3:
                    xi=(A+B*rho+C*rho**2)/den
                    assert xi.relative_norm()==a
                    return xi
        H+=1


def norm_worker(method):
    """Isolated input generation; the parent verifies every returned identity."""
    payload=json.loads(worker_sys.stdin.read());R=ZZ(payload['R'])
    K=QuadraticField(6,'s');s=K.gen();P=K['x'];x=P.gen()
    a=K(QQ(payload['a'][0])+QQ(payload['a'][1])*s)
    G=K.extension(x**3-R,'rho')
    assembly=None
    if method=='direct':
        F=G;targets=[a]
    else:
        F=K.extension(x**3-a,'ri')
        if method=='prime-blocks':
            targets=[prod(b) for b in wild_blocks(R.prime_divisors())] if a==5+2*s else list(R.prime_divisors())
        else:targets=[R]
    setup=pari.rnfisnorminit(K.pari_bnf(),F.pari_relative_polynomial(),0)
    elements=[]
    for target in targets:
        candidates=[]
        for flag in ([11,30] if method=='direct' else [0,11,30]):
            out=pari.rnfisnorm(setup,target,flag)
            if out[1]==1:
                el=F(out[0]);assert el.relative_norm()==target
                candidates.append(el)
                if method=='direct' or flag==11:break
            if flag==11 and candidates:break
        if not candidates:raise ValueError('Norm search returned a remainder')
        elements.append(min(candidates,key=lambda el:len(str(el))))
    if method=='direct':xi=elements[0]
    else:
        el=prod(elements);xi=swap_norm(K,R,a,el)
        if method=='prime-blocks':
            assembly=[dict(target=str(t),coefficients=[encode_k(c) for c in (el.list()+[K(0)]*3)[:3]]) for t,el in zip(targets,elements)]
    assert xi.relative_norm()==a
    co=(xi.list()+[K(0)]*3)[:3]
    fac=[[str(q),int(e)] for q,e in ZZ(lcm(c.denominator() for c in co)).factor()]
    print(json.dumps(dict(coefficients=[encode_k(c) for c in co],assembly=assembly,factorization=fac)))


def alternating_form(basis,logs,d3,tame,auxiliary):
    """MCT14: construct an alternating matrix and its compressed congruence form."""
    FF=GF(3);n=len(basis)
    e=vector(FF,[v[0] for v in d3]);f=vector(FF,[v[1] for v in d3])
    alpha=vector(FF,[v[0] for v in logs]);upsilon=vector(FF,[v[1] for v in logs])
    pairs=[(FF(1),e,upsilon),(FF(2),f,alpha)]
    for rec in tame:
        u=vector(FF,n);u[basis.index(ZZ(rec['prime']))]=1
        pairs.append((FF(rec['valuation']),vector(FF,rec['exponents']),u))
    for rec in auxiliary:
        pairs.append((FF(rec['degree'])**(-1),vector(FF,rec['valuations']),vector(FF,rec['residues'])))
    def wedge(u,v):return u.column()*v.row()-v.column()*u.row()
    form=sum((c*wedge(u,v) for c,u,v in pairs),matrix(FF,n))
    columns=matrix(FF,[list(v) for c,u,vv in pairs for v in (u,vv)]).transpose()
    W=columns.column_space().basis_matrix().transpose()
    r=W.ncols();J=matrix(FF,r)
    for c,u,v in pairs:
        J+=c*wedge(W.solve_right(u),W.solve_right(v))
    assert form==W*J*W.transpose() and form.rank()==J.rank()
    compressed=dict(basis=[[int(v) for v in row] for row in W.rows()],
                    matrix=[[int(v) for v in row] for row in J.rows()],rank=int(J.rank()))
    return form,compressed


def evaluate(R,inputs=None,precision=18,norm_timeout=20,quiet=False,norm_method="auto",complete=False):
    R=ZZ(R)
    assert R>1 and R.is_squarefree() and gcd(R,6)==1
    assert all(q%3==2 and kronecker(-2,q)==-1 for q in R.prime_divisors())
    t0=monotonic()
    records=[];tame_vectors=[];auxiliary_vectors=[]
    def report(*args,**kw):
        if not quiet:print(*args,**kw)
    def record(kind,prime,m):
        records.append(dict(kind=kind,prime=str(prime),matrix=[[int(v) for v in row] for row in m.rows()]))
    PR=QQ['x'];x=PR.gen()
    K=QuadraticField(6,'s');s=K.gen();PK=K['x'];xx=PK.gen()
    G=K.extension(xx**3-R,'rho');rhoG=G.gen()
    F=NumberField(x**4-6*x*x+81,'h');h=F.gen()
    d=(h**3+3*h)/18;sf=h-d;k=sf*d/3;z=(-1+d)/2
    assert sf**2==6 and d**2==-3 and k**2==-2 and z*z+z+1==0
    PF=F['y'];y=PF.gen();Mr=F.extension(y**3-R,'rho')
    M=Mr.absolute_field('m');fr,to=M.structure()
    rho=to(Mr.gen());sm=to(Mr(sf));zm=to(Mr(z));km=to(Mr(k))
    assert M.degree()==12 and rho**3==R
    nff=F.pari_nf();nfm=M.pari_nf()
    def emF(v):return to(Mr(v))
    def evK(c,ss):return sum(QQ(a)*ss**i for i,a in enumerate(c.list()))
    def evxi(co,rr,ss):return sum(evK(c,ss)*rr**i for i,c in enumerate(co))
    def vq(a,p):
        if not a:return Infinity
        return QQ(a).valuation(p)
    def inv3deg(n):return ZZ(inverse_mod(n,3))
    def lp(nf,pr,v):return bool(pari.nfislocalpower(nf,pr.pari_prime(),v,3))
    def logres(a,pr):
        pp=pr.smallest_integer();q=pp**pr.residue_class_degree()
        rr=pari.nfmodpr(nff,a,pr.pari_prime())**((q-1)//3)
        zz=pari.nfmodpr(nff,z,pr.pari_prime())
        return next(e for e in range(3) if zz**e==rr)
    def above(fpr):
        gs=fpr.gens()
        return [pr for pr in M.primes_above(fpr.smallest_integer()) if all(emF(g).valuation(pr)>0 for g in gs)]
    f3=next(pr for pr in F.primes_above(3) if (k-4).valuation(pr)>=4)
    assert f3.ramification_index()==2 and f3.residue_class_degree()==1
    m3s=above(f3)
    lam=1-z;eta3=1-lam**3
    assert lp(nff,f3,3/(lam**2*z**2))
    assert lp(nff,f3,4/((1-lam**2)*eta3))
    report('Fields ready',R,monotonic()-t0,flush=True)
    eps=R%9 in [1,8];ps=list(R.prime_divisors());basis=([0] if eps else [])+ps
    reps={}
    for p in basis:
        if p==0:a=5+2*s
        else:
            bb=1
            while True:
                aa=6*bb*bb-p
                if aa>0 and ZZ(aa).is_square():break
                bb+=1
            pi=ZZ(aa).sqrt()+bb*s
            a=p*pi
            for e in range(3):
                aa=a*(5+2*s)**e
                if ZZ(aa[1])%3==0:a=aa;break
            # Unit cubes change only the representative, and keep coefficients small.
            options=[a*(5+2*s)**(3*j) for j in range(-2,3)]
            a=min(options,key=lambda u:max(abs(c) for c in u.list()))
        reps[p]=a
        assert a.norm()==(1 if p==0 else -p**3)
        if p:
            assert ZZ(a[1])%3==0 and ZZ(a[0])%9 in (1,8)
    report('basis',basis,'reps',reps,flush=True)
    afs=[evK(reps[p],sf) for p in basis]
    if inputs is not None:
        assert list(map(int,inputs['basis']))==basis
        assert inputs['representatives']==[encode_k(reps[q]) for q in basis]
    logs=[]
    for af in afs:
        matches=[(a,b) for a,b in cartesian_product(range(3),repeat=2) if lp(nff,f3,af/(z**a*eta3**b))]
        assert len(matches)==1,matches
        logs.append(matches[0])
    report('local a coordinates',logs,flush=True)
    if not eps:assert all(a==0 for a,b in logs)
    xis=[];betas=[];aux=set();encoded_norms=[];factorizations=[];assemblies=[]
    def search_norm(a,method):
        payload=json.dumps(dict(R=int(R),a=encode_k(a)))
        proc=norm_subprocess.run(
            [worker_sys.executable,str(Path(__file__).resolve()),'--norm-worker',method],
            input=payload,text=True,capture_output=True,timeout=norm_timeout)
        if proc.returncode:
            raise ValueError(proc.stderr[-400:])
        data=json.loads(proc.stdout.strip().splitlines()[-1])
        co=[K(QQ(A)+QQ(B)*s) for A,B in data['coefficients']]
        xi=sum(c*rhoG**i for i,c in enumerate(co))
        assert xi.relative_norm()==a
        return xi,data['assembly'],[(ZZ(q),int(e)) for q,e in data['factorization']]
    for row,p in enumerate(basis):
        a=reps[p]
        assembly=None;found_fac=None
        if inputs is None:
            methods=([norm_method] if norm_method!='auto' else
                     ['prime-blocks','direct','reverse'] if R>1000 else ['direct','reverse'])
            xi=None;errors=[]
            for method in methods:
                try:
                    xi,assembly,found_fac=search_norm(a,method)
                    break
                except (norm_subprocess.TimeoutExpired,ValueError) as exc:
                    errors.append(type(exc).__name__+': '+str(exc)[:160])
            if xi is None and complete:xi=complete_norm(K,R,a,report)
            if xi is None:raise RuntimeError('Norm generation incomplete for '+str(p)+': '+str(errors))
            co=(xi.list()+[K(0)]*3)[:3]
        else:
            co=[K(QQ(A)+QQ(B)*s) for A,B in inputs['norm_coefficients'][row]]
            xi=sum(c*rhoG**i for i,c in enumerate(co))
            if 'norm_assembly' in inputs:assembly=inputs['norm_assembly'][row]
        if assembly is not None:
            Fi=K.extension(xx**3-a,'ri');ri=Fi.gen();el=Fi(1);targets=[]
            for block in assembly:
                target=ZZ(block['target']);targets.append(target)
                cc=[K(QQ(A)+QQ(B)*s) for A,B in block['coefficients']]
                assert len(cc)==3
                part=sum(c*ri**j for j,c in enumerate(cc));assert part.relative_norm()==target
                el*=part
            expected=([prod(b) for b in wild_blocks(ps)] if p==0 else ps)
            assert sorted(targets)==sorted(expected)
            assert swap_norm(K,R,a,el).list()==xi.list()
        assemblies.append(assembly)
        assert len(co)==3 and xi.relative_norm()==a
        xis.append(co);encoded_norms.append([encode_k(c) for c in co])
        den=ZZ(lcm(c.denominator() for c in co))
        if inputs is None:fac=found_fac if found_fac is not None else [(q,int(e)) for q,e in den.factor()]
        else:fac=[(ZZ(q),int(e)) for q,e in inputs['denominator_factorizations'][row]]
        assert all(e>0 and q.is_prime(proof=True) for q,e in fac)
        assert prod(q**e for q,e in fac)==den
        factorizations.append([[str(q),e] for q,e in fac])
        aux.update(q for q,e in fac)
        u=evxi(co,rho*zm**2,sm);v=evxi(co,rho*zm,sm);b=u*u*v
        ub=evxi(co,rho*zm**2,-sm);vb=evxi(co,rho*zm,-sm);bbar=ub*ub*vb
        betas.append((b/bbar)**2)
        report('norm',p,'coefficient bytes',len(str(xi)),'aux primes',len(aux),monotonic()-t0,flush=True)
    def sqrtmod3(A,N):
        mod=3
        r=ZZ(1)
        assert A%3==1
        for e in range(1,N):
            r=next(r+mod*j for j in range(3) if ((r+mod*j)**2-A)%(mod*3)==0)
            mod*=3
        return r
    def cubemod3(A,N):
        r=ZZ(1 if A%3==1 else 2)
        assert (r**3-A)%9==0
        mod=3
        for e in range(1,N):
            r=next(r+mod*j for j in range(3) if ((r+mod*j)**3-A)%(mod*9)==0)
            mod*=3
        assert (r**3-A)%(3**(N+1))==0
        return r
    def padd(P,Q):
        if P is None:return Q
        if Q is None:return P
        X,Y=P;U,V=Q
        if X==U:
            if Y==-V:return None
            m=3*X**2/(2*Y)
        else:m=(V-Y)/(U-X)
        x=m*m-X-U
        return (x,m*(X-x)-Y)
    class MorePrecision(Exception):pass
    def three_part(N,supplied):
        assert N>=6
        rc=None;approx_b=[]
        if eps:
            rc=ZZ(supplied['rho3']) if supplied is not None else cubemod3(R,N)
            assert (rc**3-R).valuation(3)>=N+1
            for co in xis:
                vals=[evxi(co,rc*z**2,ss) for ss in [sf,-sf]]+[evxi(co,rc*z,ss) for ss in [sf,-sf]]
                lower=min(2*N+min(evK(c,ss).valuation(f3) for c in co)-v.valuation(f3) for ss,v in [(sf,vals[0]),(-sf,vals[1]),(sf,vals[2]),(-sf,vals[3])])
                if lower<=4:raise MorePrecision(('rho precision',lower))
                approx_b.append((vals[0]**2*vals[2]/(vals[1]**2*vals[3]))**2)
        else:
            assert len(m3s)==1
        def certify_point(P,af):
            X,Y=P
            if not Y:return None
            res=Y*Y-X**3-54*R**2;err=vq(res,3)-vq(Y,3)
            if vq(res,3)<=2*vq(Y,3):return None
            g=Y-3*R*sf
            if 2*err-g.valuation(f3)<=4:return None
            if eps:
                t=Y+sf*d*rc*X+3*R*sf*d
                if min(2*err,2*N+(sf*d*X).valuation(f3))-t.valuation(f3)<=4:return None
            else:
                m3=m3s[0];t=Y+sm*emF(d)*rho*X+3*R*sm*emF(d)
                if err*m3.ramification_index()-t.valuation(m3)<=2*m3.ramification_index():return None
            if not lp(nff,f3,g/af):return None
            return X,Y,err
        localpts=[]
        if supplied is not None:
            assert len(supplied['points3'])==len(basis)
            for af,P in zip(afs,supplied['points3']):
                if P is None:
                    assert lp(nff,f3,af);localpts.append(None)
                else:
                    point=certify_point(tuple(map(QQ,P)),af)
                    assert point is not None,'Invalid local point or insufficient precision'
                    localpts.append(point)
        else:
            yu=sqrtmod3(1+54*R**2,N);Pu=(QQ(1),QQ(yu))
            Pw=(QQ(3),QQ(9)*sqrtmod3((1+2*R**2)//3,N)) if eps else None
            points=[None,Pu,(Pu[0],-Pu[1])]
            if eps:
                points.extend([Pw,(Pw[0],-Pw[1])])
                points.extend(padd((Pu[0],a*Pu[1]),(Pw[0],b*Pw[1])) for a,b in cartesian_product([1,-1],repeat=2))
            shifts=[None]
            for af in afs:
                if lp(nff,f3,af):localpts.append(None);continue
                found=None
                for hshift in range(6):
                    if hshift>=len(shifts):
                        if hshift==1:
                            Q=(QQ(1)/9,QQ(sqrtmod3(1+54*R**2*3**6,N))/27)
                            shift3=padd(padd(Q,Q),Q)
                        shifts.append(padd(shifts[-1],shift3))
                    for P in points:
                        try:candidate=padd(P,shifts[hshift])
                        except ZeroDivisionError:continue
                        if candidate is not None:
                            found=certify_point(candidate,af)
                            if found is not None:break
                    if found is not None:break
                if found is None:raise MorePrecision(('No certified local point',af))
                localpts.append(found)
        d3=[]
        for i,P in enumerate(localpts):
            if eps:
                t=F(1) if P is None else P[1]+sf*d*rc*P[0]+3*R*sf*d
                disc=approx_b[i]/t
                matches=[(e,f) for e,f in cartesian_product(range(3),repeat=2) if lp(nff,f3,disc/(3**e*4**f))]
            else:
                t=M(1) if P is None else P[1]+sm*emF(d)*rho*P[0]+3*R*sm*emF(d)
                disc=betas[i]/t
                matches=[(e,0) for e in range(3) if lp(nfm,m3s[0],disc/(3**e))]
            assert len(matches)==1,('3-descent matches',matches)
            d3.append(matches[0])
        return localpts,d3,rc
    N=int(inputs['precision']) if inputs is not None else int(precision)
    while True:
        try:
            localpts,d3,rc=three_part(N,inputs)
            break
        except MorePrecision as exc:
            if inputs is not None:raise AssertionError(str(exc)) from exc
            report('Increasing local precision',N,'to',2*N,str(exc),flush=True)
            N*=2
    report('local points and descended c at3',d3,monotonic()-t0,flush=True)
    mat=matrix(GF(3),len(basis))
    for i,(e,f) in enumerate(d3):
        for j,(a,b) in enumerate(logs):
            if i!=j:mat[i,j]+=2*(e*b+2*f*a)
    record('three',3,mat)
    report('3 matrix',mat.list(),flush=True)
    for q in ps:
        local=matrix(GF(3),len(basis))
        fpr=F.primes_above(q)[0];mpr=above(fpr)[0]
        assert mpr.ramification_index()==3
        # rho is a uniformizer at every q|R.
        assert rho.valuation(mpr)==1
        kappas=[0]*len(basis)
        for j,p in enumerate(basis):
            if p!=q:continue
            va=afs[j].valuation(fpr)
            assert va%3
            for i,beta in enumerate(betas):
                if i==j:continue
                v=beta.valuation(mpr);assert v%3==0
                unit=beta/rho**v
                rr=pari.nfmodpr(nfm,unit,mpr.pari_prime())**((q**mpr.residue_class_degree()-1)//3)
                zz=pari.nfmodpr(nfm,zm,mpr.pari_prime())
                log=next(e for e in range(3) if zz**e==rr)
                kappas[i]=int(log)
                local[i,j]+=inv3deg(2)*log*va
        tame_vectors.append(dict(prime=str(q),valuation=int(va%3),exponents=kappas))
        mat+=local;record('bad',q,local)
    report('after tame bad primes',mat.list(),flush=True)
    for q in sorted(aux-{2,3}-set(ps)):
        local=matrix(GF(3),len(basis))
        fpr=F.primes_above(q)[0];mpr=above(fpr)[0]
        assert mpr.ramification_index()==1
        logsq=[logres(af,fpr) for af in afs]
        vqis=[b.valuation(mpr)%3 for b in betas]
        deg=fpr.residue_class_degree()
        auxiliary_vectors.append(dict(prime=str(q),degree=int(deg),valuations=list(map(int,vqis)),residues=list(map(int,logsq))))
        for i,v in enumerate(vqis):
            for j,l in enumerate(logsq):
                if i!=j:local[i,j]-=inv3deg(deg)*v*l
        mat+=local;record('auxiliary',q,local)
        if any(vqis) and any(logsq):report('aux contribution',q,vqis,logsq,flush=True)
    report('MATRIX',R,basis,mat,monotonic()-t0,flush=True)
    assert mat+mat.transpose()==0, 'Off-diagonal entries fail skew symmetry'
    outer,compressed=alternating_form(basis,logs,d3,tame_vectors,auxiliary_vectors)
    assert outer==mat,'Alternating outer-product formula differs from direct entries'
    result=dict(R=int(R),basis=[int(p) for p in basis],representatives=[encode_k(reps[q]) for q in basis],norm_coefficients=encoded_norms,denominator_factorizations=factorizations,precision=N,points3=[None if P is None else [str(P[0]),str(P[1])] for P in localpts],rho3=str(rc) if eps else None,local_basis3=[list(map(int,v)) for v in logs],local_descent3=[list(map(int,v)) for v in d3],contributions=records,matrix=[[int(v) for v in row] for row in mat.rows()])
    if any(a is not None for a in assemblies) or (inputs is not None and 'norm_assembly' in inputs):result['norm_assembly']=assemblies
    result['local_vectors']=dict(tame=tame_vectors,auxiliary=auxiliary_vectors,compressed=compressed)
    if inputs is not None:
        expected=certificate_deepcopy(inputs)
        expected.setdefault('local_vectors',result['local_vectors'])
        assert result==expected,'Certificate output does not match exact replay'
    return result

def compare_controls(cases):
    path=Path(__file__).with_name('mordell_ray_pairing_exact.json')
    old=json.loads(path.read_text())
    count=0
    for case in cases:
        control=next((c for c in old['curves'] if c['R']==case['R']),None)
        if control is None:continue
        order=control.get('basis',control['primes'])
        perm=[case['basis'].index(p) for p in order]
        got=[[case['matrix'][i][j] for j in perm] for i in perm]
        assert got==control['normalized_pairing_matrix'],case['R']
        count+=1
    print(f'Compared {count} complete matrices with the separately constructed controls.')

def negative_controls(cases):
    K=QuadraticField(6,'s')
    assert complete_norm(K,71,K(71)).relative_norm()==71
    assert wild_blocks([5,71,29])==[[71],[5,29]]
    assert wild_blocks([5,23,239])==[[5,23,239]]
    print('Checked exhaustive-fallback witness and both wild block types')
    case=next(c for c in cases if c['R']==71)
    mutations=[]
    x=certificate_deepcopy(case)
    x['norm_coefficients'][0][0][0]=str(QQ(x['norm_coefficients'][0][0][0])+1)
    mutations.append(('corrupted norm',x))
    x=certificate_deepcopy(case)
    i=next(i for i,v in enumerate(x['points3']) if v is not None)
    x['points3'][i][1]=str(QQ(x['points3'][i][1])+1)
    mutations.append(('invalid Hensel point',x))
    x=certificate_deepcopy(case)
    x['contributions']=[v for v in x['contributions'] if v['kind']!='auxiliary']
    mutations.append(('omitted auxiliary prime',x))
    x=certificate_deepcopy(case);x['matrix'][0][1]=(x['matrix'][0][1]+1)%3
    mutations.append(('corrupted matrix',x))
    blockcase=next((c for c in cases if 'norm_assembly' in c),None)
    if blockcase is not None:
        x=certificate_deepcopy(blockcase)
        block=next(a for a in x['norm_assembly'] if a)[0]
        block['coefficients'][0][0]=str(QQ(block['coefficients'][0][0])+1)
        mutations.append(('corrupted prime-block assembly',x))
    for label,x in mutations:
        try:evaluate(x['R'],inputs=x,quiet=True)
        except AssertionError:print('Rejected '+label)
        else:raise AssertionError('Accepted '+label)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    action=parser.add_mutually_exclusive_group(required=True)
    action.add_argument('--generate',nargs='+',type=int,metavar='R')
    action.add_argument('--verify',type=Path)
    action.add_argument('--norm-worker',choices=['direct','reverse','prime-blocks'])
    parser.add_argument('--output',type=Path)
    parser.add_argument('--precision',type=int,default=18)
    parser.add_argument('--norm-timeout',type=int,default=20)
    parser.add_argument('--norm-method',choices=['auto','direct','reverse','prime-blocks'],default='auto')
    parser.add_argument('--complete',action='store_true',help='Unbounded exhaustive norm fallback after bounded searches fail')
    parser.add_argument('--compare-controls',action='store_true')
    parser.add_argument('--self-test',action='store_true')
    args=parser.parse_args()
    pari.allocatemem(500000000)
    if args.norm_worker:
        norm_worker(args.norm_worker)
        return
    if args.generate:
        if args.output is None:parser.error('--generate requires --output')
        cases=[evaluate(R,precision=args.precision,norm_timeout=args.norm_timeout,norm_method=args.norm_method,complete=args.complete)
               for R in args.generate]
        doc=dict(schema_version=1,
                 scope='Exact full local-symbol matrices conditional on the cited global descent proof; no rank or rational-point oracle is used to compute entries.',
                 software=dict(sage=sage_version,pari=list(map(int,pari.version()))),
                 cases=cases)
        args.output.write_text(json.dumps(doc,indent=2)+'\n')
    else:
        doc=json.loads(args.verify.read_text())
        assert doc['schema_version']==1
        cases=[evaluate(c['R'],inputs=c) for c in doc['cases']]
    if args.compare_controls:compare_controls(cases)
    if args.self_test:negative_controls(cases)
    print('PASS: full off-diagonal local-symbol evaluation; exact norms, Hensel bounds, '
          'local cube descent, all auxiliary primes, alternating formula and compressed rank checked. Written global proof '
          'and common Weil sign remain separate.')

if __name__=='__main__':
    main()
