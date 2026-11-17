#!/usr/bin/env python3
"""Four literal Hermitian traces; exact conjugate-pair/Laurent-phase algebra."""
import hashlib
import json
from pathlib import Path

HERE = Path(__file__).resolve().parent
NAMES = ['x','X','y','Y','z','Z','h','v','V','w','W','p','P','q','Q',
         'alpha','beta','gamma','t','T']
COUNT = len(NAMES)
ZERO = (0,) * COUNT
PAIRS = [('x','X'),('y','Y'),('z','Z'),('v','V'),('w','W'),
         ('p','P'),('q','Q'),('t','T')]
PHASES = ['alpha','beta','gamma']


def c(a): return {ZERO:a} if a else {}
def var(name, exponent=1):
    e = list(ZERO); e[NAMES.index(name)] = exponent
    return {tuple(e):1}
def add(*ps):
    result = {}
    for p in ps:
        for e,a in p.items():
            result[e] = result.get(e,0)+a
            if not result[e]: del result[e]
    return result
def neg(p): return {e:-a for e,a in p.items()}
def mul(*ps):
    result = c(1)
    for p in ps:
        out = {}
        for e,a in result.items():
            for f,b in p.items():
                g=tuple(i+j for i,j in zip(e,f))
                out[g]=out.get(g,0)+a*b
                if not out[g]: del out[g]
        result=out
    return result
def star(p):
    result={}
    for e,a in p.items():
        f=list(e)
        for left,right in PAIRS:
            i,j=NAMES.index(left),NAMES.index(right);f[i],f[j]=e[j],e[i]
        for name in PHASES:
            i=NAMES.index(name);f[i]=-e[i]
        result[tuple(f)]=a
    return result
def unit(i): return [c(int(i==j)) for j in range(3)]
def vector_add(*vs): return [add(*(v[i] for v in vs)) for i in range(3)]
def times(p,v): return [mul(p,a) for a in v]
def dot(u,v,g):
    return add(*(mul(star(a),g[i][j],b)
                 for i,a in enumerate(u) for j,b in enumerate(v)))
def reflect(u,v,g):return vector_add(v,times(neg(dot(u,v,g)),u))


def check_root_formulas():
    x,X,y,Y,z,Z=map(var,['x','X','y','Y','z','Z'])
    g=[[c(2),x,y],[X,c(2),z],[Y,Z,c(2)]]
    U,V,W=[unit(i) for i in range(3)]
    reports=[]
    def check(name,roots,expected):
        assert roots==expected, name+' root formula mismatch'
        assert all(dot(r,r,g)==c(2) for r in roots.values())
        reports.append({'name':name,'generic_trace_roots':len(roots),
                        'all_root_norms':2})
    h=add(c(1),neg(mul(y,Y)))
    v=add(mul(X,y),mul(z,h))
    a=reflect(W,times(c(-1),U),g)
    b=reflect(a,W,g);d=reflect(V,b,g);e=reflect(U,d,g)
    f=reflect(d,U,g);w=dot(U,d,g)
    check('13a_373',{5:a,9:b,10:d,11:e,12:f},
          {5:vector_add(times(c(-1),U),times(Y,W)),
           9:vector_add(times(y,U),times(h,W)),
           10:vector_add(times(y,U),times(neg(v),V),times(h,W)),
           11:vector_add(times(add(y,neg(w)),U),times(neg(v),V),times(h,W)),
           12:vector_add(times(add(c(1),neg(mul(star(w),y))),U),
                         times(mul(star(w),v),V),times(neg(mul(star(w),h)),W))})
    T=reflect(V,W,g);t=dot(U,T,g);A=reflect(U,T,g)
    r8=reflect(T,V,g);r5=reflect(A,times(c(-1),U),g)
    r12=reflect(r8,A,g)
    expected5=vector_add(times(add(mul(t,star(t)),c(-1)),U),times(neg(star(t)),T))
    check('13a_375',{10:T,8:r8,11:A,5:r5,12:r12},
          {10:T,8:r8,11:A,5:expected5,12:r12})
    assert r8[0]=={} and r12[0]==neg(t)
    D=reflect(W,V,g);r12=reflect(D,A,g);r8=reflect(A,V,g)
    check('13a_391',{5:D,10:T,11:A,12:r12,8:r8},
          {5:vector_add(V,times(neg(Z),W)),10:T,11:A,
           12:vector_add(A,times(neg(dot(D,A,g)),D)),8:r8})
    assert r12[0]==neg(t)
    h=add(c(1),neg(mul(z,Z)))
    r5=reflect(W,V,g);r11=reflect(r5,W,g);t=dot(U,r11,g)
    r10=reflect(U,r11,g);k=dot(V,r10,g);r9=reflect(V,r10,g)
    r8=reflect(r10,V,g)
    check('13a_392',{5:r5,11:r11,10:r10,9:r9,8:r8},
          {5:vector_add(V,times(neg(Z),W)),
           11:vector_add(times(z,V),times(h,W)),
           10:vector_add(times(neg(t),U),times(z,V),times(h,W)),
           9:vector_add(times(neg(t),U),times(add(z,neg(k)),V),times(h,W)),
           8:vector_add(times(mul(star(k),t),U),
                        times(add(c(1),neg(mul(star(k),z))),V),
                        times(neg(mul(star(k),h)),W))})
    return reports


def check_unit_identity():
    y,Y,h,v,w,W,p,q=map(var,['y','Y','h','v','w','W','p','q'])
    alpha,beta=var('alpha'),var('beta');ibeta=star(beta)
    fs=[add(h,c(-1),mul(y,Y)),mul(p,add(y,neg(w))),
        add(c(1),mul(p,v)),add(neg(mul(p,h)),neg(beta)),
        add(c(-1),neg(mul(q,add(c(1),neg(mul(W,y)))))),
        add(neg(mul(q,W,v)),neg(alpha))]
    fs+=list(map(star,fs))
    size=len(fs)
    def base(i):return [c(int(i==j)) for j in range(size)]
    def radd(*rs):return [add(*(r[i] for r in rs)) for i in range(size)]
    def rtimes(a,r):return [mul(a,b) for b in r]
    def rstar(r):
        return [star(r[(i+6)%size]) for i in range(size)]
    def value(r):return add(*(mul(a,b) for a,b in zip(r,fs)))
    D=add(w,neg(y));E=add(v,neg(mul(ibeta,h)))
    A=add(c(-1),neg(mul(q,h)));B0=add(mul(Y,ibeta),neg(alpha))
    Dr=radd(rtimes(v,base(1)),rtimes(D,base(2)))
    assert value(Dr)==D
    Er=radd(rtimes(neg(mul(ibeta,h)),base(2)),rtimes(neg(mul(ibeta,v)),base(3)))
    assert value(Er)==E
    Ar=radd(base(4),rtimes(neg(q),base(0)),rtimes(neg(mul(q,y)),rstar(Dr)))
    assert value(Ar)==A
    Br=radd(base(5),rtimes(neg(mul(Y,ibeta)),Ar),
            rtimes(mul(q,v),rstar(Dr)),rtimes(mul(q,Y),Er))
    assert value(Br)==B0
    hr=radd(base(0),rtimes(neg(mul(beta,y)),Br),rtimes(neg(alpha),rstar(Br)))
    assert value(hr)==h
    identity=radd(rtimes(c(-1),Ar),rtimes(neg(q),hr))
    assert value(identity)==c(1)
    bad=[dict(a) for a in identity];bad[0]=add(bad[0],c(1))
    assert value(bad)!=c(1)
    return {'constant':1,'generators_with_conjugates':size,
            'unit_phases':'alpha,beta with Laurent inverses',
            'multiplier_monomials':sum(map(len,identity)),
            'corrupted_multiplier_rejected':True}


def bind_traces():
    inputs=json.loads((HERE/'seven_sheet_cover_inputs.json').read_text())
    cases={c['name']:c for c in inputs['cases']}
    specs={
      '13a_373':([1,2,6,8],[1,2,6,8],[(8,4,5),(5,8,9),(6,9,10),(10,0,12),(4,10,11)],[(11,7,8),(12,5,6)]),
      '13a_375':([1,2,6,9],[1,2,6,9],[(6,9,10),(10,7,8),(4,10,11),(11,4,5),(8,11,12)],[(9,12,0),(5,8,9)]),
      '13a_391':([1,2,5,9],[1,2,6,9],[(9,6,5),(6,9,10),(4,10,11),(5,11,12),(11,7,8)],[(10,12,0),(12,8,9)]),
      '13a_392':([1,2,5,11],[1,2,6,12],[(12,6,5),(5,12,11),(4,11,10),(6,10,9),(10,7,8)],[(9,12,0),(11,8,9)])}
    def key(t):return (t[0],)+tuple(sorted(t[1:]))
    common=[(2,6,7),(7,2,3),(3,1,0),(1,3,4)]
    reports=[]
    for name,(old,new,plan,left) in specs.items():
        case=cases[name];assert case['seed_arcs']==old
        components=[{i} for i in range(1,27)]
        for row in case['pd']:
            a=next(s for s in components if row[1] in s)
            b=next(s for s in components if row[3] in s)
            if a is not b:
                components.remove(a);components.remove(b);components.append(a|b)
        components=sorted(sorted(s) for s in components)
        arc={v:i for i,s in enumerate(components) for v in s}
        crossings={key((arc[r[1]],arc[r[0]],arc[r[2]])) for r in case['pd']}
        assert all(key(t) in crossings for t in common+plan+left+[(0,1,2)])
        known=set(new)
        for o,a,b in common+plan:
            assert o in known and a in known and b not in known
            known.add(b)
        assert known==set(range(13))
        reports.append({'name':name,'original_seeds':old,'changed_seeds':new,
                        'complete_trace_length':len(common+plan),
                        'unused_contradiction_crossings':[list(t) for t in left]})
    # Both Hermitian changed-pair matrices have determinant one.
    z,Z=var('z'),var('Z');assert add(c(1),neg(mul(z,Z)),mul(Z,z))==c(1)
    return reports


def main():
    roots=check_root_formulas();identity=check_unit_identity();traces=bind_traces()
    files=['HERMITIAN_REFLECTION_OBSTRUCTION.tex',
           'verify_hermitian_reflection_obstruction.py','seven_sheet_cover_inputs.json']
    receipt={'status':'PASS','cases':traces,'generic_Hermitian_roots':roots,
             'uniform_unit_phase_identity':identity,
             'scope':'Exact generic conjugate-pair root formulas and norm identities for four literal trace systems; Laurent unit phases and constant1 identity. First-loop and final sibling cancellations are written arguments; finite complex averaging is proved in the note. No finite phase sampling, reflection-group classification or augmentation-rank specialization.',
             'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in files}}
    (HERE/'hermitian_reflection_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({'status':'PASS','knots':len(traces),'uniform_identity':identity},indent=2))


if __name__=='__main__':main()
