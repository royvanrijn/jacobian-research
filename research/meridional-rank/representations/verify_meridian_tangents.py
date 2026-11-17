#!/usr/bin/env python3
"""Independent rational dual-number/minor and restricted-cocycle replay.

No Sympy, seed derivative formula, or producer nullspace solver is used.
The retained H3 coefficient model and source PD parser are shared inputs.
"""
from fractions import Fraction as F
from functools import lru_cache
import hashlib
import json
from pathlib import Path
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent
Z=(F(0),F(0)); U=(F(1),F(0))
I=tuple(U if i==j else Z for i in range(3) for j in range(3))
N=(Z,)*9


def qa(x,y): return x[0]+y[0],x[1]+y[1]
def qn(x): return -x[0],-x[1]
def qm(x,y): return x[0]*y[0]+x[1]*y[1],x[0]*y[1]+x[1]*y[0]+x[1]*y[1]
def qi(x):
    a,b=x; norm=a*a+a*b-b*b
    assert norm
    return (a+b)/norm,-b/norm
def qs(xs):
    out=Z
    for x in xs: out=qa(out,x)
    return out
def am(a,b): return tuple(qa(x,y) for x,y in zip(a,b))
def nm(a): return tuple(qn(x) for x in a)
def sm(c,a): return tuple(qm(c,x) for x in a)
def mm(a,b):
    return tuple(qs(qm(a[3*i+k],b[3*k+j]) for k in range(3)) for i in range(3) for j in range(3))
def mv(a,v): return tuple(qs(qm(a[3*i+j],v[j]) for j in range(3)) for i in range(3))
def decode(x): return F(*x[0]),F(*x[1])
def read_matrix(a): return tuple((F(x),F(y)) for x,y in a)


def rref(rows):
    a=[list(r) for r in rows]; k=0; pivots=[]
    for j in range(len(a[0])):
        i=next((i for i in range(k,len(a)) if a[i][j]!=Z),None)
        if i is None: continue
        a[k],a[i]=a[i],a[k]
        c=qi(a[k][j]); a[k]=[qm(c,x) for x in a[k]]
        for i in range(len(a)):
            if i!=k and a[i][j]!=Z:
                c=a[i][j]; a[i]=[qa(x,qn(qm(c,y))) for x,y in zip(a[i],a[k])]
        pivots.append(j); k+=1
        if k==len(a): break
    return a,pivots


def qrank(rows): return len(rref(rows)[1])
def kernel(rows):
    a,pivots=rref(rows); result=[]
    for j in range(len(rows[0])):
        if j in pivots: continue
        v=[Z]*len(rows[0]);v[j]=U
        for i,p in enumerate(pivots):v[p]=qn(a[i][j])
        result.append(v)
    return result


def mod_rank(rows,p,t):
    a=[[(int(x)+int(y)*t)%p for x,y in r] for r in rows];k=0
    for j in range(len(a[0])):
        i=next((i for i in range(k,len(a)) if a[i][j]),None)
        if i is None: continue
        a[k],a[i]=a[i],a[k]; c=pow(a[k][j],-1,p)
        a[k]=[x*c%p for x in a[k]]
        for i in range(k+1,len(a)):
            c=a[i][j]
            a[i]=[(x-c*y)%p for x,y in zip(a[i],a[k])]
        k+=1
        if k==len(a):break
    return k


def certify(rows,cert):
    basis=[[decode(x) for x in r] for r in cert['basis']]
    assert qrank(basis)==len(basis)==cert['nullity']
    for v in basis:
        assert all(qs(qm(tuple(x),y) for x,y in zip(r,v))==Z for r in rows)
    assert cert['exact_rank']+len(basis)==cert['columns']==len(rows[0])
    for item in cert['modular_rank_certificates']:
        p,t=item['prime'],item['tau']
        assert (t*t-t-1)%p==0 and p in [11,19,29,31]
        assert mod_rank(rows,p,t)==item['rank']
    assert max(r['rank'] for r in cert['modular_rank_certificates'])==cert['exact_rank']
    return basis


def dp(a,b): return mm(a[0],b[0]),am(mm(a[1],b[0]),mm(a[0],b[1]))
@lru_cache(maxsize=None)
def inverse(a):
    if mm(a,a)==I:return a
    rows=[list(a[3*i:3*i+3])+[U if i==j else Z for j in range(3)] for i in range(3)]
    reduced,pivots=rref(rows)
    assert pivots==[0,1,2]
    return tuple(reduced[i][3+j] for i in range(3) for j in range(3))
def di(a):
    inv=inverse(a[0])
    return inv,nm(mm(mm(inv,a[1]),inv))
def word(images,letters):
    out=(I,N)
    for l in letters:out=dp(out,images[l-1] if l>0 else di(images[-l-1]))
    return out


def extend(row,v,crossings):
    ms=[read_matrix(m) for m in row['matrices']];ds=[None]*len(ms)
    seeds=row['seeds'];ds[seeds[0]]=sm(v[0],am(I,nm(ms[seeds[0]])))
    for j,s in enumerate(seeds[1:]): ds[s]=tuple(v[1+9*j:1+9*(j+1)])
    for o,a,b in row['trace']:
        choices=[c for c in crossings if c[0]==o and {c[1],c[2]}=={a,b}]
        c=choices[0];sign=c[3]^(a!=c[1])
        oo=(ms[o],ds[o]);aa=(ms[a],ds[a])
        result=dp(dp(di(oo) if sign else oo,aa),oo if sign else di(oo))
        assert result[0]==ms[b]
        ds[b]=result[1]
    return ms,ds


def minors(m,d):
    a=am(m,nm(I))
    for i in range(3):
        for j in range(i+1,3):
            for k in range(3):
                for l in range(k+1,3):
                    base=qa(qm(a[3*i+k],a[3*j+l]),qn(qm(a[3*i+l],a[3*j+k])))
                    deriv=qs([qm(d[3*i+k],a[3*j+l]),qm(a[3*i+k],d[3*j+l]),
                              qn(qm(d[3*i+l],a[3*j+k])),qn(qm(a[3*i+l],d[3*j+k]))])
                    assert base==deriv==Z


def tangent_replay(row,pd):
    cert=row['certificate'];basis=certify(row['matrix'],cert)
    crossings=oriented_diagram(pd)
    for v in basis:
        ms,ds=extend(row,v,crossings)
        for m,d in zip(ms,ds):
            assert mm(m,m)==I;minors(m,d)
        for o,a,b,sign in crossings:
            rel=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
            base,deriv=word(list(zip(ms,ds)),rel)
            assert base==I and deriv==N
    # Rebuild each retained crossing block via dual-number relator evaluation.
    for j in range(28):
        v=[U if k==j else Z for k in range(28)]
        ms,ds=extend(row,v,crossings)
        for k,(o,a,b,sign) in enumerate(crossings):
            rel=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
            _,deriv=word(list(zip(ms,ds)),rel)
            stored=tuple(tuple(row['matrix'][9*k+i][j]) for i in range(9))
            assert deriv==nm(mm(stored,ms[b]))
        for k,s in enumerate(row['seeds'][1:]):
            plus=am(I,ms[s]);block=mm(mm(plus,ds[s]),plus)
            assert block==tuple(tuple(row['matrix'][9*len(crossings)+9*k+i][j]) for i in range(9))
    # Centralizer of M0 has dimension5, with scalar kernel in its conjugation action.
    m0=ms[row['seeds'][0]]
    units=[tuple(U if k==j else Z for k in range(9)) for j in range(9)]
    comms=[am(mm(e,m0),nm(mm(m0,e))) for e in units]
    centralizer=kernel(list(zip(*comms)))
    assert len(centralizer)==5
    gauge=[]
    for a in centralizer:
        v=[Z]
        for s in row['seeds'][1:]:v.extend(am(mm(a,ms[s]),nm(mm(ms[s],a))))
        gauge.append(v)
        assert all(qs(qm(tuple(x),y) for x,y in zip(r,v))==Z for r in row['matrix'])
    gauge_rank=qrank(gauge)
    assert gauge_rank==4 and len(basis)==5
    assert any(v[0]!=Z for v in basis) and cert['eigenvalue_derivative_available']
    return {'name':row['name'],'orbit_index':row['orbit_index'],'normalized_tangent_dimension':5,
            'residual_conjugation_dimension':4,'modulo_conjugation_dimension':1,
            'nonzero_meridian_eigenvalue_derivative':True}


def restricted_replay(source,module,pd):
    images=[read_matrix(m) for m in source['matrices']]
    if module['module']=='Wdual':images=[tuple(m[3*j+i] for i in range(3) for j in range(3)) for m in images]
    roots=[tuple(tuple(x) for x in r) for r in module['roots']]
    basis=certify(module['matrix'],module['certificate'])
    assert len(basis)==3
    crossings=oriented_diagram(pd)
    # Direct affine multiplication, on the full oriented knot relators.
    def cocycle(v,rel):
        prefix=I;u=(Z,)*3
        for l in rel:
            j=abs(l)-1;ui=tuple(qm(v[j],x) for x in roots[j])
            if l<0:ui=tuple(qn(x) for x in mv(images[j],ui))
            u=tuple(qa(x,y) for x,y in zip(u,mv(prefix,ui)))
            prefix=mm(prefix,images[j])
        assert prefix==I
        return u
    for m,r in zip(images,roots):assert mv(m,r)==tuple(qn(x) for x in r)
    for v in basis:
        for o,a,b,sign in crossings:
            rel=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
            assert cocycle(v,rel)==(Z,)*3
    for j in range(len(images)):
        v=[U if i==j else Z for i in range(len(images))]
        for k,(o,a,b,sign) in enumerate(crossings):
            rel=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
            assert cocycle(v,rel)==tuple(tuple(module['matrix'][3*k+i][j]) for i in range(3))
    principal=[]
    for j in range(3):
        v=tuple(U if j==i else Z for i in range(3));scalars=[]
        for m,r in zip(images,roots):
            u=tuple(qa(x,qn(y)) for x,y in zip(v,mv(m,v)))
            k=next(k for k in range(3) if r[k]!=Z);c=qm(u[k],qi(r[k]))
            assert u==tuple(qm(c,x) for x in r);scalars.append(c)
        principal.append(scalars)
        assert all(qs(qm(tuple(x),y) for x,y in zip(r,scalars))==Z for r in module['matrix'])
    assert qrank(principal)==3
    return {'module':module['module'],'restricted_cocycle_dimension':3,'principal_dimension':3,'restricted_H1':0}


def replay_controls(controls):
    ms=[read_matrix(m) for m in controls['torus_matrices']]
    meridian,_=word(list(zip(ms,[N,N])),[-1,2,2])
    assert meridian==read_matrix(controls['meridian'])
    columns=[]
    for j in range(19):
        ds=[tuple(U if j==1+9*i+k else Z for k in range(9)) for i in range(2)]
        base,dr=word(list(zip(ms,ds)),[1,1,1,-2,-2,-2,-2,-2]);assert base==I
        _,dm=word(list(zip(ms,ds)),[-1,2,2])
        target=am(I,nm(meridian)) if j==0 else N
        columns.append(list(dr)+list(am(dm,nm(target))))
    matrix=list(zip(*columns))
    varying=certify(matrix,controls['varying_T35'])
    fixed=certify(matrix+[[U]+[Z]*18],controls['involutive_T35'])
    assert len(varying)==5 and len(fixed)==4
    assert any(v[0]!=Z for v in varying) and all(v[0]==Z for v in fixed)
    return {'T35_varying_nullity':5,'T35_involutive_nullity':4,'independent_control_matrix_rebuild':True}


def main():
    tangent=json.loads((HERE/'meridian_tangent_results.json').read_text())
    stable=json.loads((HERE/'stabilization_results.json').read_text())
    assert tangent['status']==stable['status']=='COMPLETE'
    for report in [tangent,stable]:
        for p,h in report['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    panel={r['name']:r for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    source={(r['name'],r['orbit_index']):r for r in tangent['results']}
    tangents=[tangent_replay(r,panel[r['name']]['pd']) for r in tangent['results']]
    modules=[]
    for row in stable['results']:
        src=source[(row['name'],row['orbit_index'])]
        modules.append({'name':row['name'],'orbit_index':row['orbit_index'],
                        'modules':[restricted_replay(src,m,panel[row['name']]['pd']) for m in row['modules']]})
    controls=replay_controls(tangent['controls'])
    paths=['verify_meridian_tangents.py','meridian_tangent_results.json','stabilization_results.json',
           'meridian_tangent_plan.json','stabilization_plan.json','parabolic_panel.json','parabolic_lifts.py']
    result={'status':'PASS','tangent_bases':16,'restricted_modules':32,'tangents':tangents,
            'restricted_cocycles':modules,'controls':controls,
            'assurance':'Independent standard-library rational dual-number/minor replay and exact kernel/rank/gauge certificates; shared source PD parser and retained H3 coefficient model. No actual curve, global exclusion or new sharp MRC case.',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'meridian_tangent_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','tangent_bases','restricted_modules','controls','assurance']},indent=2))


if __name__=='__main__':main()
