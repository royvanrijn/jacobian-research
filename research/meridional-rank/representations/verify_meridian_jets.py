#!/usr/bin/env python3
"""Replay finite jets with entrywise polynomial adjugates and all minors.

Does not use the producer's recursive matrix-series inverse or its
selected-minor equations. Positive finite jets give no char0 curve claim.
"""
import hashlib
import json
from pathlib import Path
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent
P=11;D=3
Z=(0,)*4;U=(1,0,0,0)
I=tuple(U if i==j else Z for i in range(3) for j in range(3))


def pa(a,b):return tuple((x+y)%P for x,y in zip(a,b))
def pn(a):return tuple(-x%P for x in a)
def pm(a,b):return tuple(sum(a[j]*b[k-j] for j in range(k+1))%P for k in range(4))
def pi(a):
    assert a[0]%P
    # Geometric series around the nonzero constant, rather than matrix recurrence.
    c=pow(a[0],-1,P);u=tuple((-c*x)%P if j else 0 for j,x in enumerate(a))
    out=U;power=U
    for _ in range(D):power=pm(power,u);out=pa(out,power)
    return tuple(c*x%P for x in out)
def ps(xs):
    out=Z
    for x in xs:out=pa(out,x)
    return out
def mm(a,b):return tuple(ps(pm(a[3*i+k],b[3*k+j]) for k in range(3)) for i in range(3) for j in range(3))
def mi(a):
    cof=[]
    for i in range(3):
        for j in range(3):
            r=[k for k in range(3) if k!=i];c=[k for k in range(3) if k!=j]
            v=pa(pm(a[3*r[0]+c[0]],a[3*r[1]+c[1]]),pn(pm(a[3*r[0]+c[1]],a[3*r[1]+c[0]])))
            cof.append(v if (i+j)%2==0 else pn(v))
    det=ps(pm(a[j],cof[j]) for j in range(3));invdet=pi(det)
    out=tuple(pm(cof[3*j+i],invdet) for i in range(3) for j in range(3))
    assert mm(a,out)==I
    return out
def word(ms,letters):
    out=I
    for l in letters:out=mm(out,ms[l-1] if l>0 else mi(ms[-l-1]))
    return out


def one(source,row,pd):
    variables=row['variables'];tau=row['tau']
    assert row['prime']==11 and tau in [4,8]
    assert row['first_obstruction'] is None and [r['degree'] for r in row['orders']]==[1,2,3]
    assert all(not any(r['residual_after']) for r in row['orders'])
    assert variables[0]==[0,1,0,0]
    assert all(not any(variables[j]) for j in row['gauge_coordinates'])
    base=[tuple((a+tau*b)%P for a,b in m) for m in source['matrices']]
    ms=[None]*len(base);seeds=source['seeds'];m=base[seeds[0]]
    ms[seeds[0]]=tuple(tuple((m[i] if k==0 else 0)+(int(i//3==i%3)-m[i])*variables[0][k] for k in range(4)) for i in range(9))
    ms[seeds[0]]=tuple(tuple(v%P for v in e) for e in ms[seeds[0]])
    for j,s in enumerate(seeds[1:]):
        ms[s]=tuple(tuple(((base[s][i] if k==0 else 0)+variables[1+9*j+i][k])%P for k in range(4)) for i in range(9))
    crossings=oriented_diagram(pd)
    for o,a,b in source['trace']:
        c=next(c for c in crossings if c[0]==o and {c[1],c[2]}=={a,b})
        sign=c[3]^(a!=c[1]);g=mi(ms[o]) if sign else ms[o]
        ms[b]=mm(mm(g,ms[a]),mi(g))
    for o,a,b,sign in crossings:
        rel=[-(o+1),a+1,o+1,-(b+1)] if sign else [o+1,a+1,-(o+1),-(b+1)]
        assert word(ms,rel)==I
    for m in ms:
        r=tuple(pa(a,pn(b)) for a,b in zip(m,I))
        for i in range(3):
            for j in range(i+1,3):
                for k in range(3):
                    for l in range(k+1,3):
                        assert pa(pm(r[3*i+k],r[3*j+l]),pn(pm(r[3*i+l],r[3*j+k])))==Z
    # The selected derivative matrix is nonsingular in the stored finite model.
    a=[row['jacobian'][i] for i in row['selected_rows']];b=row['selected_inverse']
    assert len(a)==len(b)==28
    assert all(sum(a[i][k]*b[k][j] for k in range(28))%P==int(i==j) for i in range(28) for j in range(28))
    return {'name':row['name'],'orbit_index':row['orbit_index'],'tau':tau,'degree':3,
            'oriented_relators':len(crossings),'all_rank_one_minors':True}


def main():
    jets=json.loads((HERE/'meridian_jet_results.json').read_text())
    assert jets['status']=='COMPLETE' and len(jets['results'])==32
    for p,h in jets['input_sha256'].items():assert hashlib.sha256((HERE/p).read_bytes()).hexdigest()==h
    source=json.loads((HERE/'meridian_tangent_results.json').read_text())
    index={(r['name'],r['orbit_index']):r for r in source['results']}
    pd={r['name']:r['pd'] for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    rows=[one(index[(r['name'],r['orbit_index'])],r,pd[r['name']]) for r in jets['results']]
    paths=['verify_meridian_jets.py','meridian_jet_results.json','meridian_jet_plan.json','parabolic_lifts.py']
    result={'status':'PASS','jets':32,'bases':16,'prime':11,'maximum_degree':3,'results':rows,
            'conclusion':'All frozen degree3 jets satisfy every oriented relator and all meridian rank-one minors modulo11 in both embeddings. No finite-order obstruction found in this test; no actual char0 curve or all-order integrability claim.',
            'replay':'Entrywise polynomial adjugates, scalar geometric-series inverses and all minors; shared source PD parser and inputs',
            'input_sha256':{p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}}
    (HERE/'meridian_jet_verification.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['status','jets','maximum_degree','conclusion']},indent=2))


if __name__=='__main__':main()
