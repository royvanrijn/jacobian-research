#!/usr/bin/env python3
"""Frozen modular implicit-curve jets at all retained H3 points."""
import hashlib
import json
from pathlib import Path
import time
from parabolic_lifts import oriented_diagram

HERE=Path(__file__).resolve().parent
P=11
I=tuple(int(i==j) for i in range(3) for j in range(3));Z=(0,)*9


def mm(a,b):
    return tuple(sum(a[3*i+k]*b[3*k+j] for k in range(3))%P for i in range(3) for j in range(3))
def ma(a,b):return tuple((x+y)%P for x,y in zip(a,b))
def mn(a):return tuple(-x%P for x in a)
def sp(a,b,n):return [sum(a[i]*b[k-i] for i in range(k+1))%P for k in range(n+1)]
def mp(a,b,n):
    return [tuple(sum(mm(a[i],b[k-i])[j] for i in range(k+1))%P for j in range(9)) for k in range(n+1)]
def mi(a,n):
    # All base meridians are involutions.
    out=[a[0]]
    assert mm(a[0],a[0])==I
    for k in range(1,n+1):
        c=Z
        for i in range(1,k+1):c=ma(c,mm(a[i],out[k-i]))
        out.append(mn(mm(a[0],c)))
    return out


def echelon(rows):
    a=[list(r) for r in rows];pivots=[];k=0
    for j in range(len(a[0])):
        i=next((i for i in range(k,len(a)) if a[i][j]%P),None)
        if i is None:continue
        a[k],a[i]=a[i],a[k];c=pow(a[k][j],-1,P)
        a[k]=[x*c%P for x in a[k]]
        for i in range(len(a)):
            if i!=k:
                c=a[i][j];a[i]=[(x-c*y)%P for x,y in zip(a[i],a[k])]
        pivots.append(j);k+=1
        if k==len(a):break
    return a,pivots
def kernel(rows):
    a,pivots=echelon(rows);out=[]
    for j in range(len(rows[0])):
        if j in pivots:continue
        v=[0]*len(rows[0]);v[j]=1
        for i,p in enumerate(pivots):v[p]=-a[i][j]%P
        out.append(v)
    return out
def invert(a):
    n=len(a);aug=[list(r)+[int(i==j) for j in range(n)] for i,r in enumerate(a)]
    reduced,pivots=echelon(aug);assert pivots==list(range(n))
    return [r[n:] for r in reduced]


def evaluate(source,ms,crossings,gauge,pivots,variables,n):
    seeds=source['seeds'];values=[None]*len(ms)
    values[seeds[0]]=[ma(ms[seeds[0]],tuple(v*variables[0][0]%P for v in ma(I,mn(ms[seeds[0]]))))]
    values[seeds[0]] += [tuple(v*variables[0][k]%P for v in ma(I,mn(ms[seeds[0]]))) for k in range(1,n+1)]
    for j,s in enumerate(seeds[1:]):
        values[s]=[tuple((ms[s][i] if k==0 else 0)+variables[1+9*j+i][k] for i in range(9)) for k in range(n+1)]
        values[s]=[tuple(x%P for x in r) for r in values[s]]
    for o,a,b in source['trace']:
        row=next(c for c in crossings if c[0]==o and {c[1],c[2]}=={a,b})
        sign=row[3]^(a!=row[1]);oo=values[o];inv=mi(oo,n)
        values[b]=mp(mp(inv if sign else oo,values[a],n),oo if sign else inv,n)
    constraints=[]
    for o,a,b,sign in crossings:
        oo=values[o];inv=mi(oo,n)
        expected=mp(mp(inv if sign else oo,values[a],n),oo if sign else inv,n)
        constraints.extend([[(values[b][k][i]-expected[k][i])%P for k in range(n+1)] for i in range(9)])
    for s,(a,b) in zip(seeds[1:],pivots):
        r=[[((values[s][k][i]-I[i]) if k==0 else values[s][k][i])%P for k in range(n+1)] for i in range(9)]
        for i in range(3):
            for j in range(3):
                if i==a or j==b:continue
                left,right=sp(r[3*i+j],r[3*a+b],n),sp(r[3*i+b],r[3*a+j],n)
                constraints.append([(x-y)%P for x,y in zip(left,right)])
    constraints.extend([variables[j][:n+1] for j in gauge])
    s=list(variables[0][:n+1])
    if n:s[1]=(s[1]-1)%P
    constraints.append(s)
    return constraints


def one(source,pd,tau):
    stamp=time.monotonic()
    ms=[tuple((a+tau*b)%P for a,b in m) for m in source['matrices']]
    crossings=oriented_diagram(pd)
    fixed=kernel([[(a+tau*b)%P for a,b in r] for r in source['matrix']]+[[1]+[0]*27])
    assert len(fixed)==4
    _,gauge=echelon(fixed);assert len(gauge)==4 and 0 not in gauge
    pivots=[]
    for s in source['seeds'][1:]:
        p=next(i for i in range(9) if (ms[s][i]-I[i])%P)
        pivots.append(divmod(p,3))
    columns=[]
    zero=[[0]*2 for _ in range(28)]
    baseline=evaluate(source,ms,crossings,gauge,pivots,zero,1)
    for j in range(28):
        v=[[0]*2 for _ in range(28)];v[j][1]=1
        out=evaluate(source,ms,crossings,gauge,pivots,v,1)
        columns.append([(a[1]-b[1])%P for a,b in zip(out,baseline)])
    jac=list(map(list,zip(*columns)))
    _,selected=echelon(list(zip(*jac)));assert len(selected)==28
    inverse=invert([jac[i] for i in selected])
    variables=[[0]*4 for _ in range(28)];orders=[]
    for degree in range(1,4):
        residual=[r[degree] for r in evaluate(source,ms,crossings,gauge,pivots,variables,degree)]
        rhs=[-residual[i]%P for i in selected]
        correction=[sum(a*b for a,b in zip(r,rhs))%P for r in inverse]
        remaining=[(b+sum(a*x for a,x in zip(r,correction)))%P for r,b in zip(jac,residual)]
        orders.append({'degree':degree,'residual_before':residual,'correction':correction,
                       'residual_after':remaining})
        for j,x in enumerate(correction):variables[j][degree]=x
        exact=[r[degree] for r in evaluate(source,ms,crossings,gauge,pivots,variables,degree)]
        assert exact==remaining and all(exact[i]==0 for i in selected)
        if any(remaining):break
    return {'name':source['name'],'orbit_index':source['orbit_index'],'prime':P,'tau':tau,
            'gauge_coordinates':gauge,'minor_pivots':pivots,'jacobian':jac,'selected_rows':selected,
            'selected_inverse':inverse,'orders':orders,'variables':variables,
            'first_obstruction':next((r['degree'] for r in orders if any(r['residual_after'])),None),
            'seconds':time.monotonic()-stamp}


def main():
    plan=json.loads((HERE/'meridian_jet_plan.json').read_text())
    assert plan['frozen_before_execution']
    source=json.loads((HERE/'meridian_tangent_results.json').read_text())
    panel={r['name']:r for r in json.loads((HERE/'parabolic_panel.json').read_text())['cases']}
    paths=['audit_meridian_jets.py','meridian_jet_plan.json','meridian_tangent_results.json',
           'meridian_tangent_verification.json','parabolic_panel.json','parabolic_lifts.py']
    hashes={p:hashlib.sha256((HERE/p).read_bytes()).hexdigest() for p in paths}
    output=HERE/'meridian_jet_results.json'
    report=json.loads(output.read_text()) if output.exists() else {'status':'RUNNING','results':[],'input_sha256':hashes}
    assert report['input_sha256']==hashes
    done={(r['name'],r['orbit_index'],r['tau']) for r in report['results']};start=time.monotonic()
    for row in source['results']:
        for tau in [4,8]:
            if (row['name'],row['orbit_index'],tau) in done:continue
            assert time.monotonic()-start<120
            result=one(row,panel[row['name']]['pd'],tau)
            report['results'].append(result);output.write_text(json.dumps(report,indent=2)+'\n')
            print(json.dumps({k:result[k] for k in ['name','orbit_index','tau','first_obstruction','seconds']}),flush=True)
            assert result['seconds']<10
    assert len(report['results'])==32
    report['status']='COMPLETE';output.write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
