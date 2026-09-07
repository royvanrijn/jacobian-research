#!/usr/bin/env python3
"""Exact ordered-root topology on the six retained cubic root curves."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import root_curve_frobenius_capacity as source
import generic_selmer_capacity as prior

PROTOCOL=Path(__file__).with_name('ROOT_CURVE_REAL_COMPONENTS_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_root_curve_real_components_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-root-curve-real-components-v1'


def graph(events,arities):
    vertices=[(i,j) for i,a in enumerate(arities) for j in range(a)];edges=[];n=len(events)
    for i,e in enumerate(events):
        left=(i-1)%n;right=i;a,b=arities[left],arities[right]
        pair=(1,2) if e['B_sign']>0 else (0,1);persistent=0 if e['B_sign']>0 else 2
        if e['type']=='I1':
            assert sorted([a,b])==[1,3];side=left if a==3 else right;single=right if a==3 else left
            edges.extend([[(side,pair[0]),(side,pair[1])],[(side,persistent),(single,0)]])
        elif e['node_tangent_discriminant_sign']>0:
            assert a==b==3
            permutation=list(range(3));permutation[pair[0]],permutation[pair[1]]=pair[1],pair[0]
            edges.extend([[(left,k),(right,permutation[k])] for k in range(3)])
        else:
            assert a==b==1;edges.append([(left,0),(right,0)])
    parent={v:v for v in vertices}
    def root(v):
        while parent[v]!=v:v=parent[v]
        return v
    for a,b in edges:parent[root(a)]=root(b)
    groups={}
    for v in vertices:groups.setdefault(root(v),[]).append(v)
    return edges,sorted(groups.values())


def worker(family):
    from sage.all import QQ,PolynomialRing,prod
    f=next(x for x in r.read(source.INPUT)['families'] if x['family']==family)
    R=PolynomialRing(QQ,'t');A=R(f['A']);B=R(f['B']);D=R(f['discriminant'])
    H=prod(R(q['coefficients_ascending']) for q in f['squarefree_factors']);assert H.gcd(B)==1
    J=H*B.squarefree_part();intervals=J.real_root_intervals();events=[]
    doubles=[R(q['coefficients_ascending']) for q in f['squarefree_factors'] if q['multiplicity']==2]
    for (a,b),multiplicity in intervals:
        assert multiplicity==1 and a<b and H(a)*H(b)!=0 and B(a)*B(b)!=0
        if H(a)*H(b)>0:continue
        assert B(a)*B(b)>0
        event={'interval':[str(a),str(b)],'B_sign':int(B(a).sign()),'type':'I1'}
        for q in doubles:
            if q(a)*q(b)>=0:continue
            assert q.degree()==1;t0=-q[0]/q[1];x0=-3*B(t0)/(2*A(t0))
            delta=A.derivative()(t0)**2-6*x0*(A.derivative(2)(t0)*x0+B.derivative(2)(t0))
            assert delta and A.derivative()(t0)*x0+B.derivative()(t0)==0
            event.update({'type':'I2','node_parameter':str(t0),'node_double_root':str(x0),
                'node_tangent_discriminant':str(delta),'node_tangent_discriminant_sign':int(delta.sign())})
        events.append(event)
    samples=[];arities=[]
    for i,e in enumerate(events):
        sample=(QQ(e['interval'][1])+QQ(events[i+1]['interval'][0]))/2 if i+1<len(events) else QQ(e['interval'][1])+1
        assert D(sample)!=0;samples.append(str(sample));arities.append(3 if D(sample)>0 else 1)
    if events:edges,groups=graph(events,arities);count=len(groups)
    else:edges=[];groups=[];count=3 if D.leading_coefficient()>0 else 1
    g=f['root_curve_genus'];assert 1<=count<=g+1
    result={'family':family,'status':'PASS','genus':g,'events':events,'interval_samples':samples,
        'interval_real_root_counts':arities,'edges':edges,'components':groups,'real_components':count,
        'real_jacobian_two_torsion_dimension':g+count-1,
        'bindings':prior.bind([Path(__file__),PROTOCOL,source.INPUT,Path(r.__file__)])}
    r.write_new(WORK/(family+'.json'),result)


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for f in r.read(source.INPUT)['families']:
        family=f['family'];path=WORK/(family+'.json')
        if not path.exists():
            error=None
            with (WORK/(family+'.log')).open('x') as log:
                try:
                    p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--family',family],stdout=log,stderr=log,timeout=30)
                    if p.returncode:error='worker failure'
                except subprocess.TimeoutExpired:error='bounded timeout'
            if error:r.write_new(path,{'family':family,'status':'UNKNOWN','reason':error})
        row=r.read(path);rows.append(row);print(family,row['status'],row.get('real_components'),row.get('real_jacobian_two_torsion_dimension'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.root-curve-real-components.v1','rows':rows,
        'bindings':prior.bind([Path(__file__),PROTOCOL,source.INPUT]),
        'boundary':'Exact real topology of normalized root curves; arithmetic rational 2-torsion may be smaller. No exceptional points or rank labels entered workers.'})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--family');args=p.parse_args()
    if args.mode=='worker':worker(args.family)
    else:capture()
