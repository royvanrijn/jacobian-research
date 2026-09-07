#!/usr/bin/env python3
"""Exact generic trace halves over the retained quadratic branch fields."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import branch_blocks
import branch_divisibility_capacity as branch
import verify_branch_divisibility_capacity as prior

PROTOCOL=Path(__file__).with_name('NATIVE_BRANCH_HALF_SECTION_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_native_branch_half_section_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_native_branch_half_section_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-native-branch-half-section-v1'


def export():
    d=r.read(branch.INPUT);old={x['label']:x for x in r.read(branch_blocks.INPUT)['covers']};rows=[]
    for row in d['covers']:
        c=old[row['label']];assert c['residual_chord']['q_coefficients']==row['q']
        rows.append({**row,'word':c['published_basis_w'],'half_x':c['lifted_section']['x0_coefficients'],
            'half_y':c['lifted_section']['y0_coefficients']})
    r.write_new(INPUT,{'schema':'rank-jump.native-branch-half-section-inputs.v1','A':d['A'],'B':d['B'],
        'sections':d['sections'],'rows':rows,'bindings':branch.bindings([Path(__file__),PROTOCOL,branch.INPUT,branch_blocks.INPUT]),
        'boundary':'Only equation-defined generic native section terms and trace words. Retrospective support selection is unchanged.'})


def compute():
    from sage.all import QQ,PolynomialRing,NumberField,EllipticCurve
    d=r.read(INPUT);previous={x['label']:x for x in r.read(prior.OUTPUT)['rows']};R=PolynomialRing(QQ,'t');rows=[]
    for c in d['rows']:
        q=R(c['q']);K=NumberField(q,'b');b=K.gen();ev=lambda a:K(R(a)(b))
        A,B=ev(d['A']),ev(d['B']);E=EllipticCurve(K,[A,B])
        points=[E(ev(s['x']),ev(s['y'])) for s in d['sections']]
        P=E(ev(c['half_x']),ev(c['half_y']));W=sum((int(n)*p for n,p in zip(c['word'],points)),E(0))
        assert 2*P==W and W!=E(0)
        mask=sum((int(n)%2)<<i for i,n in enumerate(c['word']))
        assert mask==previous[c['label']]['finite_kernel_mask'] and mask
        enc=lambda z:[str(z[i]) for i in range(2)]
        rows.append({'label':c['label'],'kernel_mask':mask,'exact_branch_kernel_dimension':1,
            'half_point':[enc(P[0]),enc(P[1])],'trace_point':[enc(W[0]),enc(W[1])]})
    return {'schema':'rank-jump.native-branch-half-section.v1','status':'PASS','rows':rows,
        'bindings':branch.bindings([Path(__file__),PROTOCOL,INPUT,prior.OUTPUT,Path(r.__file__)]),
        'boundary':'Exact doubling and generic-word identities over quadratic branch fields; no additional global-pool class or twist point beyond the known native section is constructed.'}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        reason=None
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=60)
                if p.returncode:reason='worker failure'
            except subprocess.TimeoutExpired:reason='bounded timeout'
        if reason:r.write_new(OUTPUT,{'status':'UNKNOWN','reason':reason});return
    out=r.read(path);r.write_new(OUTPUT,out);print('PASS',len(out['rows']),'exact branch halves')


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','check']);a=p.parse_args()
    if a.mode=='worker':r.write_new(WORK/'worker.json',compute())
    elif a.mode=='check':assert compute()==r.read(OUTPUT);print('PASS branch halves replay')
    else:globals()[a.mode]()
