#!/usr/bin/env python3
"""Exact maximal-order norm-form transfer from the reference to302."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r

PROTOCOL=Path(__file__).with_name('CURVE302_MAXIMAL_NORM_FORM_PROTOCOL.json')
SOURCES=[r.OUT/'rank_jump_reference_strict_class_construction_inputs_v1.json',r.OUT/'rank_jump_curve302_strict_constructor_arithmetic_v1.json']
CALIBRATION=r.OUT/'rank_jump_reference_class_targeted_relations_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_curve302_maximal_norm_form_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-curve302-maximal-norm-form-v1'

def compute():
    from sage.all import QQ,ZZ,PolynomialRing,pari
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    import prepare_small_conductor_norm_form as forms
    pari.allocatemem(64000000,r.read(PROTOCOL)['bounds']['pari_stack_bytes'],silent=True)
    cases=[]
    for source in SOURCES:
        data=r.read(source);R=PolynomialRing(QQ,'z');f=R(data['cubic_ascending']);nf=pari.nfinit([pari(f),data['S_finite']])
        assert str(nf.disc())==data['field_discriminant']
        zk=list(nf.nf_get_zk());assert zk[0]==1
        pair=list(map(ZZ,pari.nfalgtobasis(nf,zk[1]*zk[2])));w=zk[1]-pair[2];t=zk[2]-pair[1]
        def coords(value):
            v=list(map(ZZ,pari.nfalgtobasis(nf,value)))
            return [v[0]+v[1]*pair[2]+v[2]*pair[1],v[1],v[2]]
        table=[[coords(u*v) for v in [1,w,t]] for u in [1,w,t]]
        a,b,c,d=-table[1][1][2],table[1][1][1],-table[2][2][2],table[2][2][1]
        assert table[1][1]==[-a*c,b,-a] and table[1][2]==[-a*d,0,0] and table[2][2]==[-b*d,d,-c]
        initial=list(map(int,[a,b,c,d]));reduced,M,history=forms.reduce_form(initial)
        assert forms.discriminant(initial)==nf.disc()==forms.discriminant(reduced)
        assert M[0]*M[3]-M[1]*M[2]==1 and forms.transform(initial,M)==reduced
        for m,n in [(1,0),(0,1),(1,1),(-1,1),(2,1)]:
            alpha=a*(M[0]*m+M[1]*n)+(M[2]*m+M[3]*n)*w
            expected=a*a*sum(ZZ(coef)*m**(3-i)*n**i for i,coef in enumerate(reduced))
            assert pari.nfeltnorm(nf,alpha)==expected
        result={'source':str(source.relative_to(r.ROOT)),'cubic_ascending':data['cubic_ascending'],
            'S_finite':data['S_finite'],'field_discriminant':str(nf.disc()),'defining_order_index':str(nf[3]),
            'maximal_order_basis':list(map(str,zk)),'normal_basis_GP':list(map(str,[1,w,t])),
            'multiplication_table':[[list(map(str,row)) for row in block] for block in table],
            'initial_binary_cubic_descending':list(map(str,initial)),'binary_cubic_descending':list(map(str,reduced)),
            'sl2_matrix':[[M[0],M[1]],[M[2],M[3]]],'reduction_steps':history,
            'fixed_a':str(a),'w_power_basis':[str(pari.lift(w).polcoef(i)) for i in range(3)],
            'reduced_hessian':list(map(str,forms.hessian(reduced))),'slope_scale':str(2**96),
            'maximum_reduced_coefficient_bits':max(abs(v).bit_length() for v in reduced),
            'field_discriminant_bits':abs(int(nf.disc())).bit_length(),
            'norm_identity':'Norm(a*(M00*m+M01*n)+(M10*m+M11*n)*w)=a^2*F(m,n)'}
        if not cases:
            cal=r.read(CALIBRATION)
            for key in ['cubic_ascending','field_discriminant','S_finite','fixed_a','w_power_basis','sl2_matrix','binary_cubic_descending','slope_scale']:
                assert result[key]==cal[key],key
        cases.append(result);r.write_new(WORK/('case_%d.json'%len(cases)),result)
        print('CASE',len(cases),'FIELD_BITS',result['field_discriminant_bits'],'FORM_BITS',result['maximum_reduced_coefficient_bits'],flush=True)
    return {'schema':'rank-jump.curve302-maximal-norm-form.v1','status':'PASS','cases':cases,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,*SOURCES,CALIBRATION,Path(forms.__file__)]},
        'boundary':'Exact constructor coordinates and reference replay only. No new302 principal relation, strict class, rational point or rank feature.'}

def capture():
    WORK.mkdir(parents=True,exist_ok=True)
    with (WORK/'worker.log').open('x') as log:
        try:
            p=subprocess.run([sys.executable,__file__,'worker'],stdout=log,stderr=log,timeout=r.read(PROTOCOL)['bounds']['worker_seconds'])
            error='worker failure' if p.returncode else None
        except subprocess.TimeoutExpired:error='bounded timeout'
    if error and not OUTPUT.exists():r.write_new(OUTPUT,{'status':'UNKNOWN','reason':error})
    print(r.read(OUTPUT)['status'],flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);a=p.parse_args()
    if a.mode=='worker':r.write_new(OUTPUT,compute())
    else:capture()
