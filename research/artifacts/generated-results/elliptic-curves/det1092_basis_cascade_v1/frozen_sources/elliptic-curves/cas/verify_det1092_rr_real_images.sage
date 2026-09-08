#!/usr/bin/env sage-python
"""Independent Newton-trace and explicit congruence replay at the real place."""
import hashlib,json,signal
from pathlib import Path
from sage.all import QQ,PolynomialRing,matrix,diagonal_matrix
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUT=ART/'det1092_rr_real_images_v2'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def verify():
    proto=read(OUT/'protocol.json');old=read(ART/'det1092_rr_real_images_v1/protocol.json')
    assert {k:v for k,v in proto.items() if k!='script_sha256'}=={k:v for k,v in old.items() if k!='script_sha256'}
    for p,h in proto['inputs'].items():assert sha(ROOT/p)==h
    original=ROOT/'elliptic-curves/cas/construct_det1092_rr_real_images.sage'
    assert old['script_sha256']==sha(original)
    dependency=ART/'det1092_rr_good_local_images_v1/replay.json'
    for p,h in read(dependency)['inputs'].items():assert sha(ROOT/p)==h
    paths=[OUT/'protocol.json',original,dependency,ART/'det1092_rr_real_images_v1/failure.json']
    R=PolynomialRing(QQ,'T');T=R.gen();results=[]
    for i in range(10):
        path=OUT/('case-%02d.json'%i);d=read(path);paths.append(path)
        for p,h in d['inputs'].items():assert sha(ROOT/p)==h
        source=read(ROOT/d['source']);q=R(source['q']).monic();x0=QQ(source['base_x'])
        assert QQ(source['scale'])*R(source['q'])(x0)!=0
        assert (QQ(source['scale'])*R(source['q'])(x0)).is_square()
        seq=[q,q.derivative()]
        while seq[-1].degree()>0:seq.append(-(seq[-2]%seq[-1]))
        def variations(side):
            signs=[p.leading_coefficient().sign()*side**p.degree() for p in seq]
            return sum(a!=b for a,b in zip(signs,signs[1:]))
        count=variations(-1)-variations(1);assert count==d['real_branch_points']
        dim=max(0,count//2-1);assert dim==d['real_Kummer_dimension']<=1
        # Newton recurrence for traces, independent of the constructor's
        # multiplication matrices and matrix powers.
        a=[q[6-j] for j in range(7)];moments=[QQ(6)]
        for k in range(1,16):
            value=-sum(a[j]*moments[k-j] for j in range(1,min(k,7)))
            if k<=6:value-=k*a[k]
            moments.append(value)
        tests=[];rank=0
        for j,item in enumerate(d['divisor_tests']):
            assert j==item['generic_divisor_index']<17
            g=R(source['generic_divisors'][j]);h=R(item['h_mod_q'])
            assert h==((-1)**g.degree()*g*(x0-T)**g.degree())%q and h.gcd(q)==1
            H=matrix(QQ,6,6,[sum(h[k]*moments[k+u+v] for k in range(6)) for u in range(6) for v in range(6)])
            V=matrix(QQ,[[QQ(a) for a in row] for row in item['congruence']])
            ds=[QQ(a) for a in item['diagonal']];D=diagonal_matrix(QQ,ds)
            assert len(ds)==6 and all(ds) and V.det()
            assert V.transpose()*H*V==D
            signature=sum(a.sign() for a in ds);mixed=abs(signature)<count
            assert signature==item['signature'] and mixed==item['mixed_real_signs']
            tests.append({'generic_divisor_index':j,'signature':int(signature),'mixed_signs':mixed})
            if mixed:
                rank=1;assert j==len(d['divisor_tests'])-1
        assert rank==dim==d['generic_image_rank']
        result={'classification':'verified application and independent replay',
                'case_index':i,'real_branch_points':count,'complete_real_Kummer_dimension':dim,
                'generic_image_rank':rank,'divisor_tests':tests,
                'status':'PASS_INDEPENDENT_COMPLETE_REAL_IMAGE',
                'checker_sha256':sha(Path(__file__))}
        retain(OUT/('case-%02d-replay.json'%i),result);results.append(result)
    return {'classification':'verified application and new deduction',
            'status':'PASS_INDEPENDENT_TEN_COMPLETE_REAL_KUMMER_IMAGES',
            'cases':results,'limits':{**proto['limits'],'independent_panel_wall_seconds':25},
            'scope':'All ten real local images are spanned by inherited generic divisors. No global Selmer or point search.',
            'inputs':{str(p.relative_to(ROOT)):sha(p) for p in paths},'checker_sha256':sha(Path(__file__))}
if __name__=='__main__':
    signal.alarm(25);result=verify();retain(OUT/'replay.json',result)
    print(result['status'],flush=True)
