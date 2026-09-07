#!/usr/bin/env python3
"""A fixed relation-root constructor, tested by exact good-prime valuations."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import early_relation_pool as pool
import retained_norm_inherited_hit as hit
import retained_norm_batch_capacity as batch
import bounded_gain_reference as ref

PROTOCOL=Path(__file__).with_name('RELATION_ROOT_CLASS_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_relation_root_class_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_relation_root_class_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-relation-root-class-v1'


def bindings(paths):
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def export():
    from sage.all import QQ
    d=r.read(pool.INPUT);h=r.read(hit.OUTPUT)
    row=next(x for x in d['relations'] if [x['m'],x['n']]==h['address'])
    primes=sorted({d['columns'][i]['p'] for i,e in row['ideal_factorization']}-set(d['S_finite']))[:32]
    ref.configure();f,pts,scale=ref.base.model_data(ref.TOKEN)
    b,a=map(QQ,d['affine_map_old_root_from_masked'])
    assert a.is_square()
    generic=[[str(a*x+b),'-1','0'] for x,y in pts]
    result={'schema':'rank-jump.relation-root-class-inputs.v1',
        'cubic_ascending':d['cubic_ascending'],'field_discriminant':d['field_discriminant'],
        'S_finite':d['S_finite'],'probe_primes':primes,
        'root_ascending':h['square_root_ascending'],'alpha_ascending':h['alpha_ascending'],
        'generic_classes_ascending':generic,'generic_product_mask':h['generic_product_mask'],
        'parent_address':h['address'],
        'bindings':bindings([Path(__file__),PROTOCOL,pool.INPUT,hit.OUTPUT,ref.INPUT,
            Path(ref.__file__),Path(ref.base.__file__),Path(r.__file__)])}
    r.write_new(INPUT,result)
    print('Exported one candidate and',len(primes),'fixed probe primes',flush=True)


def compute():
    from sage.all import QQ,ZZ,pari,PolynomialRing,lcm
    from sage.env import SAGE_VERSION
    pari.allocatemem(64000000,268435456,silent=True)
    d=r.read(INPUT);R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending'])
    nf=pari.nfinit([pari(f),d['S_finite']]);assert str(nf.disc())==d['field_discriminant']
    th=pari.Mod('z',pari(f));elt=lambda coeff:pari(R(coeff))(th)
    w=elt(d['root_ascending']);alpha=elt(d['alpha_ascending'])
    gamma=list(map(elt,d['generic_classes_ascending']))
    parent=pari.nfeltnorm(nf,alpha)*alpha
    for i,g in enumerate(gamma):
        assert QQ(str(pari.nfeltnorm(nf,g))).is_square()
        if d['generic_product_mask']>>i&1:parent*=g
    assert w*w==parent
    N=pari.nfeltnorm(nf,w);beta=N*w
    assert pari.nfeltnorm(nf,beta)==N**4
    enc=lambda z:[str(pari.lift(z).polcoef(i)) for i in range(3)]
    basis=[enc(z) for z in nf.nf_get_zk()]
    rows=[]
    for p in d['probe_primes']:
        assert p not in d['S_finite'] and ZZ(p).is_prime(proof=True)
        assert ZZ(f.discriminant()).valuation(p)==0
        for P in pari.idealprimedec(nf,p):
            vr=int(pari.idealval(nf,w,P));vb=int(pari.idealval(nf,beta,P))
            gv=[int(pari.idealval(nf,g,P)) for g in gamma]
            assert all(v%2==0 for v in gv)
            assert vb==vr+int(P[2])*QQ(str(N)).valuation(p)
            def membership_data(z,v):
                coords=list(pari.nfalgtobasis(nf,z));den=ZZ(lcm([QQ(str(c)).denominator() for c in coords]))
                k=v+int(P[2])*den.valuation(p);assert k>=0
                hnf=lambda j:[[str(pari.idealpow(nf,P,j)[a,b]) for b in range(3)] for a in range(3)]
                return {'denominator':str(den),'integral_coordinates':[str(den*QQ(str(c))) for c in coords],
                    'integral_valuation':k,'ideal_power_hnf':hnf(k),'next_ideal_power_hnf':hnf(k+1)}
            rows.append({'p':p,'e':int(P[2]),'f':int(P[3]),'prime_hnf':str(pari.idealhnf(nf,P)),
                'root_valuation':vr,'projection_valuation':vb,'generic_valuations':gv,
                'root_membership':membership_data(w,vr),'projection_membership':membership_data(beta,vb)})
    witnesses=[i for i,row in enumerate(rows) if row['projection_valuation']%2]
    return {'schema':'rank-jump.relation-root-class.v1','status':'PASS',
        'candidate_status':'EXCLUDED_AT_GOOD_PRIME' if witnesses else 'UNKNOWN',
        'norm_root':str(N),'projection_ascending':enc(beta),'maximal_order_basis':basis,
        'rows':rows,'odd_valuation_witness_rows':witnesses,
        'generic_corrections_excluded':2**len(gamma) if witnesses else 0,
        'boundary':'Any witness excludes this represented class and every product with G from Selmer. No witness is not a membership certificate. The root of an inherited relation need not remain in its parent squareclass span.',
        'software':{'sage':SAGE_VERSION,'pari':str(pari('version()'))},
        'bindings':bindings([Path(__file__),PROTOCOL,INPUT,Path(r.__file__)])}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);path=WORK/'worker.json'
    if not path.exists():
        reason=None
        with (WORK/'worker.log').open('x') as log:
            try:
                p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker'],stdout=log,stderr=log,timeout=30)
                if p.returncode:reason='Worker failure; see retained checkpoint'
            except subprocess.TimeoutExpired:reason='Bounded worker timeout'
        if reason:r.write_new(path,{'status':'UNKNOWN','reason':reason})
    result=r.read(path);r.write_new(OUTPUT,result)
    print(result.get('candidate_status',result['status']),len(result.get('odd_valuation_witness_rows',[])),'witnesses',flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','check']);a=p.parse_args()
    if a.mode=='worker':r.write_new(WORK/'worker.json',compute())
    elif a.mode=='check':assert compute()==r.read(OUTPUT);print('PASS relation-root replay')
    else:globals()[a.mode]()
