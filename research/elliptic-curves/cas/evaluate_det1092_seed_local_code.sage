#!/usr/bin/env sage-python
"""Evaluate only the historical first seed against the already frozen code.

The constructor's generic output is retained. Its initial evaluation attempt
failed a Sage matrix/vector coercion before writing any seed result.
"""
import hashlib,json,signal,sys
from pathlib import Path
from sage.all import QQ,GF,PolynomialRing,matrix,vector,pari,prod
CAS=Path(__file__).resolve().parent
sys.path.insert(0,str(CAS))
from research_runtime.local_kummer import LocalSquareclasses
ROOT=CAS.parents[1]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
DIR=ART/'det1092_seed_local_code_v3'
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def read(p):return json.loads(p.read_text())
def retain(p,d):
    payload=json.dumps(d,indent=2,sort_keys=True)+'\n'
    if p.exists():assert p.read_text()==payload
    else:p.write_text(payload)
def signs_at_roots(f,x):
    seq=[f,f.derivative()]
    while seq[-1].degree()>0:seq.append(-(seq[-2]%seq[-1]))
    def variations(signs):
        signs=[s for s in signs if s]
        return sum(a!=b for a,b in zip(signs,signs[1:]))
    minus=variations([int(g.leading_coefficient().sign())*(-1)**g.degree() for g in seq])
    plus=variations([int(g.leading_coefficient().sign()) for g in seq])
    at=variations([int(g(x).sign()) for g in seq])
    assert f(x)
    return [0]*(minus-at)+[1]*(at-plus)
def run():
    generic_path=DIR/'case-08-generic.json'
    generic=read(generic_path)
    for path,digest in generic['inputs'].items():assert sha(ROOT/path)==digest
    frozen_hash=sha(generic_path)
    inp=read(DIR/'case-08-input.json')
    R=PolynomialRing(QQ,'x');x=R.gen();f=R(list(map(QQ,inp['polynomial'])))
    pari.allocatemem(64000000,268435456,silent=True)
    nf=pari.nfinit([pari(f),inp['places']]);theta=pari(x).Mod(pari(f))
    assert [list(map(str,R(a).list())) for a in nf.nf_get_zk()]==generic['local_order_basis']
    betas=[pari(QQ(X))-theta for X,Y in inp['generic_cubic_points']]
    # First exceptional point is loaded only after the generic code hash.
    source=ART/'det1092_first_centre_rr_net_replay_v1.json'
    xp,yp=map(QQ,read(source)['reconstructed_point_literal302'])
    X,Y=4*xp,8*yp+4*xp+4
    assert f(X)==Y*Y
    beta=pari(X)-theta;coords=[];tests=[]
    for row in generic['local'][:-1]:
        engine=LocalSquareclasses(nf,row['place'])
        signatures=matrix(GF(2),[list(engine.signature(b)) for b in betas]).transpose().matrix_from_columns(row['basis_indices'])
        cc=list(map(int,signatures.solve_right(vector(GF(2),engine.signature(beta)))))
        quotient=beta/prod(betas[j]**c for j,c in zip(row['basis_indices'],cc))
        vals=[int(pari.idealval(nf,quotient,P)) for P in engine.primes]
        factors=pari.matrix(len(vals),2,[a for P,v in zip(engine.primes,vals) for a in [P,v//2]])
        scale=pari.nfbasistoalg(nf,pari.idealappr(nf,factors))
        assert all(v%2==0 for v in vals)
        tests.append({'place':row['place'],'coordinates':cc,'scale':list(map(str,R(scale.lift()).list())),'expected_square':True})
        coords.extend(cc)
    sr=signs_at_roots(f,X)
    real=generic['local'][-1]
    signatures=matrix(GF(2),real['root_signs']).transpose().matrix_from_columns(real['basis_indices'])
    coords.extend(map(int,signatures.solve_right(vector(GF(2),sr))))
    checks=matrix(GF(2),generic['compatibility_checks']);col=vector(GF(2),coords)
    syndrome=list(map(int,checks*col));assert any(syndrome)
    assert checks.nrows()<=5
    offsets=[];start=0
    for row in generic['local']:
        offsets.append((row['place'],start,start+row['generic_rank']))
        start+=row['generic_rank']
    separators=[]
    for mask in range(1,1<<checks.nrows()):
        word=sum((checks[j] for j in range(checks.nrows()) if mask>>j&1),checks[0]*0)
        if word*col:
            support=[p for p,lo,hi in offsets if any(word[lo:hi])]
            separators.append({'mask':mask,'word':list(map(int,word)),'places':support,'place_count':len(support)})
    separators.sort(key=lambda r:(r['place_count'],r['mask']))
    assert sha(generic_path)==frozen_hash
    result={'classification':'retrospective evaluation of a generic-only frozen code',
        'status':'PASS_FIRST_SEED_LOCAL_COMPATIBILITY_DEFECT',
        'point_cubic':[str(X),str(Y)],'coordinates':coords,'syndrome':syndrome,
        'local_square_tests':tests,'real_signs':sr,'all_separating_checks':separators,
        'minimal_place_count':separators[0]['place_count'],'minimal_separator':separators[0],
        'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [source,generic_path,Path(__file__)]}}
    retain(DIR/'first-seed-evaluation.json',result)
    print('PASS_FIRST_SEED_LOCAL_COMPATIBILITY_DEFECT',syndrome,separators[0]['places'],flush=True)
if __name__=='__main__':
    signal.alarm(25)
    run()
