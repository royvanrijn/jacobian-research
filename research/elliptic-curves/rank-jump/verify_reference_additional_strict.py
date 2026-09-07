#!/usr/bin/env python3
"""Certify explicit factored strict classes by independent Hilbert tests."""
import argparse
from math import prod
from pathlib import Path
import subprocess
import sys
import retrospective as r

PROTOCOL=Path(__file__).with_name('REFERENCE_ADDITIONAL_STRICT_VERIFICATION_PROTOCOL.json')
SOURCE=r.OUT/'rank_jump_reference_large_support_extraction_v1.json'
REFERENCE=r.OUT/'rank_jump_reference_strict_class_construction_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_reference_additional_strict_verification_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-reference-additional-strict-verification-v1'

def compute():
    from sage.all import QQ,ZZ,AA,GF,PolynomialRing,pari,matrix
    d=r.read(SOURCE);ref=r.read(REFERENCE);policy=r.read(PROTOCOL)
    for path,sha in d['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    pari.allocatemem(64000000,policy['bounds']['pari_stack_bytes'],silent=True)
    R=PolynomialRing(QQ,'z');f=R(ref['cubic_ascending']);nf=pari.nfinit([pari(f),ref['S_finite']])
    roots=f.roots(AA,multiplicities=False);S=set(ref['S_finite']);columns=d['columns']
    gammas=[pari.Mod(pari(R(g['beta_ascending'])),pari(f)) for g in ref['generic_classes']]
    chosen=d['additional_independence_witnesses'];assert len(chosen)==2
    records=[d['strict_candidates'][j] for j in chosen]
    used_ids={label['index'] for c in records for label in c['factor_labels'] if label['kind']=='projected_atom'}
    atoms={x['index']:x for x in d['used_atoms'] if x['index'] in used_ids}
    prime_cache={};factors={('generic',i):g for i,g in enumerate(gammas)};idealrows={}
    valuation_checks=0
    for index,atom in atoms.items():
        alpha=pari.Mod(pari(R(atom['alpha_ascending'])),pari(f));N=ZZ(pari.nfeltnorm(nf,alpha));assert str(N)==atom['norm']
        valuations=dict(atom['ideal_factorization']);assert all(e>=0 for e in valuations.values())
        assert prod(ZZ(columns[i]['p'])**(columns[i]['f']*e) for i,e in valuations.items())==abs(N)
        normvals={}
        for i,e in valuations.items():
            c=columns[i];p=c['p'];normvals[p]=normvals.get(p,0)+c['f']*e
            if p not in prime_cache:prime_cache[p]={str(pari.idealhnf(nf,P)):P for P in pari.idealprimedec(nf,p)}
            P=prime_cache[p][c['hnf']];assert int(pari.idealval(nf,alpha,P))==e;valuation_checks+=1
        projected={}
        for i,c in enumerate(columns):
            p=c['p']
            if p in normvals:projected[i]=valuations.get(i,0)+c['e']*normvals[p]
        idealrows[index]=projected;factors[('projected_atom',index)]=N*alpha
    print('ATOMS',len(atoms),'EXACT_VALUATIONS',valuation_checks,flush=True)
    local=[]
    for p in sorted(S):
        for P in pari.idealprimedec(nf,p):
            pi=pari.nfbasistoalg(nf,pari.idealappr(nf,P));assert int(pari.idealval(nf,pi,P))==1
            gens=[pi]
            if p==2:
                assert int(P[3])==1
                gens += [1+pi**k for k in range(1,2*int(P[2])+1)]
                expected=int(P[2])*int(P[3])+2
            else:
                mod=pari.nfmodprinit(nf,P);unit=None
                for b in nf.nf_get_zk():
                    for k in range(17):
                        u=b+k
                        if pari.idealval(nf,u,P)==0 and not pari.issquare(pari.nfmodpr(nf,u,mod)):
                            unit=u;break
                    if unit is not None:break
                assert unit is not None;gens.append(unit);expected=2
            pairing=[[int(pari.nfhilbert(nf,a,b,P)==-1) for b in gens] for a in gens]
            assert matrix(GF(2),pairing).rank()==expected
            signatures={key:[int(pari.nfhilbert(nf,a,b,P)==-1) for b in gens] for key,a in factors.items()}
            for record in records:
                total=[0]*len(gens)
                for label in record['factor_labels']:
                    bits=signatures[(label['kind'],label['index'])];total=[a^b for a,b in zip(total,bits)]
                assert not any(total)
            local.append({'p':p,'hnf':str(pari.idealhnf(nf,P)),'generators_GP':list(map(str,gens)),
                          'hilbert_gram':pairing,'squareclass_dimension':expected,'strict_pairing_rows':[[0]*len(gens) for record in records]})
        print('LOCAL',p,'PASS',flush=True)
    polys={key:R([QQ(pari.lift(a).polcoef(i)) for i in range(3)]) for key,a in factors.items()}
    class_records=[]
    for j,record in zip(chosen,records):
        exponents={};signs=[0]*3
        for label in record['factor_labels']:
            key=(label['kind'],label['index']);b=polys[key]
            signs=[a^int(b(x)<0) for a,x in zip(signs,roots)]
            if label['kind']=='projected_atom':
                for i,e in idealrows[label['index']].items():exponents[i]=exponents.get(i,0)+e
        assert not any(signs)
        # Every generic beta has square norm and primitive linear point form,
        # so its outside-S valuations are even independently of factorization.
        assert all(e%2==0 for i,e in exponents.items() if columns[i]['p'] not in S)
        class_records.append({'candidate_index':j,'factor_labels':record['factor_labels'],
            'projected_atom_square_ideal_outside_S':[[i,e//2] for i,e in sorted(exponents.items()) if columns[i]['p'] not in S],
            'generic_correction_indices':[x['index'] for x in record['factor_labels'] if x['kind']=='generic'],
            'real_sign_parities':signs})
    for g,b in zip(ref['generic_classes'],gammas):
        pol=R(g['beta_ascending']);assert pol.degree()==1 and pol[1]<0 and ZZ(-pol[1]).is_square()
        assert ZZ(pol[0]).gcd(ZZ(-pol[1]).sqrt())==1
        assert ZZ(pari.nfeltnorm(nf,b))==ZZ(g['norm']) and ZZ(g['norm']).is_square()
    fullchars=list(map(int,d['factor_proof_characters'][:16]))+[int(c['proof_character']) for c in records]
    assert r.rank(fullchars)==18
    coordinates=[];selected=[]
    for block,p in enumerate(d['proof_primes']):
        rr=f.change_ring(GF(p)).roots(multiplicities=False);assert len(rr)==3
        for k,x in enumerate(rr):
            coordinate=3*block+k;v=r.pack([(c>>coordinate)&1 for c in fullchars])
            if r.rank(selected+[v])==len(selected):continue
            selected.append(v);coordinates.append((p,int(x),coordinate))
            if len(selected)==18:break
        if len(selected)==18:break
    assert len(selected)==18
    rows=[[] for _ in range(18)]
    for p,root,coordinate in coordinates:
        bits={}
        for key,b in polys.items():
            value=int(b.change_ring(GF(p))(root));assert value
            bits[key]=int(pow(value,(p-1)//2,p)==p-1)
        col=[bits[('generic',j)] for j in range(16)]
        for record in records:col.append(sum(bits[(x['kind'],x['index'])] for x in record['factor_labels'])%2)
        assert col==[(c>>coordinate)&1 for c in fullchars]
        for row,bit in zip(rows,col):row.append(bit)
    assert matrix(GF(2),rows).rank()==18
    return {'schema':'rank-jump.reference-additional-strict-verification.v1','status':'PASS',
        'additional_independent_strict_classes':2,'generic_mod2_rank':16,'combined_mod2_rank':18,
        'equation':d['equation'],'class_representatives':class_records,
        'representation':'For each factor labelled projected_atom i use Norm(alpha_i)*alpha_i; for each generic i use beta_i from the frozen generic projection. Multiply the listed factors in Q[z]/(f). This is an explicit factored squareclass representative.',
        'atoms':[atoms[i] for i in sorted(atoms)],'columns':columns,'local_hilbert_certificates':local,
        'independence_coordinates':[{'p':p,'root':root} for p,root,c in coordinates],'independence_matrix':rows,
        'exact_principal_atom_valuation_checks':valuation_checks,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in [Path(__file__),PROTOCOL,SOURCE,REFERENCE]},
        'rational_or_Sha':'UNKNOWN','specialization_event':'An explicit fibre-specific principal-ideal parity dependency has been constructed. Its variation in t and the geometric/arithmetic event creating it have not yet been identified.',
        'boundary':'Unconditional individually verified strict Selmer classes outside G. No exceptional points were construction inputs. No full Selmer dimension, rank increase, rational solubility, or class-group completeness is asserted.'}

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
