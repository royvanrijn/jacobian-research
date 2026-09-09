#!/usr/bin/env sage-python
"""Fixed-bad-prime valuation audit of two existing squareclasses, <=25s."""
import hashlib, json
from pathlib import Path
from sage.all import QQ, ZZ, PolynomialRing, NumberField, pari, matrix

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT/'artifacts/generated-results/elliptic-curves'
OUT = ART/'det1092_seed_half_ideal_v1'
ARITH = ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
CASES = [OUT/(label+'.json') for label in ['302-first-unlock','generic-section-0']]
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
def save(name,data):
    p=OUT/name
    if p.exists(): assert json.loads(p.read_text()) == data
    else:
        with p.open('x') as stream:
            json.dump(data,stream,indent=2,sort_keys=True);stream.write('\n')
def rows(A): return [[str(x) for x in row] for row in A.rows()]
save('ideal-parity-protocol.json', {
    'classification':'retrospective exact ideal-parity diagnostic',
    'inputs':{str(p.relative_to(ROOT)):sha(p) for p in [ARITH,*CASES,Path(__file__)]},
    'limits':{'seconds':25,'classes':2,'rational_primes':20,'class_groups':0,
              'unrestricted_factorizations':0,'point_searches':0},
    'method':'Factor only the already frozen bad rational primes. Correct the factor-free half-ideal there if and only if all valuations are even.'})
d=json.loads(ARITH.read_text());R=PolynomialRing(QQ,'z');f=R(d['cubic_ascending'])
K=NumberField(f,'theta');nf=pari.nfinit([pari(f),d['S_finite']])
assert ZZ(nf.disc())==ZZ(d['field_discriminant']) and ZZ(nf[3])==ZZ(d['defining_order_index'])
primes=[(p,pari.idealprimedec(nf,p)) for p in d['S_finite']]
results=[]
for path in CASES:
    c=json.loads(path.read_text());alpha=K(R(c['alpha']));pa=pari.Mod(pari(R(c['alpha'])),pari(f))
    I=pari(matrix(QQ,c['half_ideal_basis']));J=I;local=[];odd=[]
    for p,above in primes:
        for index,P in enumerate(above):
            v=int(pari.idealval(nf,pa,P));w=int(pari.idealval(nf,I,P))
            local.append({'p':p,'index':index,'e':int(P[2]),'f':int(P[3]),
                          'alpha_valuation':v,'I_valuation':w,
                          'prime_HNF':rows(matrix(QQ,pari.idealhnf(nf,P).sage()))})
            assert 2*w-v >= 0
            if v%2: odd.append([p,index])
            J=pari.idealmul(nf,J,pari.idealpow(nf,P,v//2-w))
    parity_ideal=pari.idealdiv(nf,pa,pari.idealpow(nf,J,2))
    expected=pari.idealhnf(nf,1)
    for row in local:
        if row['alpha_valuation']%2:
            expected=pari.idealmul(nf,expected,pari(matrix(QQ,row['prime_HNF'])))
    assert parity_ideal==expected
    item={'label':c['label'],'local':local,'odd_valuation_primes':odd,
          'corrected_ideal':rows(matrix(QQ,J.sage())),
          'corrected_ideal_norm':str(pari.idealnorm(nf,J)),
          'parity_ideal':rows(matrix(QQ,parity_ideal.sage())),
          'identity':'(alpha)=J^2*product(odd-valuation bad prime ideals)',
          'is_global_square_ideal':not odd,
          'boundary':'A square principal ideal yields a 2-torsion ideal class, possibly trivial. Nontrivial squareclass does not prove nontrivial ideal class.'}
    results.append(item)
    print(c['label'],'odd ideal primes',odd,'square ideal',not odd,flush=True)
save('ideal-parity.json',{'status':'PASS_TWO_EXACT_IDEAL_PARITY_AUDITS','cases':results,
    'classification':'verified application','checker_sha256':sha(Path(__file__))})
