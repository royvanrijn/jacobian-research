#!/usr/bin/env sage-python
"""One reduction for each of the seven nonzero words in the retained 3D kernel.

This finite diagnostic accepts a unit only on an exact ideal norm-one witness.
All other outcomes remain UNKNOWN. No BNF, ideal enumeration, or point search.
"""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,vector,pari,prod
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
INPUT=ART/'curve302_even_ideal_extension_complete_v1.json'
BASE=ART/'curve302_relative_ideal_anatomy_v1.json'
EXTRA=ART/'curve302_even_ideal_extension_v1.json'
STRICT=ART/'curve302_descent_anatomy_v1.json'
ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
WORK=ROOT/'artifacts/local/elliptic-curves/curve302-even-artin-kernel-v1'
def read(p):return json.loads(p.read_text())
def save(p,d):
    if p.exists():assert read(p)==d
    else:
        with p.open('x') as out:json.dump(d,out,indent=2,sort_keys=True);out.write('\n')
def rows(M):return [[str(x)for x in row]for row in matrix(QQ,M)]
def compute():
    WORK.mkdir(parents=True,exist_ok=True)
    save(WORK/'protocol.json',{'inputs':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),INPUT,BASE,EXTRA,STRICT,ARITH]},
                             'limits':{'wall_seconds':60,'kernel_words':7,'reductions_per_word':1,'point_searches':0,'bnf_calls':0}})
    d,b,e,s,a=map(read,[INPUT,BASE,EXTRA,STRICT,ARITH]);R=PolynomialRing(QQ,'z');f=R(a['cubic_ascending']);nf=pari.nfinit([pari(f),a['S_finite']])
    gs=[pari.Mod(pari(R([QQ(p['a']),-QQ(p['d'])**2])),pari(f))for p in s['half_ideal_packet']['points']]
    ideals=[pari(matrix(QQ,c['ideal']))for c in b['columns']]+[pari(matrix(QQ,c['ideal']))for c in e['extra_columns']]
    squaregens=[]
    for i in range(8):squaregens.append(pari.Mod(pari(R(read(ART/('det1092_generic_virtual_units_v1/basis-%02d.json'%i))['beta'])),pari(f)))
    for i in range(10):
        c=s['initial_artin_packet']['columns'][i];alpha=pari.Mod(pari(R(c['multiplier'])),pari(f))
        beta=prod(g for bit,g in zip(s['half_ideal_packet']['words'][i],gs)if bit)
        squaregens.append(beta/alpha**2)
    for c in e['extra_columns']:
        alpha=pari.Mod(pari(R(c['multiplier'])),pari(f));beta=prod(g for bit,g in zip(c['word'],gs)if bit)
        squaregens.append(beta/alpha**2)
    assert len(ideals)==len(squaregens)==22
    for J,beta in zip(ideals,squaregens):assert pari.idealpow(nf,J,2)==pari.idealhnf(nf,beta)
    kernel=matrix(GF(2),d['kernel_words']);assert kernel.nrows()==3
    results=[]
    for mask in range(1,8):
        word=vector(GF(2),[(mask>>i)&1 for i in range(3)])*kernel
        J=pari.idealhnf(nf,1);beta=pari.Mod(1,pari(f))
        for bit,I,g in zip(word,ideals,squaregens):
            if bit:J=pari.idealmul(nf,J,I);beta*=g
        H,alpha=pari.idealred(nf,[J,1]);assert pari.idealmul(nf,H,alpha)==J
        alpha=pari.nfbasistoalg(nf,alpha);u=beta/alpha**2;assert pari.idealpow(nf,H,2)==pari.idealhnf(nf,u)
        N=ZZ(pari.idealnorm(nf,H));row={'mask':mask,'word':list(map(int,word)),'norm':str(N),'reduced_ideal':rows(H),
                                     'unit_found':bool(N==1),'square_generator':[str(pari.lift(u).polcoef(i))for i in range(3)]}
        save(WORK/('word-%d.json'%mask),row);results.append(row)
        print('KERNEL',mask,'NORM DIGITS',len(str(N)),'UNIT',N==1,flush=True)
    save(WORK/'summary.json',{'unit_masks':[r['mask']for r in results if r['unit_found']],
                            'boundary':'One reduced representative of each of seven fixed kernel words only. Nonunit ideals are not nonprincipality certificates.'})
if __name__=='__main__':compute()
