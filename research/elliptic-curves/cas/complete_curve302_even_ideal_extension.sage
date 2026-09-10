#!/usr/bin/env sage-python
"""Repair the sole fixed noncyclic quotient by separating its 47-primary part."""
import hashlib,json
from pathlib import Path
from sage.all import QQ,ZZ,GF,PolynomialRing,matrix,pari,prod
ROOT=Path(__file__).resolve().parents[2]
ART=ROOT/'artifacts/generated-results/elliptic-curves'
BASE=ART/'curve302_relative_ideal_anatomy_v1.json'
EXTRA=ART/'curve302_even_ideal_extension_v1.json'
STRICT=ART/'curve302_descent_anatomy_v1.json'
ARITH=ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json'
OUT=ART/'curve302_even_ideal_extension_complete_v1.json'
def read(p):return json.loads(p.read_text())
def rows(M):return [[str(x) for x in row] for row in matrix(QQ,M).rows()]
def compute():
    b,d,s,a=map(read,[BASE,EXTRA,STRICT,ARITH])
    c=d['extra_columns'][1]
    assert all(e['status']=='UNKNOWN_NONCYCLIC' for e in c['entries'])
    R=PolynomialRing(QQ,'z');f=R(a['cubic_ascending']);nf=pari.nfinit([pari(f),a['S_finite']])
    H=pari(matrix(QQ,c['ideal']));good=H;terms=[]
    for p in a['S_finite']+[47]:
        for j,P in enumerate(pari.idealprimedec(nf,p)):
            e=int(pari.idealval(nf,H,P))
            if e:
                assert e>0;terms.append((p,j,P,e));good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e))
    N=ZZ(pari.idealnorm(nf,good));assert N>0 and all(N%p for p in a['S_finite']+[47])
    assert good[0,0]==N and good[1,1]==good[2,2]==1
    points=s['half_ideal_packet']['points'];gs=[pari.Mod(pari(R([QQ(p['a']),-QQ(p['d'])**2])),pari(f)) for p in points]
    entries=[]
    for word in b['characters']['words']:
        beta=prod(g for bit,g in zip(word,gs) if bit);coords=pari.nfalgtobasis(nf,beta)
        res=ZZ(coords[0]-good[0,1]*coords[1]-good[0,2]*coords[2])%N
        assert res.gcd(N)==1
        symbol=int(pari.kronecker(res,N));assert symbol in [-1,1]
        local=[{'p':p,'j':j,'exponent':e,'degree':int(P[3]),'valuation':int(pari.idealval(nf,beta,P)),
                'bit':int(not pari.nfislocalpower(nf,P,beta,2))} for p,j,P,e in terms]
        bit=(int(symbol==-1)+sum(t['bit']*t['exponent'] for t in local))%2
        entries.append({'residue':str(res),'symbol':symbol,'local_terms':local,'bit':bit})
    extra_cols=[]
    for i,c0 in enumerate(d['extra_columns']):extra_cols.append([e['bit'] for e in (entries if i==1 else c0['entries'])])
    M=matrix(GF(2),b['artin_matrix']).augment(matrix(GF(2),extra_cols).transpose())
    T=matrix(GF(2),d['basis_words']).transpose()
    U=T.solve_right(matrix(GF(2),b['characters']['words']).transpose())
    rank=int(M.rank())
    result={'schema':'curve302.even-ideal-extension.complete.v1',
            'bindings':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__),BASE,EXTRA,STRICT,ARITH]},
            'repair':{'column':1,'ideal':rows(H),'good_ideal':rows(good),'norm':str(N),'entries':entries},
            'artin_matrix':[list(map(int,row)) for row in M], 'detected_rank':rank,
            'kernel_words':[list(map(int,row)) for row in M.right_kernel().basis()],
            'unramified_in_even_basis':[list(map(int,row)) for row in U.transpose()],
            'unramified_half_ideal_detected_rank':int((M*U).rank()),
            'full_even_ideal_rank_interval':[rank,22], 'full_even_unit_kernel_interval':[0,22-rank],
            'relative_ideal_rank_interval':[rank-8,14],
            'boundary':'Character rank gives lower bounds. Kernel ideals are not proved principal or units.'}
    if OUT.exists():assert read(OUT)==result
    else:
        with OUT.open('x') as stream:json.dump(result,stream,indent=2,sort_keys=True);stream.write('\n')
    print('PASS FULL EVEN',rank,'UNRAMIFIED',result['unramified_half_ideal_detected_rank'],
          'RELATIVE',result['relative_ideal_rank_interval'],'KERNEL',result['kernel_words'],flush=True)
if __name__=='__main__':compute()
