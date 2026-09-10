#!/usr/bin/env sage-python
"""Test the two fixed historical seed norm ideals with the20 known characters.

No factorization, reductions or new characters. Nonunit residues remain UNKNOWN.
The replay uses local Hilbert symbols and an independent integer Jacobi routine.
"""
import argparse,hashlib,json,runpy,sys
from pathlib import Path
from sage.all import QQ,ZZ,NumberField,PolynomialRing,matrix,vector,pari,prod
ROOT=Path(__file__).resolve().parents[2];ART=ROOT/'artifacts/generated-results/elliptic-curves'
OUTPUT=ART/'curve302_seed_norm_artin_obstruction_v1.json'
sys.path.insert(0,str(ROOT/'elliptic-curves/rank-jump'))
from verify_half_ideal_artin import jacobi
def read(p):return json.loads(p.read_text())
def compute(replay=False):
    paths=[ART/'curve302_relative_ideal_anatomy_v1.json',ART/'curve302_descent_anatomy_v1.json',
           ART/'rank_jump_curve302_strict_constructor_arithmetic_v1.json']
    paths += [ART/('det1092_seed_artin_v2/norm-form-%02d.json'%i)for i in range(2)]
    d,s,a=map(read,paths[:3]);R=PolynomialRing(QQ,'z');f=R(a['cubic_ascending'])
    K=NumberField(f,'theta');theta=K.gen();nf=pari.nfinit([pari(f),a['S_finite']])
    basis=[K(R(b))for b in a['maximal_order_basis']]
    B=matrix(QQ,[list(b)for b in basis]).transpose();Bi=B.inverse()
    gs=[ZZ(p['a'])-ZZ(p['d'])**2*theta for p in s['half_ideal_packet']['points']]
    betas=[prod(g for bit,g in zip(w,gs)if bit)for w in d['characters']['words']]
    def pa(g):return pari.Mod(pari(R(list(g))),pari(f))
    places=[(p,j,P)for p in a['S_finite']for j,P in enumerate(pari.idealprimedec(nf,p))]
    rows=[]
    for index,path in enumerate(paths[3:]):
        source=read(path);H=matrix(QQ,source['ideal_basis']);assert H.det()==ZZ(source['ideal_norm'])
        good=pari(H);bad=[]
        for p,j,P in places:
            e=int(pari.idealval(nf,pari(H),P));assert e>=0
            if e:good=pari.idealmul(nf,good,pari.idealpow(nf,P,-e));bad.append((p,j,P,e))
        J=matrix(QQ,good);N=ZZ(J.det());assert N>0 and all(N%p for p in a['S_finite'])
        assert J[0,0]==N and J[1,1]==J[2,2]==1
        residues=vector(QQ,[1,-J[0,1],-J[0,2]])
        assert all(J[i,j]==0 for i in range(3)for j in range(i))
        for i in range(3):
            for j in range(3):
                assert (residues*Bi*vector(QQ,list(basis[i]*basis[j]))-residues[i]*residues[j])%N==0
        entries=[]
        for k,beta in enumerate(betas):
            value=ZZ(residues*Bi*vector(QQ,list(beta)))%N
            if value.gcd(N)!=1:
                entries.append({'character':k,'status':'UNKNOWN_NONUNIT_RESIDUE'});continue
            symbol=jacobi(int(value),int(N))if replay else int(pari.kronecker(value,N))
            assert symbol in [-1,1];bits=[]
            for p,j,P,e in bad:
                if replay:
                    pi=pari.nfbasistoalg(nf,pari.idealappr(nf,P));assert pari.idealval(nf,pi,P)==1
                    bit=int(pari.nfhilbert(nf,pa(beta),pi,P)==-1)
                else:bit=int(not pari.nfislocalpower(nf,P,pa(beta),2))
                bits.append({'p':p,'j':j,'exponent':e,'bit':bit})
            bit=(int(symbol==-1)+sum(t['exponent']*t['bit']for t in bits))%2
            entries.append({'character':k,'status':'PASS','residue':str(value),'symbol':symbol,'bad_terms':bits,'bit':bit})
        witnesses=[e['character']for e in entries if e.get('bit')==1]
        rows.append({'index':index,'good_ideal':[[str(x)for x in row]for row in J],'good_norm':str(N),
                     'entries':entries,'nonprincipal_witness_characters':witnesses,
                     'integral_norm_equation':'NO_INTEGER_SOLUTION'if witnesses else 'UNKNOWN'})
        print('NORM EQUATION',index,rows[-1]['integral_norm_equation'],'WITNESSES',witnesses,flush=True)
    return {'schema':'curve302.seed-norm-artin-obstruction.v1','cases':rows,
            'bindings':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest()for p in [Path(__file__),*paths]},
            'limits':{'wall_seconds':60,'ideals':2,'characters':20,'new_factorization':0,'bnf_calls':0},
            'boundary':'Only nonzero certified ordinary Artin characters certify nonprincipality. These characters use the known31-point group, not generic sections alone.'}
if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    result=compute(args.check)
    if args.check:assert result==read(OUTPUT)
    else:
        with OUTPUT.open('x')as out:json.dump(result,out,indent=2,sort_keys=True);out.write('\n')
