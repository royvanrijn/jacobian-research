#!/usr/bin/env python3
"""Full generic half-ideal images; factorization-free regularized Artin evaluation."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import strict_half_ideals as half
import remaining_bad_primes as rem
import half_ideal_artin_completion as old_artin
import strict_class_blocks as blocks

PROTOCOL=Path(__file__).with_name('GENERIC_SUNIT_CARRIER_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_generic_sunit_carrier_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-generic-sunit-carrier-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
        (Path(__file__),PROTOCOL,half.INPUT,rem.INPUT,old_artin.OUTPUT,blocks.OUTPUT)}


def worker(index):
    from sage.all import QQ,ZZ,PolynomialRing,pari,matrix
    pari.allocatemem(64000000,268435456,silent=True)
    data=next(c for c in r.read(half.INPUT)['cases'] if c['case_index']==index)
    old=next(c for c in r.read(old_artin.OUTPUT)['rows'] if c['case_index']==index)
    block=r.read(blocks.OUTPUT)['rows'][index]
    factor=r.read(rem.INPUT)['cases'][index]['factor'];assert factor['factorization_complete']
    ps=[p for p,e in factor['factors']]
    R=PolynomialRing(QQ,'z');f=R(data['integral_cubic_ascending'])
    nf=pari.nfinit([pari(f),ps]);assert str(nf.disc())==data['field_discriminant']
    assert list(map(str,nf.nf_get_zk()))==data['maximal_order_basis']
    theta=pari.Mod('z',pari(f));points=data['points'];n=len(points)
    masks=[c['point_mask'] for c in data['half_ideals']]
    m=int(block['generic_dimension']);assert m in (16,17)
    assert all(P['input_index']<m for P in points[:m]) and all(P['input_index']>=m for P in points[m:])
    ideals=[(p,j,P) for p in ps for j,P in enumerate(pari.idealprimedec(nf,p))]
    mat=lambda M:[[str(M[i,j]) for j in range(3)] for i in range(3)]
    columns=[]
    for j,point in enumerate(points):
        a,b,d=[ZZ(point[k]) for k in ('a','b','d')];gamma=pari(a)-pari(d*d)*theta
        assert pari.nfeltnorm(nf,gamma)==b*b
        I=pari.idealadd(nf,pari(b),gamma);assert mat(I)==point['gcd_ideal_hnf']
        J=I;removed=[]
        for p,k,P in ideals:
            e=int(pari.idealval(nf,I,P));assert e>=0
            removed.append([p,k,e])
            if e:J=pari.idealmul(nf,J,pari.idealpow(nf,P,-e))
        J=pari.idealhnf(nf,J);N=ZZ(pari.idealnorm(nf,J))
        assert N>0 and N%2 and all(N%p for p in ps)
        assert J[0,0]==N and J[1,1]==J[2,2]==1
        coords=pari.nfalgtobasis(nf,theta)
        root=ZZ(coords[0]-J[0,1]*coords[1]-J[0,2]*coords[2])%N
        assert f(root)%N==0 and ZZ(f.derivative()(root)).gcd(N)==1
        # Check that J is a half ideal outside exactly the retained support.
        quotient=pari.idealdiv(nf,pari.idealpow(nf,J,2),gamma)
        for p,k,P in ideals:
            e=int(pari.idealval(nf,quotient,P))
            if e:quotient=pari.idealmul(nf,quotient,pari.idealpow(nf,P,-e))
        assert pari.idealhnf(nf,quotient)==pari.idealhnf(nf,1)
        evaluations=[];bits=[]
        for P in points:
            ai,di=ZZ(P['a']),ZZ(P['d']);value=(ai-di*di*root)%N
            remaining=N
            while True:
                g=remaining.gcd(value)
                if g==1:break
                remaining//=g
            exceptional=N//remaining
            assert remaining.gcd(exceptional)==1
            derivative=3*ai*ai+ZZ(f[1])*di**4
            assert value.gcd(remaining)==1 and derivative.gcd(exceptional)==1 and di.gcd(exceptional)==1
            symbol=int(pari.kronecker(value,remaining))*int(pari.kronecker(derivative,exceptional))
            assert symbol in (-1,1);bits.append(int(symbol==-1))
            evaluations.append({'residue':str(value),'unit_modulus':str(remaining),'derivative_modulus':str(exceptional),
                'derivative_numerator':str(derivative),'symbol':symbol})
        strict=[(r.pack(bits)&mask).bit_count()%2 for mask in masks]
        columns.append({'point_position':j,'coprime_hnf':mat(J),'norm':str(N),'root_residue':str(root),
            'removed_S_factors':removed,'point_evaluations':evaluations,'strict_character_bits':strict,'strict_column':r.pack(strict)})
    pointcols=[c['strict_column'] for c in columns]
    def xor_word(mask):
        ans=0
        for j,c in enumerate(pointcols):
            if mask>>j&1:ans^=c
        return ans
    strictcols=[xor_word(mask) for mask in masks];assert strictcols==old['columns']
    generic=pointcols[:m];rg=r.rank(generic);rv=r.rank(strictcols);joint=r.rank(generic+strictcols)
    return {'case_index':index,'id':data['id'],'generic_dimension':m,'known_point_dimension':n,
        'strict_dimension':len(masks),'strict_masks':masks,'columns':columns,
        'strict_half_ideal_columns':strictcols,'generic_half_ideal_columns':generic,
        'generic_detected_rank':rg,'strict_detected_rank':rv,'joint_detected_rank':joint,
        'strict_dimensions_outside_generic_plus_S_units_lower_bound':joint-rg,
        'status':'PASS','bindings':bindings(),
        'boundary':r.read(PROTOCOL)['boundary']}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[]
    for index in r.read(PROTOCOL)['cases']:
        path=WORK/f'case-{index}.json';log=WORK/f'case-{index}.log'
        if not path.exists():
            with log.open('x') as stream:
                try:
                    proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--index',str(index)],stdout=stream,stderr=stream,timeout=60)
                    error=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:error='60-second timeout'
            if error:r.write_new(path,{'case_index':index,'status':'UNKNOWN','reason':error,'bindings':bindings()})
        row=r.read(path);assert row['bindings']==bindings();rows.append(row)
        print(index,row['status'],row.get('generic_detected_rank'),row.get('strict_dimensions_outside_generic_plus_S_units_lower_bound'),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.generic-sunit-carrier.v1','rows':rows,'bindings':bindings(),
        'provenance':r.read(PROTOCOL)['provenance']})


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['capture','worker']);parser.add_argument('--index',type=int);args=parser.parse_args()
    if args.mode=='worker':r.write_new(WORK/f'case-{args.index}.json',worker(args.index))
    else:capture()
