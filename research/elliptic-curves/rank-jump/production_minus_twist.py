#!/usr/bin/env python3
"""Bounded local -1 twist images, using a previously complete global boundary."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import local_collision as lc
import strict_class_blocks as strict
import remaining_bad_primes as rem
import derivative_local_duality as derivative
import scalar_cup

PROTOCOL=Path(__file__).with_name('PRODUCTION_MINUS_TWIST_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_production_minus_twist_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_production_minus_twist_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-production-minus-twist-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
            (Path(__file__), PROTOCOL, derivative.INPUT, derivative.OUTPUT,
             strict.OUTPUT, rem.INPUT, rem.OUTPUT, rem.bad.INPUT, rem.dy.INPUT,
             r.INPUT, rem.bad.LOCAL_SOURCE, scalar_cup.OUTPUT)}


def intersection(a,b):
    """Basis of the intersection in the same ambient bit coordinates."""
    a=list(r.basis(a).values()); bb=r.basis(b)
    residual=[r.reduce(x,bb) for x in a]
    relations=lc.orthogonal([r.pack((x>>j)&1 for x in residual)
                             for j in range(max(residual,default=0).bit_length())],len(a))
    out=[lc.lift(v,a) for v in relations]
    assert r.rank(out)==len(a)+len(bb)-r.rank(a+list(bb.values()))
    return out


def local_square(value,p):
    """Exact Q_p criterion; output retains a valuation and unit witness."""
    from sage.all import QQ
    value=QQ(value)
    if not value:return None
    valuation=int(value.valuation(p)); unit=value/QQ(p)**valuation
    modulus=8 if p==2 else p
    residue=r.mod(r.F(str(unit)),modulus)
    square=valuation%2==0 and (residue==1 if p==2 else pow(residue,(p-1)//2,p)==1)
    return {'valuation':valuation,'unit_residue':int(residue),'modulus':int(modulus),'square':square}


def compute(index):
    from sage.all import QQ, ZZ, AA, GF, PolynomialRing, pari
    from sage.version import version
    sys.path.insert(0,str(rem.bad.LOCAL_SOURCE.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    spec=r.read(PROTOCOL)['limits']
    old=r.read(rem.bad.INPUT)['cases'][index]
    source=rem.bad.cases()[index]
    block=r.read(strict.OUTPUT)['rows'][index]
    assert block['all_bad_places_complete']
    primes=[p for p in block['tested_places'] if p!='infinity']
    _,allpoints=r.short(source['model'],source['generic_points']+source['points'])
    points=[allpoints[i] for i in old['selected_input_indices']]
    n,m=len(points),old['generic_dimension']
    R=PolynomialRing(QQ,'z'); z=R.gen()
    f=R(list(map(QQ,old['integral_cubic_ascending'])))
    assert f[2]==0 and f[3]==1
    factor=r.read(rem.INPUT)['cases'][index]['factor']
    assert factor['factorization_complete'] and 16*f.discriminant()==ZZ(factor['model_discriminant'])
    assert {p for p,e in factor['factors']}<=set(primes)
    nf=pari.nfinit([pari(f),primes]); theta=pari.Mod('z',pari(f))
    scale=QQ(old['elliptic_scaling_d'])
    betas=[pari(QQ(P[0])*scale**2)-theta for P in points]
    derivative_beta=-pari(f.discriminant())*pari(f.derivative())(theta)
    betas.append(derivative_beta)
    twist=z**3+f[1]*z-f[0]
    E=pari.ellinit([0,0,0,f[1],-f[0]])
    roots=f.roots(AA,multiplicities=False)
    locals=[];joint=[0]*(n+1);twist_product=[];original_product=[];offset=0
    for local in strict.local_rows(index):
        p=local['prime'];dim=local['point_kummer_dimension']; witnesses=[];trials=0
        if p=='infinity':
            signatures=[r.pack(int(QQ(P[0])*scale**2<a) for a in roots) for P in points]
            signatures.append(r.pack(int(-f.discriminant()*f.derivative()(a)<0) for a in roots))
            width=len(roots); basis=[3] if dim else []
            method='ordered real roots: twist bounded component has signs (-,-,+)'
        else:
            chars=LocalSquareclasses(nf,p)
            signatures=[r.pack(chars.signature(beta)) for beta in betas]
            width=len(chars.signature(derivative_beta))
            assert dim==chars.point_kummer_dimension
            basis=[]
            if p%4==1:
                basis=list(r.basis(signatures[:m]).values())
                method='-1 is a square in Q_p; labelled local conditions coincide'
            else:
                method='exact local-square x witnesses spanning the dimension bound'
                red=pari.elllocalred(E,p); u,shift=map(QQ,list(red[2])[:2])
                g=R(twist(u*u*z+shift)/u**6)
                def candidates():
                    for t in range(spec['small_integer_radius']+1):
                        yield QQ(t)
                        if t:yield QQ(-t)
                    if p==2:
                        for t in range(1,2*spec['small_integer_radius'],2):yield QQ(t)/4
                    else:
                        residue_poly=PolynomialRing(GF(p),'t')([GF(p)(c) for c in g])
                        for root in residue_poly.roots(multiplicities=False):
                            center=ZZ(root)
                            for depth in range(1,spec['root_depth']+1):
                                for shift_bit in spec['root_offsets']:
                                    yield QQ(center+shift_bit*ZZ(p)**depth)
                seen=set()
                for t in candidates():
                    if len(basis)==dim:break
                    if t in seen:continue
                    seen.add(t);trials+=1
                    if trials>spec['max_candidates_per_place']:break
                    x=u*u*t+shift; value=twist(x); square=local_square(value,p)
                    if square is None or not square['square']:continue
                    signature=r.pack(chars.signature(pari(x)+theta))
                    if r.rank(basis+[signature])>len(basis):
                        basis.append(signature)
                        witnesses.append({'minimal_x':str(t),'short_x':str(x),
                            'cubic_value':str(value),'square_witness':square,'signature':signature})
        original=list(r.basis(signatures[:m]).values())
        assert len(original)==dim
        complete=len(basis)==dim
        for j,sig in enumerate(signatures):joint[j]|=sig<<offset
        original_product.extend(v<<offset for v in original)
        twist_product.extend(v<<offset for v in basis)
        locals.append({'place':p,'width':width,'point_dimension':dim,'complete':complete,
            'method':method,'trials':trials,'witnesses':witnesses,'original_basis':original,
            'twist_basis':basis,'global_signatures':signatures,
            'intersection_dimension':len(intersection(original,basis)) if complete else None})
        offset+=width
        print('place',p,'dimension',len(basis),'/',dim,'trials',trials,flush=True)
    ell=sum(x['point_dimension'] for x in locals)
    assert r.rank(joint[:-1])==ell-1 and r.rank(joint)==ell
    assert len(intersection(joint,original_product))==ell-1
    return {'schema':'rank-jump.production-minus-twist-case.v1','bindings':bindings(),
        'case_index':index,'id':old['id'],'software':{'sage':version,'pari':str(pari.version())},
        'integral_cubic_ascending':old['integral_cubic_ascending'],
        'generic_dimension':m,'witness_dimension':n,
        'known_strict_dimension':len(block['witness_strict_kernel_masks']),
        'local':locals,'relaxed_boundary_dimension':ell,'joint_global_signatures':joint,
        'original_product_basis':original_product,'twist_product_basis':twist_product,
        'complete':all(x['complete'] for x in locals)}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);cases=[]
    for i in r.read(PROTOCOL)['cases']:
        path=WORK/f'case-{i}.json'
        if not path.exists():
            with (WORK/f'case-{i}.log').open('x') as log:
                try:
                    proc=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker',
                        '--index',str(i),'--destination',str(path)],cwd=r.ROOT,stdout=log,stderr=log,
                        timeout=r.read(PROTOCOL)['limits']['seconds_per_case'])
                    failure=None if proc.returncode==0 else 'worker failure'
                except subprocess.TimeoutExpired:failure='60-second timeout'
                if failure and not path.exists():
                    r.write_new(path,{'bindings':bindings(),'case_index':i,'status':'UNKNOWN','reason':failure})
        row=r.read(path);assert row['bindings']==bindings();cases.append(row)
        print('checkpoint',i,row.get('complete',row.get('status')),flush=True)
    r.write_new(INPUT,{'schema':'rank-jump.production-minus-twist-inputs.v1','bindings':bindings(),'cases':cases})


def build(check=False):
    data=r.read(INPUT);assert data['bindings']==bindings();rows=[]
    for row in data['cases']:
        if not row.get('complete'):
            rows.append({'case_index':row['case_index'],'status':'UNKNOWN'});continue
        boundary=row['joint_global_signatures'];ell=row['relaxed_boundary_dimension']
        assert r.rank(boundary)==ell
        b=len(intersection(boundary,row['twist_product_basis']))
        common=len(intersection(intersection(boundary,row['original_product_basis']),row['twist_product_basis']))
        n=row['witness_dimension'];k=row['known_strict_dimension']
        assert n==k+ell-1
        rows.append({'case_index':row['case_index'],'id':row['id'],'status':'PASS',
            'generic_dimension':row['generic_dimension'],'known_rank_lower_bound':n,
            'known_strict_dimension':k,'original_boundary_dimension':ell-1,
            'twist_boundary_dimension':b,'common_boundary_dimension':common,
            'full_original_Selmer_dimension':f'{n} + epsilon',
            'full_twist_Selmer_dimension':f'{k+b} + epsilon',
            'full_common_Selmer_dimension':f'{k+common} + epsilon',
            'full_Selmer_dimension_change':b-(ell-1),
            'epsilon_definition':f'dim Cl(O_K,S)/2 - {k} >= 0, unchanged under this twist',
            'changed_local_places':[x['place'] for x in row['local']
                                    if x['intersection_dimension']<x['point_dimension']]})
    result={'schema':'rank-jump.production-minus-twist.v1','bindings':bindings(),
            'input_sha256':r.digest(INPUT.read_bytes()),'rows':rows}
    if check:assert r.read(OUTPUT)==result;print('PASS boundary intersections')
    else:r.write_new(OUTPUT,result)
    for row in rows:print(row)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['capture','worker','build','check'])
    parser.add_argument('--index',type=int);parser.add_argument('--destination',type=Path);args=parser.parse_args()
    if args.mode=='capture':capture()
    elif args.mode=='worker':r.write_new(args.destination,compute(args.index))
    else:build(args.mode=='check')
