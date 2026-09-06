#!/usr/bin/env python3
"""Equation-defined boundary obstruction and reduced class-engine control."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import fresh_governing_panel as base

PROTOCOL=Path(__file__).with_name('MATCHED103B2_CLASS_BOUNDARY_PROTOCOL.json')
OUTPUT=r.OUT/'rank_jump_matched103b2_class_boundary_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-matched103b2-class-boundary-v1'


def bindings():
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (Path(__file__),PROTOCOL,base.INPUT,base.OUTPUT,Path(base.__file__),Path(r.__file__),base.LOCAL)}


def setup(token):
    from sage.all import ZZ,pari
    row=next(x for x in r.read(base.OUTPUT)['rows'] if x['token']==token)
    f,pts,_=base.model_data(token);factor=row['factor'];assert factor['status']=='PASS'
    product=ZZ(1)
    for p,e in factor['factors']:assert ZZ(p).is_prime(proof=True);product*=ZZ(p)**e
    assert product==abs(16*f.discriminant())
    primes=[p for p,e in factor['factors']];nf=pari.nfinit([pari(f),primes])
    return row,f,pts,primes,nf


def boundary(token):
    from sage.all import AA,pari,QQ,EllipticCurve,GF,matrix
    sys.path.insert(0,str(base.LOCAL.parents[1]))
    from research_runtime.local_kummer import LocalSquareclasses
    old,f,pts,primes,nf=setup(token);theta=pari.Mod('z',pari(f));delta=f.discriminant()
    beta=-pari(delta)*pari(f.derivative())(theta)
    assert pari.nfeltnorm(nf,beta)==delta**4 and delta>0
    E=EllipticCurve([0,0,0,f[1],f[0]])
    S=[2]+[p for p in primes if p!=2 and E.local_data(p,proof=True).conductor_valuation()>0]
    # No unsupported ramification assumption: all omitted polynomial-discriminant
    # primes have to be checked for derivative unramifiedness directly.
    for p in set(primes)-set(S):
        assert p!=2
        for P in pari.idealprimedec(nf,p):assert int(pari.idealval(nf,beta,P))%2==0
    joint=[0]*len(pts);offset=0;ell=0;locals=[]
    for p in S:
        chars=LocalSquareclasses(nf,p);sigs=[r.pack(chars.signature(pari(x)-theta)) for x,y in pts]
        ds=list(chars.signature(beta));width=len(ds);dim=chars.point_kummer_dimension
        assert r.rank(sigs)<=dim;ell+=dim
        for i,s in enumerate(sigs):joint[i]|=s<<offset
        offset+=width
        locals.append({'place':p,'width':width,'generic_signatures':sigs,'derivative_signature':r.pack(ds),'point_dimension':dim,'generic_local_dimension':r.rank(sigs)})
    roots=f.roots(AA,multiplicities=False);assert len(roots)==3
    signs=[int(-delta*f.derivative()(a)<0) for a in roots];assert signs==[1,0,1]
    sigs=[r.pack(int(x<a) for a in roots) for x,y in pts];assert set(sigs)<={0,6}
    for i,s in enumerate(sigs):joint[i]|=s<<offset
    locals.append({'place':'infinity','width':3,'generic_signatures':sigs,'derivative_signature':5,'point_dimension':1,'complete_real_point_basis':[6]})
    ell+=1;g=r.rank(joint);m=len(pts)
    assert g==m==17 and g==matrix(GF(2),[[v>>j&1 for j in range(offset+3)] for v in joint]).rank()
    # Real Hilbert pairing is the dot product of sign masks: 5 dot 6 =1.
    assert (5&6).bit_count()%2==1
    return {'status':'PASS','S_finite':S,'polynomial_discriminant':str(delta),
       'field_discriminant':str(nf.disc()),'field_discriminant_bits':int(abs(int(nf.disc())).bit_length()),
       'derivative_coefficients':[str(pari.lift(beta).polcoef(i)) for i in range(3)],
       'derivative_norm':str(delta**4),'local':locals,'generic_joint_signatures':joint,
       'generic_dimension':m,'generic_local_dimension':g,'generic_strict_dimension':m-g,
       'local_point_product_dimension':ell,'derivative_constraint_rank_lower_bound':1,
       'Selmer_boundary_dimension_lower_bound':g,'Selmer_boundary_dimension_upper_bound':ell-1,
       'additional_boundary_dimension_upper_bound':ell-1-g,
       'full_Selmer_dimension_formula':f'c_S + 17 + e, 0 <= e <= {ell-1-g}',
       'additional_Selmer_dimension_formula':f'c_S + e, 0 <= e <= {ell-1-g}',
       'c_S':'UNKNOWN: dim_F2 Cl(O_K,S)/2',
       'boundary':'Exact symbolic bounds from generic sections and a nonzero equation-defined reciprocity functional. No additional class representative supplied.'}


def reduction(token):
    from sage.all import pari,QQ,PolynomialRing,matrix
    old,f,pts,primes,nf=setup(token)
    reduced,image=pari.polredbest(nf,1);inverse=pari.modreverse(image)
    assert pari(f)(image)==0 and reduced(inverse)==0
    # PARI polynomials are callable even for constants; normalize every entry.
    R=PolynomialRing(QQ,'z')
    basis=[pari(R([QQ(pari.lift(b).polcoef(i)) for i in range(3)]))(image) for b in nf.nf_get_zk()]
    vectors=[[QQ(pari.lift(b).polcoef(i)) for i in range(3)] for b in basis]
    det=matrix(QQ,vectors).determinant();assert det*det*QQ(pari.poldisc(reduced))==QQ(nf.disc())
    assert pari.lift(image)(inverse)==pari.Mod('z',pari(f))
    return {'status':'PASS','reduced_cubic_ascending':[str(reduced.polcoef(i)) for i in range(4)],
        'original_root_in_reduced':[str(pari.lift(image).polcoef(i)) for i in range(3)],
        'reduced_root_in_original':[str(pari.lift(inverse).polcoef(i)) for i in range(3)],
        'transported_maximal_order_basis':[[str(x) for x in row] for row in vectors],
        'field_discriminant':str(nf.disc()),'original_coefficient_bits':max(abs(int(c)).bit_length() for c in f.list()),
        'reduced_coefficient_bits':max(abs(int(reduced.polcoef(i))).bit_length() for i in range(4)),
        'boundary':'Exact same cubic field. A smaller defining polynomial does not reduce the field discriminant or prove a class-group change.'}


def class_worker(token):
    from sage.all import pari,QQ,PolynomialRing
    rec=r.read(WORK/f'{token}-reduction.json')
    if rec['status']!='PASS':return {'status':'UNKNOWN','reason':'reduction failed'}
    R=PolynomialRing(QQ,'z');f=R(list(map(QQ,rec['reduced_cubic_ascending'])))
    basis=[pari(R(list(map(QQ,row)))) for row in rec['transported_maximal_order_basis']]
    nf=pari.nfinit([pari(f),basis]);assert str(nf.disc())==rec['field_discriminant']
    pari.setrand(20260906)
    try:
        bnf=pari.bnfinit(nf,1)
        assert pari.bnfcertify(bnf)==1
        return {'status':'PASS_CERTIFIED','class_group_cyclic':list(map(str,bnf.bnf_get_cyc())),
                'units':list(map(str,bnf.bnf_get_fu())),'class_generators':list(map(str,bnf.bnf_get_gen()))}
    except Exception as exc:
        return {'status':'UNKNOWN','exception_type':type(exc).__name__,'reason':str(exc)}


def capture():
    WORK.mkdir(parents=True,exist_ok=True);rows=[];lim=r.read(PROTOCOL)['limits']
    for token in r.read(PROTOCOL)['cases']:
        row={'token':token}
        for stage,key in [('boundary','boundary_seconds_per_case'),('reduction','reduction_seconds_per_case'),('class','class_seconds_per_case')]:
            path=WORK/f'{token}-{stage}.json'
            if not path.exists():
                with (WORK/f'{token}-{stage}.log').open('x') as log:
                    try:
                        proc=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--token',token,'--stage',stage],stdout=log,stderr=log,timeout=lim[key])
                        failure=None if proc.returncode==0 else 'worker failure'
                    except subprocess.TimeoutExpired:failure=f'{lim[key]}-second timeout'
                if failure:r.write_new(path,{'bindings':bindings(),'status':'UNKNOWN','reason':failure})
            record=r.read(path);assert record['bindings']==bindings();row[stage]=record
            print(token,stage,record['status'],record.get('additional_boundary_dimension_upper_bound'),flush=True)
        rows.append(row)
    r.write_new(OUTPUT,{'schema':'rank-jump.matched103b2-class-boundary.v1','bindings':bindings(),'rows':rows})


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['capture','worker']);p.add_argument('--token');p.add_argument('--stage');args=p.parse_args()
    if args.mode=='capture':capture()
    else:
        from sage.all import pari
        pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
        result={'boundary':boundary,'reduction':reduction,'class':class_worker}[args.stage](args.token)
        r.write_new(WORK/f'{args.token}-{args.stage}.json',{'bindings':bindings(),**result})
