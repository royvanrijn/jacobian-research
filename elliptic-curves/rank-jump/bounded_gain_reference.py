#!/usr/bin/env python3
"""Supplement the frozen panel with one conditionally uncensored reference."""
import argparse
from pathlib import Path
import subprocess
import sys
from math import prod
import retrospective as r
import fresh_governing_panel as base
import fresh_governing_octics as octics
import fresh_strict_boundary_coordinates as coords

PROTOCOL=Path(__file__).with_name('BOUNDED_GAIN_REFERENCE_PROTOCOL.json')
INPUT=r.OUT/'rank_jump_bounded_gain_reference_inputs_v1.json'
PROVENANCE=r.OUT/'rank_jump_bounded_gain_reference_provenance_v1.json'
OUTPUT=r.OUT/'rank_jump_bounded_gain_reference_v1.json'
REPORT=r.OUT/'rank_jump_bounded_gain_reference_comparison_v1.json'
WORK=r.ROOT/'artifacts/local/rank-jump-bounded-gain-reference-v1'
LABEL=r.OUT/'small_conductor_class_completion_v1.json'
EQUATION=r.OUT/'small_conductor_rank22_proof_v1.json'
TOKEN='reference-00'


def export():
    row=next(x for x in r.read(r.INPUT)['rows'] if x['id']=='a1-fibration-05:3/17')
    assert len(row['generic_points'])==16
    r.short(row['model'],row['generic_points'])
    # This projection deliberately discards all exceptional coordinates and labels.
    hints=sorted({2,3,*[int(p) for p,e in r.read(EQUATION)['discriminant_factorization']]})
    r.write_new(INPUT,{'schema':'rank-jump.bounded-gain-reference-inputs.v1',
        'cases':[{'token':TOKEN,'model':row['model'],'generic_sections':row['generic_points']}],
        'equation_factor_hints':hints})
    r.write_new(PROVENANCE,{'schema':'rank-jump.bounded-gain-reference-provenance.v1',
        'id':row['id'],'family':row['family'],'parameter':row['parameter'],
        'projection':{'generic':'rows[id=a1-fibration-05:3/17]:model,generic_points',
                      'factor_hints':'small_conductor_rank22_proof_v1.json:discriminant_factorization.primes'},
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in (r.INPUT,EQUATION,INPUT,PROTOCOL)},
        'boundary':'Retrospective selection and source whitelist only. Not read by arithmetic workers.'})


def configure():
    base.INPUT=INPUT; base.WORK=WORK; base.PROTOCOL=PROTOCOL
    octics.PROTOCOL=PROTOCOL


def bindings():
    paths=(Path(__file__),PROTOCOL,INPUT,Path(base.__file__),Path(octics.__file__),
           Path(coords.__file__),Path(coords.prior.__file__),Path(r.__file__),base.LOCAL,
           Path(__file__).with_name('local_collision.py'))
    return {str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in paths}


def factor(token):
    from sage.all import ZZ
    f,pts,scale=base.model_data(token);D=ZZ(16*f.discriminant());rem=abs(D);fac=[]
    for p in r.read(INPUT)['equation_factor_hints']:
        assert ZZ(p).is_prime(proof=True)
        e=0
        while rem%p==0:rem//=p;e+=1
        if e:fac.append([p,e])
    assert prod(p**e for p,e in fac)*rem==abs(D)
    return {'status':'PASS' if rem==1 else 'UNKNOWN','factors':fac,'unresolved_cofactor':str(rem),
            'elliptic_discriminant':str(D),'integral_cubic_ascending':list(map(str,f.list())),'scale':str(scale)}


def setup(token):
    from sage.all import pari,ZZ
    fact=r.read(WORK/f'{token}-factor.json');loc=r.read(WORK/f'{token}-local.json')
    assert fact['status']=='PASS' and loc['status'] in ('PASS','PARTIAL')
    f,pts,scale=base.model_data(token);primes=[p for p,e in fact['factors']]
    assert prod(p**e for p,e in fact['factors'])==abs(16*f.discriminant())
    assert all(ZZ(p).is_prime(proof=True) for p in primes)
    return {'factor':fact,'local':loc},f,pts,primes,pari.nfinit([pari(f),primes])


def worker(stage):
    from sage.all import pari
    configure();pari.allocatemem(64000000,r.read(PROTOCOL)['limits']['pari_stack_bytes'],silent=True)
    coords.prior.setup=setup
    result={'factor':factor,'local':base.local_worker,'boundary':coords.worker,'octic':octics.worker}[stage](TOKEN)
    r.write_new(WORK/f'{TOKEN}-{stage}.json',{'bindings':bindings(),**result})


def capture():
    WORK.mkdir(parents=True,exist_ok=True);stages={}
    for stage in ('factor','local','boundary','octic'):
        path=WORK/f'{TOKEN}-{stage}.json'
        if not path.exists():
            reason=None
            if stage in ('local','boundary') and stages['factor']['status']!='PASS':reason='Incomplete factor coverage'
            if stage=='boundary' and stages['local']['status'] not in ('PASS','PARTIAL'):reason='Incomplete local coverage'
            if reason is None:
                with (WORK/f'{stage}.log').open('x') as log:
                    try:
                        p=subprocess.run([sys.executable,str(Path(__file__).resolve()),'worker','--stage',stage],
                            stdout=log,stderr=log,timeout=r.read(PROTOCOL)['limits'][stage+'_seconds'])
                        if p.returncode:reason='Worker failure; see retained checkpoint log'
                    except subprocess.TimeoutExpired:reason='Bounded worker timeout'
            if reason:r.write_new(path,{'status':'UNKNOWN','reason':reason,'bindings':bindings()})
        value=r.read(path);assert value['bindings']==bindings();stages[stage]=value
        print(stage,value['status'],value.get('strict_generic_dimension',value.get('additional_boundary_capacity_upper_bound')),flush=True)
    r.write_new(OUTPUT,{'schema':'rank-jump.bounded-gain-reference.v1','token':TOKEN,'stages':stages,'bindings':bindings()})


def report():
    d=r.read(OUTPUT)['stages'];label=r.read(LABEL);provenance=r.read(PROVENANCE)
    assert label['status']=='PASS'
    s=label['conditional_on_grh_selmer_dimension'];R=label['conditional_on_grh_exact_rank']
    assert s==R==label['unconditional_rank_lower_bound']==22
    assert label['conditional_on_grh_class_two_rank']==16
    result={'id':provenance['id'],'assumption':label['assumption'],
        'label_rank_lower_bound_unconditional':R,'label_rank_exact_conditional':R,
        'label_selmer_dimension_exact_conditional':s,'generic_dimension':16,
        'quotient_rank_exact_conditional':R-16,'label_Sha2_dimension_conditional':s-R,
        'additional_quotient_CT_computed_point_blind':'UNKNOWN',
        'original_full_CT_zero_conditional_inference':True,
        'class_group_label_used_as_arithmetic_input':False,
        'same_family_matched_panel_replacement':False,
        'inherited_governing_field_degree':d['octic'].get('governing_field_degree'),
        'inherited_minus_twist_CT_rank':d['local'].get('minus_twist_CT_rank')}
    if d['boundary']['status']=='PASS':
        b=d['boundary'];k=b['generic_strict_dimension'];h=b['Selmer_boundary_dimension_interval'][1];a=b['additional_boundary_capacity_upper_bound']
        # Complete conditional Selmer dimension determines c_S+dim(boundary), not c_S alone.
        cmin=max(k,s-h);cmax=min(label['conditional_on_grh_class_two_rank'],s-(16-k))
        assert cmin<=cmax
        result.update(generic_strict_dimension=k,local_product_dimension=b['local_point_product_dimension'],
            boundary_upper_bound=h,additional_boundary_capacity_upper_bound=a,
            additional_strict_rational_dimension_unconditional_lower_bound=max(0,R-h-k),
            localized_class_dimension_conditional_interval=[cmin,cmax],
            additional_strict_rational_dimension_conditional_interval=[cmin-k,cmax-k],
            exact_conditional_identity=f'(c_S - {k}) + e = {s-16}, 0 <= e <= {a}')
    files=(Path(__file__),PROTOCOL,INPUT,OUTPUT,PROVENANCE,LABEL,
           r.ROOT/'elliptic-curves/notes/SMALL_CONDUCTOR_CLASS_COMPLETION_PROOF_2026-09-06.md')
    r.write_new(REPORT,{'schema':'rank-jump.bounded-gain-reference-comparison.v1','result':result,
        'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in files},
        'boundary':'The complete rank/Selmer label uses point-derived proof data and GRH. This joins an external theorem; it does not independently compute the additional Selmer basis or its CT matrix. Sha[2]=0 is conditional; no assertion about other primes.'})
    print(result,flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['export','capture','worker','report']);p.add_argument('--stage');args=p.parse_args()
    if args.mode=='worker':worker(args.stage)
    else:globals()[args.mode]()
