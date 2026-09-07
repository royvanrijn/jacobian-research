#!/usr/bin/env python3
"""Exact dyadic and real completion of the retained small-place comparison."""
import argparse
from pathlib import Path
import subprocess
import sys
import retrospective as r
import bad_prime_support as bad

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'DYADIC_REAL_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_dyadic_real_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_dyadic_real_support_v1.json'
CHECKPOINT=r.ROOT/'artifacts/local/rank-jump-dyadic-real-v1'


def bindings():
    return {'odd_prime_input_sha256':r.digest(bad.INPUT.read_bytes()),
      'odd_prime_report_sha256':r.digest(bad.OUTPUT.read_bytes()),
      'odd_prime_script_sha256':r.digest(Path(bad.__file__).read_bytes()),
      'panel_input_sha256':r.digest(r.INPUT.read_bytes()),
      'local_squareclass_source_sha256':r.digest(bad.LOCAL_SOURCE.read_bytes()),
      'protocol_sha256':r.digest(PROTOCOL.read_bytes()),'script_sha256':r.digest(Path(__file__).read_bytes())}


def real_bit(A,B,x):
    A,B,x=(r.F(str(v)) for v in (A,B,x));disc=-4*A**3-27*B*B
    if not disc:raise ValueError('singular cubic')
    if disc<0:return 0
    assert A<0
    return int(x<0 or 3*x*x<-A)


def compute_case(index):
    from sage.all import QQ,PolynomialRing,pari
    from sage.version import version
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    from research_runtime.local_kummer import LocalSquareclasses
    old=r.read(bad.INPUT)['cases'][index];row=bad.cases()[index]
    model,allpoints=r.short(row['model'],row['generic_points']+row['points'])
    assert old['id']==row['id'] and old['short_model']==model
    points=[allpoints[i] for i in old['selected_input_indices']]
    R=PolynomialRing(QQ,'z');pol=pari(R(list(map(QQ,old['integral_cubic_ascending']))))
    nf=pari.nfinit([pol,[2]]);theta=pari.Mod('z',pol);d=QQ(old['elliptic_scaling_d'])
    betas=[pari(QQ(P[0])*d*d)-theta for P in points]
    chars=LocalSquareclasses(nf,2);sigs=[list(map(int,chars.signature(beta))) for beta in betas]
    data=[]
    for P,uniformizer,bid in chars.data:
        e,f=int(P[2]),int(P[3])
        data.append({'ramification_index':e,'residue_degree':f,'unit_modulus_exponent':2*e+1,
          'unit_cyclic_factors':[int(x) for x in bid.bid_get_cyc()],'uniformizer':str(uniformizer)})
    assert sum(x['ramification_index']*x['residue_degree'] for x in data)==3
    assert all(len(s)==3+2*len(data) for s in sigs)
    E=pari.ellinit([0,0,0,QQ(old['integral_cubic_ascending'][1]),QQ(old['integral_cubic_ascending'][0])])
    red=pari.elllocalred(E,2);minimal=pari.ellchangecurve(E,red[2])
    two={'prime':2,'point_signature_rows':sigs,'point_kummer_dimension':chars.point_kummer_dimension,
      'prime_decomposition':data,'character_order':'valuation parity followed by unit discrete logs modulo 2 for even cyclic factors, in PARI order',
      'local_reduction':{'conductor_exponent':int(red[0]),'kodaira_code':int(red[1]),
        'tamagawa_number':int(red[3]),'minimal_change':[str(x) for x in red[2]],
        'minimal_discriminant_valuation':int(pari.valuation(minimal.disc(),2))}}
    A,B=map(r.F,model[3:]);disc=-4*A**3-27*B*B;dim=int(disc>0)
    real={'prime':'infinity','point_kummer_dimension':dim,
      'point_signature_rows':[[real_bit(A,B,P[0])] if dim else [] for P in points],
      'local_reduction':{'cubic_discriminant_sign':1 if disc>0 else -1},
      'character_order':'bounded real component bit when the cubic has three real roots'}
    return {'schema':'rank-jump.dyadic-real-case.v1','bindings':bindings(),'case_index':index,
      'id':row['id'],'generic_dimension':old['generic_dimension'],'witness_dimension':old['witness_dimension'],
      'selected_input_indices':old['selected_input_indices'],'software':{'sage':version,'pari':str(pari('version()'))},
      'local_order_basis':[str(x) for x in nf.nf_get_zk()],
      'local_order_scope':'maximal at 2 only; no global maximality claim','local':[two,real]}


def capture():
    CHECKPOINT.mkdir(parents=True,exist_ok=True);records=[]
    for i,row in enumerate(bad.cases()):
        path=CHECKPOINT/f'case-{i}.json'
        if not path.exists():
            result=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--case',str(i),
              '--destination',str(path)],cwd=r.ROOT,capture_output=True,text=True,timeout=30)
            log=CHECKPOINT/f'case-{i}.log'
            with log.open('x') as f:f.write(result.stdout+result.stderr)
            if result.returncode or not path.exists():raise RuntimeError(f'case {i} incomplete; see {log}')
        record=r.read(path);assert record['bindings']==bindings()
        records.append(record);print('checkpoint',i,row['id'],flush=True)
    r.write_new(INPUT,{'schema':'rank-jump.dyadic-real-inputs.v1','bindings':bindings(),'cases':records})


def combine(record):
    old=r.read(bad.INPUT)['cases'][record['case_index']]
    assert old['id']==record['id'] and old['selected_input_indices']==record['selected_input_indices']
    assert old['generic_dimension']==record['generic_dimension'] and old['witness_dimension']==record['witness_dimension']
    return {**record,'local':[record['local'][0]]+old['local']+[record['local'][1]]}


def build(check=False):
    inp=r.read(INPUT);assert inp['bindings']==bindings()
    rows=[]
    for record in inp['cases']:
        added=bad.characterize(record);combined=bad.characterize(combine(record))
        rows.append({'id':record['id'],'new_place_comparison':added,'seven_place_comparison':combined})
    result={'schema':'rank-jump.dyadic-real-support.v1','input_sha256':r.digest(INPUT.read_bytes()),
      'bindings':bindings(),'rows':rows,
      'summary':{'cases':len(rows),'new_local_cases':sum(len(x['new_place_comparison']['local']) for x in rows),
        'new_places_with_quotient_support':sum(p['quotient_support_dimension']>0 for x in rows for p in x['new_place_comparison']['local']),
        'seven_place_cases_with_quotient_support':sum(x['seven_place_comparison']['joint_quotient_support_dimension']>0 for x in rows),
        'generic_full_seven_place_products':sum(x['seven_place_comparison']['generic_image_is_full_product_point_image'] for x in rows)},
      'boundary':'Only the declared seven places and known witness subgroups. No full Selmer calculation, additional point, rank upper bound, or prospective selector.'}
    if check:
        if r.read(OUTPUT)!=result:raise ValueError('dyadic-real binary replay mismatch')
        print('PASS dyadic-real binary replay')
    else:r.write_new(OUTPUT,result)
    print(result['summary'])
    for x in rows:
        print(x['id'],'new',[(p['prime'],p['generic_image_dimension'],p['local_point_dimension'],p['quotient_support_dimension']) for p in x['new_place_comparison']['local']],
          'joint',x['seven_place_comparison']['joint_generic_image_dimension'],x['seven_place_comparison']['joint_witness_image_dimension'],
          'product',x['seven_place_comparison']['full_product_point_image_dimension'])


def verify():
    from sage.all import pari,QQ,AA,PolynomialRing
    for record in r.read(INPUT)['cases']:
        if compute_case(record['case_index'])!=record:raise ValueError('dyadic arithmetic replay mismatch')
        old=r.read(bad.INPUT)['cases'][record['case_index']];row=bad.cases()[record['case_index']]
        model,allpoints=r.short(row['model'],row['generic_points']+row['points'])
        points=[allpoints[i] for i in record['selected_input_indices']];A,B=map(QQ,model[3:])
        R=PolynomialRing(QQ,'z');f=R([B,A,0,1]);roots=f.roots(AA,multiplicities=False)
        for P in points:
            x=QQ(P[0]);bit=int(len(roots)==3 and any(x<t for t in roots))
            assert bit==real_bit(A,B,x)
        pol=pari(R(list(map(QQ,old['integral_cubic_ascending']))))
        nf=pari.nfinit([pol,[2,3,5,7,11,13]]);th=pari.Mod('z',pol);d=QQ(old['elliptic_scaling_d'])
        betas=[pari(QQ(P[0])*d*d)-th for P in points]
        prime_ideals=[P for p in (2,3,5,7,11,13) for P in pari.idealprimedec(nf,p)]
        profile=bad.characterize(combine(record))
        for c in profile['simultaneous_generic_corrections']:
            mask=c['generic_correction_mask']
            if mask is None:continue
            j=record['generic_dimension']+c['quotient_index'];beta=betas[j]
            for i in range(record['generic_dimension']):
                if mask>>i&1:beta*=betas[i]
            den=pari.denominator(pari.nfalgtobasis(nf,beta));beta*=den**2
            assert all(pari.nfislocalpower(nf,P,beta,2)==1 for P in prime_ideals)
            for root in roots:
                signs=int(QQ(points[j][0])<root)
                for i in range(record['generic_dimension']):
                    if mask>>i&1:signs^=int(QQ(points[i][0])<root)
                assert signs==0
        print('PASS dyadic/real arithmetic and corrected local powers',record['id'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('worker','capture','build','check','verify'))
    p.add_argument('--case',type=int);p.add_argument('--destination',type=Path);a=p.parse_args()
    if a.mode=='worker':r.write_new(a.destination,compute_case(a.case))
    elif a.mode=='capture':capture()
    elif a.mode=='verify':verify()
    else:build(a.mode=='check')
