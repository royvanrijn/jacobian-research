#!/usr/bin/env python3
"""Bounded exact local squareclass comparison at fixed odd primes."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import retrospective as r
import local_collision as lc

HERE=Path(__file__).resolve().parent
PROTOCOL=HERE/'BAD_PRIME_PROTOCOL.json'
INPUT=r.OUT/'rank_jump_bad_prime_inputs_v1.json'
OUTPUT=r.OUT/'rank_jump_bad_prime_support_v1.json'
CHECKPOINT=r.ROOT/'artifacts/local/rank-jump-bad-prime-v1'
LOCAL_SOURCE=r.ROOT/'elliptic-curves/cas/research_runtime/local_kummer.py'


def cases():
    inp=r.read(r.INPUT);indexed={x['id']:x for x in inp['rows']}
    return [indexed[pair[side]] for pair in inp['pairs'] for side in ('high','low')]


def bindings():
    return {'panel_input_sha256':r.digest(r.INPUT.read_bytes()),
      'panel_report_sha256':r.digest(r.RESULT.read_bytes()),
      'protocol_sha256':r.digest(PROTOCOL.read_bytes()),
      'script_sha256':r.digest(Path(__file__).read_bytes()),
      'local_squareclass_source_sha256':r.digest(LOCAL_SOURCE.read_bytes())}


def compute_case(index):
    from sage.all import QQ,ZZ,PolynomialRing,pari,lcm
    from sage.version import version
    sys.path.insert(0,str(r.ROOT/'elliptic-curves/cas'))
    from research_runtime.local_kummer import LocalSquareclasses
    row=cases()[index];profile,_,_=r.characterize(row);primes=r.read(PROTOCOL)['primes']
    model,allpoints=r.short(row['model'],row['generic_points']+row['points'])
    selected=profile['independent_input_indices'];points=[allpoints[i] for i in selected]
    A,B=map(QQ,model[3:]);d=lcm(A.denominator(),B.denominator())
    Ai,Bi=A*d**4,B*d**6
    assert Ai.denominator()==Bi.denominator()==1
    R=PolynomialRing(QQ,'z');theta=R.gen();pol=theta**3+Ai*theta+Bi
    # Deliberately no prepared_nf/SageArithmetic adapter: partial local order,
    # no full factorization, shared cache, or global maximal-order assertion.
    nf=pari.nfinit([pari(pol),primes])
    th=pari.Mod(pari(theta),pari(pol));xs=[QQ(P[0])*d**2 for P in points]
    betas=[pari(x)-th for x in xs]
    E=pari.ellinit([0,0,0,Ai,Bi]);local=[]
    for p in primes:
        chars=LocalSquareclasses(nf,p)
        sigs=[list(map(int,chars.signature(beta))) for beta in betas]
        decomposition=[{'ramification_index':int(P[2]),'residue_degree':int(P[3])} for P in chars.primes]
        assert sum(x['ramification_index']*x['residue_degree'] for x in decomposition)==3
        red=pari.elllocalred(E,p)
        minimal=pari.ellchangecurve(E,red[2])
        local.append({'prime':p,'prime_decomposition':decomposition,
          'character_order':'valuation parity, residue-square bit, for each prime ideal in PARI order',
          'point_signature_rows':sigs,'point_kummer_dimension':chars.point_kummer_dimension,
          'local_reduction':{'conductor_exponent':int(red[0]),'kodaira_code':int(red[1]),
             'minimal_change':[str(x) for x in red[2]],'tamagawa_number':int(red[3]),
             'minimal_discriminant_valuation':int(pari.valuation(minimal.disc(),p))},
          'uniformizers':[str(x[1]) for x in chars.data]})
    return {'schema':'rank-jump.bad-prime-case.v1','bindings':bindings(),'case_index':index,'id':row['id'],
       'software':{'sage':version,'pari':str(pari('version()'))},
       'selected_input_indices':selected,'generic_dimension':profile['generic_kummer_dimension_exact'],
       'witness_dimension':len(points),'short_model':model,'elliptic_scaling_d':str(d),
       'integral_cubic_ascending':[str(pol[i]) for i in range(4)],
       'local_order_basis':[str(x) for x in nf.nf_get_zk()],
       'local_order_scope':'maximal at 3,5,7,11,13 only; no global maximality claim','local':local}


def capture():
    CHECKPOINT.mkdir(parents=True,exist_ok=True);records=[]
    for i,row in enumerate(cases()):
        path=CHECKPOINT/f'case-{i}.json'
        if not path.exists():
            result=subprocess.run(['sage','-python',str(Path(__file__).resolve()),'worker','--case',str(i),
              '--destination',str(path)],cwd=r.ROOT,capture_output=True,text=True,timeout=30)
            log=CHECKPOINT/f'case-{i}.log'
            with log.open('x') as f:f.write(result.stdout+result.stderr)
            if result.returncode or not path.exists():raise RuntimeError(f'case {i} incomplete; see {log}')
        record=r.read(path);assert record['bindings']==bindings()
        records.append(record);print('checkpoint',i,row['id'],flush=True)
    r.write_new(INPUT,{'schema':'rank-jump.bad-prime-inputs.v1','bindings':bindings(),'cases':records})


def characterize(record):
    m=record['generic_dimension'];n=record['witness_dimension'];rows=[];joint=[0]*n;width=0
    for p in record['local']:
        sigs=list(map(r.pack,p['point_signature_rows']))
        dimension=p['point_kummer_dimension'];g=r.rank(sigs[:m]);allrank=r.rank(sigs)
        assert 0<=g<=allrank<=dimension
        for i,s in enumerate(sigs):joint[i]|=s<<width
        width+=len(p['point_signature_rows'][0])
        generic_basis=r.basis(sigs[:m]);residuals=[r.reduce(s,generic_basis) for s in sigs[m:]]
        assert r.rank(residuals)==allrank-g
        rows.append({'prime':p['prime'],'local_reduction':p['local_reduction'],
          'local_point_dimension':dimension,'generic_image_dimension':g,'witness_image_dimension':allrank,
          'quotient_support_dimension':allrank-g,'generic_image_is_full_point_image':g==dimension,
          'exceptional_residual_signatures':residuals})
    gr=r.rank(joint[:m]);ar=r.rank(joint)
    residual=[r.reduce(s,r.basis(joint[:m])) for s in joint[m:]]
    kernel=lc.orthogonal([r.pack([(s>>j)&1 for s in residual]) for j in range(width)],n-m)
    indices=[];images=[]
    for i,s in enumerate(joint[:m]):
        if r.rank(images+[s])>len(images):indices.append(i);images.append(s)
    corrections=[]
    for j,s in enumerate(joint[m:]):
        try:word=lc.coordinates(s,images)
        except ValueError:
            corrections.append({'quotient_index':j,'generic_correction_mask':None});continue
        mask=sum((1<<i) for bit,i in enumerate(indices) if word>>bit&1)
        assert lc.lift(mask,joint[:m])==s
        corrections.append({'quotient_index':j,'generic_correction_mask':mask})
    return {'id':record['id'],'generic_dimension':m,'witness_quotient_dimension':n-m,'local':rows,
       'joint_generic_image_dimension':gr,'joint_witness_image_dimension':ar,
       'joint_quotient_support_dimension':ar-gr,
       'full_product_point_image_dimension':sum(p['local_point_dimension'] for p in rows),
       'generic_image_is_full_product_point_image':gr==sum(p['local_point_dimension'] for p in rows),
       'generic_right_inverse_indices':indices,'simultaneous_generic_corrections':corrections,
       'joint_exceptional_residual_signatures':residual,
       'joint_exceptional_kernel_in_quotient_coordinates':kernel,
       'boundary':'Local image in this fixed prime dictionary only; quotient coordinates refer to the specified independent point basis.'}


def build(check=False):
    inp=r.read(INPUT);assert inp['bindings']==bindings()
    rows=[characterize(x) for x in inp['cases']]
    local=[p for x in rows for p in x['local']];bad=[p for p in local if p['local_reduction']['conductor_exponent']]
    result={'schema':'rank-jump.bad-prime-support.v1','input_sha256':r.digest(INPUT.read_bytes()),
      'bindings':bindings(),'rows':rows,
      'summary':{'cases':len(rows),'local_cases':len(local),'bad_reduction_cases':len(bad),
         'bad_cases_with_new_quotient_support':sum(p['quotient_support_dimension']>0 for p in bad),
         'all_cases_with_new_quotient_support':sum(p['quotient_support_dimension']>0 for p in local),
         'generic_full_local_point_images':sum(p['generic_image_is_full_point_image'] for p in local),
         'generic_full_product_point_images':sum(x['generic_image_is_full_product_point_image'] for x in rows),
         'simultaneous_exceptional_generic_corrections':sum(sum(y['generic_correction_mask'] is not None for y in x['simultaneous_generic_corrections']) for x in rows),
         'cases_with_joint_quotient_support':sum(x['joint_quotient_support_dimension']>0 for x in rows)},
      'boundary':'Exact local class incidence for known retrospective point subgroups. No full Selmer space, global solubility criterion, rank upper bound or prospective score.'}
    if check:
        if r.read(OUTPUT)!=result:raise ValueError('bad-prime replay mismatch')
        print('PASS bad-prime binary replay')
    else:r.write_new(OUTPUT,result)
    print(result['summary'])
    for x in rows:
        print(x['id'],'local',[(p['prime'],p['generic_image_dimension'],p['local_point_dimension'],p['quotient_support_dimension']) for p in x['local']],
              'joint',x['joint_generic_image_dimension'],x['joint_witness_image_dimension'])


def verify():
    from sage.all import pari,QQ,PolynomialRing
    inp=r.read(INPUT)
    for record in inp['cases']:
        if compute_case(record['case_index'])!=record:raise ValueError('local arithmetic replay mismatch')
        profile=characterize(record);d=QQ(record['elliptic_scaling_d'])
        row=cases()[record['case_index']]
        _,allpoints=r.short(row['model'],row['generic_points']+row['points'])
        points=[allpoints[i] for i in record['selected_input_indices']]
        R=PolynomialRing(QQ,'z');pol=pari(R(list(map(QQ,record['integral_cubic_ascending']))))
        nf=pari.nfinit([pol,r.read(PROTOCOL)['primes']]);th=pari.Mod('z',pol)
        betas=[pari(QQ(P[0])*d*d)-th for P in points]
        prime_ideals=[P for p in r.read(PROTOCOL)['primes'] for P in pari.idealprimedec(nf,p)]
        for c in profile['simultaneous_generic_corrections']:
            mask=c['generic_correction_mask']
            if mask is None:continue
            beta=betas[record['generic_dimension']+c['quotient_index']]
            for i in range(record['generic_dimension']):
                if mask>>i&1:beta*=betas[i]
            # PARI 2.17.3 can reject nonintegral basis columns here. Clear
            # their denominators by a rational square, preserving the class.
            den=pari.denominator(pari.nfalgtobasis(nf,beta));beta*=den**2
            assert all(pari.nfislocalpower(nf,P,beta,2)==1 for P in prime_ideals)
        print('PASS local arithmetic',record['id'],flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=('worker','capture','build','check','verify'))
    p.add_argument('--case',type=int);p.add_argument('--destination',type=Path);a=p.parse_args()
    if a.mode=='worker':r.write_new(a.destination,compute_case(a.case))
    elif a.mode=='capture':capture()
    elif a.mode=='verify':verify()
    else:build(a.mode=='check')
