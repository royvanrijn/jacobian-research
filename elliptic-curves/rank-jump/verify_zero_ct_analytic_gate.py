#!/usr/bin/env python3
"""Independent coefficients, Tate data and MPFR derivative-error replay."""
import argparse
import json
from pathlib import Path
import retrospective as r
import zero_ct_analytic_gate as source

OUTPUT=r.OUT/'rank_jump_zero_ct_analytic_gate_verification_v1.json'


def compute():
    from sage.all import QQ,ZZ,EllipticCurve
    import inspect
    from sage.schemes.elliptic_curves.lseries_ell import Lseries_ell
    data=r.read(source.OUTPUT);prior=r.read(source.PRIOR);spec=r.read(source.PROTOCOL)
    assert data['status']==prior['status']=='PASS'
    for doc in (data,prior):
        for path,sha in doc['bindings'].items():assert r.digest((r.ROOT/path).read_bytes())==sha
    maximum=max(x['terms'] for x in data['rows'])
    base=EllipticCurve(QQ,spec['base_model']);an=list(map(int,base.anlist(maximum,python_ints=True)))
    small_trace_checks=[]
    for p in r.primes(199):
        def character(v):
            v%=p
            if not v:return 0
            return 1 if pow(v,(p-1)//2,p)==1 else -1
        ap=-sum(character(x**3-11*x*x-14*x-1) for x in range(p))
        assert an[p]==ap
        small_trace_checks.append({'prime':p,'a_p':ap})
    records=[]
    for row in data['rows']:
        d=row['twist'];K=row['terms'];assert d in spec['twists']
        if d==1:characters=[1]
        else:
            assert ZZ(d).is_prime() and d%4==1
            characters=[0]+[1 if pow(x,(d-1)//2,d)==1 else -1 for x in range(1,d)]
        twisted=[an[n]*characters[n%d] for n in range(K+1)] if d!=1 else an[:K+1]
        assert r.digest(json.dumps(twisted,separators=(',',':')).encode())==row['coefficient_sha256']
        E=EllipticCurve(QQ,row['model']);assert list(map(int,E.minimal_model().ainvs()))==row['minimal_model']
        local=[]
        for saved in row['local_data']:
            place=saved['prime'];datum=E.local_data(place,algorithm='generic')
            assert int(datum.conductor_valuation())==saved['conductor_exponent']
            assert str(datum.kodaira_symbol())==saved['kodaira']
            assert int(datum.tamagawa_number())==saved['tamagawa_number']
            local.append({'prime':place,'conductor_exponent':int(datum.conductor_valuation()),'kodaira':str(datum.kodaira_symbol())})
        assert E.root_number()==row['root_number']==-1
        # Independent implementation: vector E1 evaluation in MPFR, with its
        # documented explicit truncation and floating-point error budget.
        value,error=E.lseries().deriv_at1(K,prec=spec['limits']['precision_bits'])
        lo=value.exact_rational()-error.exact_rational();hi=value.exact_rational()+error.exact_rational()
        arb=row['L_derivative_interval'];alo,ahi=QQ(arb['lower']),QQ(arb['upper'])
        assert max(lo,alo)<=min(hi,ahi)
        independent_nonzero=bool(lo>0 or hi<0)
        assert independent_nonzero==row['nonzero_derivative_certified']
        ct=prior['original_full_CT_matrix'] if d==1 else next(x['full_CT_matrix'] for x in prior['joint_prime_table'] if x['prime']==d)
        ct_rank=r.rank(map(r.pack,ct));radical=[v for v in range(8) if all((v&r.pack(line)).bit_count()%2==0 for line in ct)]
        record={'twist':d,'all_coefficient_twist_identities_verified':K,'generic_Tate_checks':local,
                'independent_derivative_interval':{'lower':str(lo),'upper':str(hi),'value':str(value),'error':str(error)},
                'overlaps_Arb_enclosure':True,'nonzero_derivative_certified':independent_nonzero,
                'full_CT_rank':ct_rank,'full_CT_radical_masks':radical}
        if independent_nonzero:
            assert ct_rank==2 and len(radical)==2
            record.update({'exact_MW_rank':1,'Sha_finite':True,'rational_Kummer_line_nonzero_mask':radical[1],
                           'Sha_2_primary_structure':'C2 x C2','Sha_4_equals_Sha_2':True})
        else:
            record.update({'exact_MW_rank':'UNKNOWN','Sha_finite':'UNKNOWN',
                           'rational_Kummer_subspace':'UNKNOWN','Sha_2_primary_structure':'UNKNOWN'})
        records.append(record);print('verified',d,independent_nonzero,'CT',ct_rank,flush=True)
    assert next(x for x in records if x['twist']==97)['rational_Kummer_line_nonzero_mask']==5
    independent_source=inspect.getsourcefile(Lseries_ell)
    return {'schema':'rank-jump.zero-ct-analytic-gate-verification.v1','status':'PASS',
            'bindings':{str(p.relative_to(r.ROOT)):r.digest(p.read_bytes()) for p in
                        (Path(__file__),source.OUTPUT,source.PRIOR,source.PROTOCOL,Path(r.__file__))},
            'independent_lseries_source':{'module':'sage.schemes.elliptic_curves.lseries_ell','sha256':r.digest(Path(independent_source).read_bytes())},
            'base_small_prime_character_checks':small_trace_checks,'records':records,
            'external_rank_theorem':{'source':'https://annals.math.princeton.edu/wp-content/uploads/annals-v191-n2-p01-s.pdf',
                                    'location':'Introduction, page330, recalled Gross-Zagier/Kolyvagin implication; not the converse Theorem A'},
            'boundary':'Derivative nonvanishing gives exact rank one and finite Sha. Zero-containing intervals at 41 and113 remain inconclusive. The rational class at97 is certified without a point coordinate.'}


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('mode',choices=['build','check']);args=p.parse_args();result=compute()
    if args.mode=='build':r.write_new(OUTPUT,result)
    else:assert r.read(OUTPUT)==result
    print('PASS independent analytic, coefficient, local and rational-line verification')
