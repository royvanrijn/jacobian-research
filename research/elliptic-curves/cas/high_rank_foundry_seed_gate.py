"""Bounded seed admissibility is an input outcome, never a rank upper bound."""
from fractions import Fraction as F

import compact_atlas_specialization as atlas
from future_point_admission import FinitePointAdmission
from mod2_reduction_independence import short_curve_has_no_rational_2_torsion_modular_certificate
from v3_warm_support import read, require, sha


def assess(row):
    family=next(f for f in read(atlas.ATLAS)['families'] if f['family']==row['family'])
    model,points=atlas.specialize(family,row['parameter'])
    require(model==tuple(map(F,row['model'])),'selected equation differs')
    # Empty initialization lets us inspect a bounded admission miss explicitly.
    # Invalid points, specialization mismatches and backend exceptions still fail.
    admission=FinitePointAdmission(model,(),prime_bound=1000)
    results=[admission.consider(point) for point in points]
    full=all(r['status']=='INDEPENDENT_FINITE_COLUMN' for r in results)
    torsion=next((p for p in admission.primes if p<=200 and
                  short_curve_has_no_rational_2_torsion_modular_certificate(model,p)),None)
    return {'status':'PASS_GENERIC_SEED_GATE' if full and torsion else 'UNRESOLVED_GENERIC_SEED',
            'id':row['id'],'family':row['family'],'parameter':row['parameter'],
            'curve':list(map(str,model)),'specialized_points':[list(map(str,p)) for p in points],
            'required_points':17,'admitted_columns':len(admission.points),
            'admissions':results,'prime_bound':1000,'primes':admission.primes,
            'no_rational_2_torsion_prime':torsion,'torsion_prime_bound':200,
            'point_searches':0,'atlas_sha256':sha(atlas.ATLAS),
            'claim_boundary':'Bounded finite-column seed test only. An incomplete witness is UNKNOWN; no elliptic rank upper bound or dependence conclusion.'}
