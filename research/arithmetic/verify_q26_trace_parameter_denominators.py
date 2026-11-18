#!/usr/bin/env python3
"""Exact leading-row identities, joint root-free certificate and F43 census."""
from __future__ import annotations
import argparse,hashlib,json,platform
from pathlib import Path
import flint
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]

def require(test,message):
    if not test:raise AssertionError(message)

def evaluate(a,x,p=None):
    value=0
    for c in reversed(a):
        value=value*x+int(c)
        if p:value%=p
    return value

def factor_product(factors):
    value=1
    for p,e in factors.items():value*=int(p)**e
    return value

def prime43_poly_certificate(name,f):
    g=flint.nmod_poly(list(f),43)
    require(g.degree()==f.degree(),'Degree drop modulo43')
    require(g.gcd(g.derivative()).degree()==0,'Nonseparable polynomial modulo43')
    x=flint.nmod_poly([0,1],43)
    require(g.gcd(pow(x,43,g)-x).degree()==0,'A residue root modulo43')
    values=[evaluate(list(f),a,43) for a in range(43)]
    require(all(values),'Root-free evaluation certificate')
    return {'polynomial':name,'degree':g.degree(),
            'leading_residue':int(g[g.degree()]),'values_on_F43':values,
            'gcd_with_derivative_degree':0,'gcd_with_X43_minus_X_degree':0}

def affine_census(coeff):
    rows=[]
    for t in range(43):
        values=[evaluate(a,t,43) for a in coeff]
        ds=[i*values[i]%43 for i in range(1,len(values))]
        dt=[evaluate([j*a[j] for j in range(1,len(a))],t,43) for a in coeff]
        for s in range(43):
            if evaluate(values,s,43):continue
            fs=evaluate(ds,s,43);ft=evaluate(dt,s,43)
            rows.append({'s':s,'t':t,'F_s':fs,'F_t':ft,'smooth':bool(fs or ft)})
    return rows

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-receipt',type=Path,
                        default=HERE/'q26_trace_parameter_denominator_verification.json')
    args=parser.parse_args()
    data=json.loads((HERE/'q26_trace_parameter_denominator_inputs.json').read_text())
    curve=json.loads((HERE/'q26_high_genus_affine_trace_curve_inputs.json').read_text())['curve']
    coeff=curve['coefficients_ascending_s']
    require([len(a)-1 for a in coeff]==list(range(144,131,-1)),
            'Newton row-degree hypothesis')
    A=flint.fmpz_poly([a[-1] for a in coeff])
    lc=flint.fmpz_poly(coeff[-1]);row=data['leading_s_factorization']
    B=flint.fmpz_poly(row['B56_coefficients_ascending'])
    G=flint.fmpz_poly(row['G2_coefficients_ascending'])
    require(B.degree()==56 and G==flint.fmpz_poly([-1,-1,1]),'Factor degree/identity')
    require(lc==25*B**2*G**10,'Entire leading coefficient factor product')
    require(int(A[0])==factor_product(data['leading_constant_factorizations']['A0']),
            'A0 integer factorization')
    require(int(A[12])==factor_product(data['leading_constant_factorizations']['Ad']),
            'Ad integer factorization')
    primes=sorted(set(map(int,data['leading_constant_factorizations']['A0']))|
                  set(map(int,data['leading_constant_factorizations']['Ad'])))
    require(primes==data['excluded_leading_primes'],'Excluded leading primes')
    certs=[prime43_poly_certificate(name,f) for name,f in [('A12',A),('B56',B),('G2',G)]]
    require([a['leading_residue'] for a in certs]==[36,16,1],'Leading-unit residues')
    require(int(A[0])%43==22,'A0 unit residue')
    points=affine_census(coeff)
    require(points==data['prime43']['affine_points'],'Complete affine census input')
    require(len(points)==52 and sum(p['smooth'] for p in points)==47,'Affine census counts')
    paths=[Path(__file__).resolve(),HERE/'Q26_UNIFORM_TRACE_PARAMETER_DENOMINATORS.md',
           HERE/'q26_trace_parameter_denominator_inputs.json',
           HERE/'q26_high_genus_affine_trace_curve_inputs.json',
           HERE/'Q26_HIGH_GENUS_AFFINE_TRACE_CURVE.md',
           HERE/'Q26_SMALL_PRIME_LOCAL_BRANCHES.md',HERE/'ODD_AFFINE_CRITICAL_TRACE_CURVES.md',
           HERE/'sources/milne-ant300.pdf']
    result={'schema_version':1,'result':'PASS','python':platform.python_version(),
            'python_flint':flint.__version__,
            'Newton_rows':{'d':12,'D':144,'degrees':[len(a)-1 for a in coeff]},
            'leading_s_factorization':'25*B56(t)^2*(t^2-t-1)^10',
            'excluded_leading_primes':primes,'prime43_polynomial_certificates':certs,
            'prime43_affine_census':{'all_1849_pairs_checked':True,'total':52,'smooth':47,
                                     'first_smooth_point':next(p for p in points if p['smooth'])},
            'uniform_conclusion':'At every prime with leading units and joint root-free A/LC reductions, all affine Qp accessory parameters s,t are integral. A positive-natural-density set of such primes exists by the joint Frobenius43 class and Chebotarev.',
            'boundary_conclusion':'The smooth projective normalization has no Q43 point at a pole of s or t, hence no Q-rational point outside this affine coordinate chart.',
            'scope':'Trace parameters only. Not finite-S integrality, not integrality of reconstructed q/N coefficients, not rational nonexistence, not a low-genus quotient or all-prime local classification.',
            'source_sha256':{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in paths},
            'formal_verification':False,'independent_external_replay':False,'external_human_review':False}
    args.write_receipt.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: exact Q26 Newton denominator criterion, joint43 certificate and52/47 census.')

if __name__=='__main__':main()
