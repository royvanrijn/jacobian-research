"""Internal exact replay: whole H3 certificates at fixed indices, plus uniform initial pivots."""
from pathlib import Path
import argparse, copy, hashlib, json, platform, time
import sympy as s
import gap_three_hurwitz as h
from gap_three_wronskian import centered_coefficients
from higher_gap_adjoint_rank import auxiliary, gap_three_W

HERE = Path(__file__).resolve().parent
SOURCE = 'sources/holtz-khrushchev-kushel-2015-v1.pdf'
SOURCE_HASH = '05aaeaeddc525d095461d161bab6417415c36e22537f9322e5760a69c7289db5'
counts = {}
def check(label, value):
    assert value, label
    counts[label] = counts.get(label, 0)+1

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--max-index', type=int, default=16)
    args = parser.parse_args()
    check('index_bound_at_least_two', args.max_index >= 2)
    check('HKK_primary_source_pin', hashlib.sha256((HERE/SOURCE).read_bytes()).hexdigest() == SOURCE_HASH)
    beta, z, x = s.symbols('beta z x')
    n, b, C, D = h.uniform_pivot_numerators()
    n, b, C = h.uniform_top_coefficients(11)
    for j in range(1, len(C)):
        rhs = (4*n+b+2)*C[j-1].subs(b,b+1)
        if j>=2:
            rhs += (4*n+b+2)*(2*n+b-1+j)*C[j-2].subs(b,b+1)
        check('uniform_top_coefficient_difference_identity',
              s.expand(C[j].subs(n,n+1)-C[j]-rhs) == 0)
        check('uniform_top_coefficient_initial_value', C[j].subs(n,0) == 0)
    uniform = []
    for j, polynomial in D.items():
        shift = (j+1)//2
        P = s.Poly(s.expand(polynomial.subs(n,x+shift)),x,b)
        check('uniform_pivot_numerator_positive_coefficients', all(v>0 for v in P.coeffs()))
        check('uniform_pivot_numerator_positive_constant', P.TC()>0)
        uniform.append({'pivot':j, 'n_lower_bound':shift, 'variables':['x','beta'],
                        'substitution':'n=x+'+str(shift),
                        'terms':[[list(monomial),str(value)] for monomial,value in P.terms()]})
    delta = {j:h.special_minor(j,C) for j in range(1,9)}
    check('special_minors_match_initial_pivot_numerators',
          s.expand(delta[3]-D[3])==0 and s.expand(delta[4]-D[4])==0 and
          s.expand(delta[5]-D[5])==0 and s.expand(C[1]*delta[6]-D[6])==0)
    for j in range(1,9):
        shift=(j+1)//2
        P=s.Poly(s.expand(delta[j].subs(n,x+shift)),x,b)
        check('uniform_later_special_minor_positive_coefficients',all(v>0 for v in P.coeffs()))
        check('uniform_later_special_minor_positive_constant',P.TC()>0)
        uniform.append({'special_minor':j,'n_lower_bound':shift,'variables':['x','beta'],
                        'substitution':'n=x+'+str(shift),
                        'terms':[[list(monomial),str(value)] for monomial,value in P.terms()]})
    G4=C[2]*C[7]-C[1]*C[8]
    G5=C[8]*D[3]-C[2]*(C[1]*C[9]-C[10])
    G6=D[4]*(C[1]*C[9]-C[10])-D[3]*(C[2]*C[10]-C[1]*C[11])
    check('uniform_h7_special_minor_elimination_identity',
          s.expand(G4*D[5]-D[4]*G5-C[2]*delta[7])==0)
    check('uniform_h8_special_minor_elimination_identity',
          s.expand(G5*D[6]-D[5]*G6-C[1]*D[3]*delta[8])==0)
    formula = {0:s.Integer(1), 1:C[1], 2:C[2], 3:D[3]/C[1],
               4:D[4]/C[2], 5:D[5]/D[3], 6:D[6]/(C[1]*D[4]),7:delta[7]/delta[5],8:delta[8]/delta[6]}
    profiles = []
    for index in range(2,args.max_index+1):
        started = time.monotonic()
        c = centered_coefficients(index,beta)
        rows = h.euclidean_rows(index,beta)
        check('complete_generalized_Euclidean_identity_replay', h.check_euclidean_rows(rows,c,index))
        leading = [row[2*index-i] for i,row in enumerate(rows)]
        signs = [h.positive_fraction_certificate(v,beta) for v in leading]
        check('all_generalized_Euclidean_leading_coefficients_positive',
              len(leading)==2*index+1 and all(v['positive_on_nonnegative_axis'] for v in signs))
        check('constant_endpoint_product',
              s.expand(leading[-1]-s.prod((index+beta+3*j)*(index+beta+3*j+1)
                                       for j in range(index))) == 0)
        if index<=5:
            for j in range(min(8,2*index)+1):
                check('uniform_initial_pivot_formula_matches_full_algorithm',
                      s.cancel(leading[j]-formula[j].subs({n:index,b:beta})) == 0)
            B=auxiliary(index,index+beta-1,3,z)
            W=gap_three_W(B.as_expr(),z)
            check('independent_auxiliary_Wronskian_input',
                  s.expand(W.as_expr()-sum(c[k]*z**k for k in range(2*index+1))) == 0)
        if index==3:
            tampered = copy.deepcopy(rows)
            tampered[4][2] += 1
            check('tampered_Euclidean_row_rejected', not h.check_euclidean_rows(tampered,c,index))
        profiles.append({'n':index,'beta_scope':'every real beta>=0',
                         'all_infinite_H3_minors_nonnegative_via_HKK_16':True,
                         'leading_coefficients':signs})
        print(json.dumps({'n':index,'status':'pass','seconds':round(time.monotonic()-started,3)}),flush=True)
    check('negative_fraction_control_rejected',
          not h.positive_fraction_certificate(1-beta,beta)['positive_on_nonnegative_axis'])
    check('negative_denominator_control_rejected',
          not h.positive_fraction_certificate(1/(1-beta),beta)['positive_on_nonnegative_axis'])
    certificates = {'arithmetic':'exact Q[n,beta] and Q(beta)',
                    'uniform_initial_pivot_numerators':uniform,
                    'fixed_index_profiles':profiles,
                    'fixed_index_input_generation':'centered_coefficients and full HKK algorithm; no sampled beta values',
                    'all_index_TN_proved':False}
    certificate_name = 'gap_three_hurwitz_certificates.json'
    (HERE/certificate_name).write_text(json.dumps(certificates,indent=2)+'\n')
    names=['gap_three_hurwitz.tex','gap_three_hurwitz.py','gap_three_hurwitz_verify.py',
           certificate_name,'gap_three_wronskian.py','higher_gap_adjoint_rank.py',SOURCE]
    receipt = {'status':'pass','canonical_proof':'gap_three_hurwitz.tex',
               'command':'python research/factorial/image-rigidity/gap_three_hurwitz_verify.py --max-index '+str(args.max_index),
               'python':platform.python_version(),'sympy':s.__version__,
               'HKK_external_dependency':{'url':'https://arxiv.org/pdf/1506.07379v1',
                                        'theorems':[4,16,31],'primary_statements_and_algorithm_inspected':True,
                                        'independent_full_paper_review':False},
               'whole_infinite_H3_TN_fixed_indices':list(range(2,args.max_index+1)),
               'fixed_index_parameter_scope':'every real beta>=0',
               'uniform_initial_pivots':'h0,...,h8 positive for n>=4,beta>=0; h0,...,h6 for n=3; h0,...,h4 for n=2',
               'all_index_gap_three_SFC_closed':False,'unrestricted_binomial_SFC_closed':False,
               'finite_index_extrapolation':False,'independent_replay':False,
               'formal_verification':False,'external_review':False,'novelty_claim':False,
               'total':sum(counts.values()),'counts':counts,
               'files':{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in names}}
    (HERE/'gap_three_hurwitz_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps({k:v for k,v in receipt.items() if k!='files'},indent=2))
if __name__=='__main__':main()
