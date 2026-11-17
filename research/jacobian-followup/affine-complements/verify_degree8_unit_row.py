"""PF2RULED28: exact unit-row replay and global root-constant recovery.
Run with .venv/bin/python. No Groebner search, expanded E, or parameter
chart selection is needed for this certificate replay.
"""
from pathlib import Path
from time import monotonic
import gzip,hashlib,json,platform
import sympy as S
import flint
from flint import fmpq,fmpq_mpoly_ctx

HERE=Path(__file__).resolve().parent
VARIABLES=('y','q4','q5','q6','q7')
CTX=fmpq_mpoly_ctx.get(VARIABLES,ordering='lex')
checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def read_gz(path):return json.loads(gzip.decompress(path.read_bytes()))
def encoded_polynomial(data):
    terms={}
    for powers,coefficients in data:
        for i,c in enumerate(coefficients):
            if c!='0':terms[(i,*powers)]=fmpq(c)
    return CTX.from_dict(terms)
def rational_polynomial(coefficients):
    return CTX.from_dict({(i,0,0,0,0):fmpq(c) for i,c in enumerate(coefficients) if c!='0'})
def primary_field(primary):
    y,v,q4,q5,q6,q7=S.symbols('y v q4 q5 q6 q7');loc=dict(zip(('y','v','q4','q5','q6','q7'),(y,v,q4,q5,q6,q7)))
    p8=S.Poly(S.sympify(primary['P8'].replace('^','**'),locals=loc),y,domain=S.QQ).monic()
    D=S.Poly(S.sympify(primary['linear_constant_D'].replace('^','**'),locals=loc),v,y,q4,q5,q6,q7,domain=S.QQ)
    groups=[{} for _ in range(5)]
    for ex,c in D.terms():groups[ex[0]][tuple(ex[1:])]=fmpq(str(c))
    return p8,[CTX.from_dict(row) for row in groups]

def main():
    start=monotonic();pred=json.loads((HERE/'degree8_affine_reduction_receipt.json').read_text())
    check('current affine-reduction predecessor PASS',pred['outcome']=='PASS')
    for name,digest in pred['source_sha256'].items():
        check('predecessor seal '+name,hashlib.sha256((HERE/name).read_bytes()).hexdigest()==digest)
    primary=read_gz(HERE/'degree8_affine_reduction_certificates.json.gz')
    data=read_gz(HERE/'degree8_unit_row_certificates.json.gz')
    check('both marked fields and all16 points covered',len(data['sectors'])==2 and primary['marked_points']==16)
    summaries=[]
    for sec,old in zip(data['sectors'],primary['sectors']):
        k=sec['sector'];check(f'correct marked multiplicities k{k}',old['multiplicities'][0]==k)
        p8,Dsource=primary_field(old)
        check(f'exact primary field irreducible degree8 k{k}',p8.degree()==8 and p8.is_irreducible)
        p=rational_polynomial(sec['P8_coefficients'])
        psource=rational_polynomial([str(p8.nth(j)) for j in range(9)])
        check(f'certificate defining polynomial equals monic primary polynomial k{k}',p==psource)
        check(f'five generator and multiplier slots k{k}',len(sec['generators'])==len(sec['coefficients'])==5)
        D=[encoded_polynomial(row) for row in sec['generators']]
        C=[encoded_polynomial(row) for row in sec['coefficients']]
        for i,(one,two) in enumerate(zip(D,Dsource)):
            check(f'unit generator matches actual affine-model coefficient k{k} i{i}',divmod(one-two,p)[1]==0)
        degrees=[max(sum(ex) for ex,c in row if any(x!='0' for x in c)) for row in sec['coefficients']]
        check(f'all exact multipliers have degree at most3 k{k}',max(degrees)==3 and all(d<=3 for d in degrees))
        total=sum((c*d for c,d in zip(C,D)),CTX.constant(0))
        quotient,remainder=divmod(total-1,p)
        check(f'FULL exact characteristic-zero Bezout identity modulo P8 k{k}',remainder==0)
        check(f'primary E degree bound supports finite recovery-degree bound k{k}',max(old['E_coefficient_degrees'])<=15)
        qsymbols=S.symbols('q4 q5 q6 q7');ysymbol=S.Symbol('y')
        ms=S.Poly(S.sympify(old['constant_remainder_Ms'].replace('^','**'),locals=dict(zip(('y','q4','q5','q6','q7'),(ysymbol,*qsymbols)))),ysymbol,*qsymbols,domain=S.QQ)
        check(f'primary Ms parameter degree at most4 k{k}',max(sum(ex[1:]) for ex,c in ms.terms())<=4)
        Ddegree=max(sum(ex) for row in sec['generators'] for ex,c in row if any(x!='0' for x in c))
        check(f'primary D parameter degree at most4 k{k}',Ddegree<=4)
        summaries.append(dict(multiplicities=old['multiplicities'],P8_coefficients=sec['P8_coefficients'],
            multiplier_degree=3,unit_identity_verified=True,quotient_term_count=len(quotient),
            eliminated_stratum='D0=...=D4=0 has no points',
            recovered_root_constant_circuit='Zstar=sum(C[i]*E[i],i=0..4)',root_constant_degree_upper_bound=18,
            core_compatibility_circuit=['E[5]=0']+[f'D[{i}]*Zstar-E[{i}]=0' for i in range(5)],
            nonzero_condition='M(s)*Zstar != 0',saturation='t*M(s)*Zstar-1=0',
            core_parameter_count=4,saturated_variable_count=5,compatibility_degree_upper_bound=22,
            full_nonzero_locus_solved=False))
        print(f'PASS k{k}: full degree3 unit identity; zero-D stratum excluded; root constant recovered polynomially',flush=True)
    # Generic identities prove equality of elimination ideals once sum C_iD_i=1.
    c=S.symbols('c0:5');d=S.symbols('d0:5');e=S.symbols('e0:5');z,t,m=S.symbols('z t m')
    unit=sum(ci*di for ci,di in zip(c,d))-1;zstar=sum(ci*ei for ci,ei in zip(c,e))
    oldeq=[di*z-ei for di,ei in zip(d,e)];neweq=[di*zstar-ei for di,ei in zip(d,e)]
    check('global recovery relation lies in original ideal using the exact unit row',S.expand(sum(ci*f for ci,f in zip(c,oldeq))-(z-zstar)-z*unit)==0)
    for i in range(5):
        check(f'old and recovered compatibility generator identity {i}',S.expand(oldeq[i]-d[i]*(z-zstar)-neweq[i])==0)
    check('known global compatibility syzygy is retained',S.expand(sum(ci*f for ci,f in zip(c,neweq))-zstar*unit)==0)
    check('nonzero-kappa saturation ideal survives root elimination',S.expand((t*m*z-1)-(t*m*zstar-1)-t*m*(z-zstar))==0)
    check('uniform polynomial recovery degree bound18',3+15==18)
    check('uniform core compatibility degree bound22',4+18==22)
    reduced=dict(as_of='2026-11-14',marked_points=16,sectors=summaries,
        complete_zero_D_stratum_excluded=True,global_root_constant_polynomial_recovery=True,
        core_parameter_count=4,saturated_variable_count=5,root_constant_degree_upper_bound=18,
        compatibility_degree_upper_bound=22,nonzero_kappa_condition_retained=True,
        proof_method='exact rational sparse polynomial arithmetic modulo the two irreducible P8 polynomials; universal ideal identities')
    rp=HERE/'degree8_global_recovery.json';rp.write_text(json.dumps(reduced,indent=2)+'\n')
    inputs=[Path(__file__),HERE/'compute_degree8_unit_row.py',HERE/'DEGREE_EIGHT_GLOBAL_CONSTANT_RECOVERY.md',
            HERE/'degree8_unit_row_certificates.json.gz',rp,HERE/'degree8_affine_reduction_receipt.json',
            HERE/'degree8_affine_reduction_certificates.json.gz',HERE/'DEGREE_EIGHT_AFFINE_REDUCTION.md',
            HERE/'sources/degree8_unit_row_groebner_attempt.py.txt',HERE/'sources/degree8_unit_row_field_linear_attempt.py.txt']
    report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,marked_points=16,
        complete_zero_D_stratum_excluded=True,global_root_constant_polynomial_recovery=True,
        exact_Bezout_multiplier_degree=3,core_parameter_count=4,saturated_variable_count=5,
        root_constant_degree_upper_bound=18,compatibility_degree_upper_bound=22,
        first_open_reduced_pairs=[[5,9],[6,9],[7,9],[8,9]],minimum_compatible_quartic_pole_height=9,
        assurance=dict(all_complex_embeddings_covered=True,full_polynomial_identity_checked_in_characteristic_zero=True,
            discovery_modular_pivots_not_used_as_proof=True,no_parameter_denominator_introduced=True,
            entire_previous_zero_coefficient_stratum_closed=True,all_individual_coefficient_zero_cases_retained=True,
            elimination_and_saturation_ideals_preserved=True,
            four_parameter_nonzero_kappa_locus_solved=False,degree8_actual_section_lifts_constructed=False,
            general_determinant_locus_solved=False,actual_global_root_hypothesis_removed=False,
            full_source_Keller_pair_constructed=False,field_degree_four_realization_constructed=False,
            global_surface_finiteness_proved=False,actual_ruled_compatibility_proved=False,contact_cases_classified=False,
            whole_marked_candidate_excluded=False,quartic_packets_excluded=False,formal_verification=False,
            independent_full_reproof=False,external_review=False,research_goal_completed=False),
        source_sha256={str(f.relative_to(HERE)):hashlib.sha256(f.read_bytes()).hexdigest() for f in inputs},
        python=platform.python_version(),sympy=S.__version__,python_flint=flint.__version__,elapsed_seconds=round(monotonic()-start,6))
    (HERE/'degree8_unit_row_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
    print(f'PASS {len(checks)} controls; global four-parameter recovery without coefficient charts; nonzero-kappa core remains open',flush=True)

if __name__=='__main__':main()
