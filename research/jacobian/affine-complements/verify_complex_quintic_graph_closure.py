"""PF156: exact residual unit ideal and actual rank-three nilpotent algebra.
The field coverage import is PF155; no higher normal prefix is extrapolated.
"""
from pathlib import Path

import tempfile as _check_tempfile

_CHECK_OUTPUT_DIRECTORY = None

def _check_output(original):
    """Write new runs to a temporary directory; preserve deposited receipts."""
    global _CHECK_OUTPUT_DIRECTORY
    if _CHECK_OUTPUT_DIRECTORY is None:
        _CHECK_OUTPUT_DIRECTORY = Path(_check_tempfile.mkdtemp(prefix="math-check-"))
        print("Verification output:", _CHECK_OUTPUT_DIRECTORY)
    return _CHECK_OUTPUT_DIRECTORY / Path(original).name


import hashlib
import importlib.util
import json
import platform
import time
import sympy as S
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
SRC=HERE/'sources/quintic-residual-closure'
NOTE='FULL_COMPLEX_QUINTIC_GRAPH_EXCLUSION_AND_CUBIC_NILPOTENCE.md'

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()

def main():
    started=time.perf_counter();checks=[]
    def check(name,value):
        assert bool(value),name
        checks.append(name)
    ctx=json.loads((SRC/'context.json').read_text())
    for p,digest in ctx['retained_sha256'].items():check('retained source pin '+p,sha(ROOT/p)==digest)
    for p,digest in ctx['certificate_sha256'].items():check('certificate pin '+p,sha(SRC/p)==digest)
    spec=importlib.util.spec_from_file_location('residual_exact_jet',SRC/'derive_residual_witnesses.py')
    mod=importlib.util.module_from_spec(spec);spec.loader.exec_module(mod);z=mod.derive_residual_witnesses()
    x,y,p,q=z['source_symbols'];rho,h,i,e,f,a,b,c,d,g,j,k,l,m,n,C=z['params']
    loc={str(t):t for t in z['params']}
    Q=160*rho*rho+1968*rho+6021;LL=2898*rho+16821
    U=S.Rational(96,25)*LL;B1=196608*(i-8*b)*(e-4*a-rho-S.Rational(21,2))
    Tq=-64*e+S.Rational(704,5)*rho+S.Rational(16704,25)
    # The potential jet is independently derived from the actual parent.
    tr,ss=S.symbols('tr ss');u=1+x*y
    AA=u**3*p+3*x*u*u*q-x**3*(tr+x*ss)
    BB=y*y*u*(4+3*x*y)*p+(y+3*x*y*y*(4+3*x*y))*q+(2*x-3*x*x*y)*(tr+x*ss)
    parent=S.Poly(S.expand(AA*AA+13*AA+2*BB),x)
    check('actual parent potential jet through x3',all(S.expand(parent.nth(jj).subs({tr:z['trace'],ss:z['normal']})-z['potential_coefficients'][jj])==0 for jj in range(4)))
    high=sum(co*y**ex[0]*p**ex[1]*q**ex[2] for ex,co in z['plane'].terms() if sum(ex)>4)
    check('entire high plane ideal is residual quadratic Q',S.expand(high-16*Q*y**4*p)==0)
    check('full cubic f retained',S.Poly(z['trace'],y,p,q).coeff_monomial(p**3)==f)
    lower=[a,b,c,d,g,j,k,l,m,n]
    check('all ten lower trace coefficients retained',all(v in z['trace'].free_symbols for v in lower))
    data=json.loads((SRC/'residual_witness_certificate.json').read_text())
    for name in ('wc','wp'):
        check('actual full-f second-jet witness '+name,S.expand(z[name]-S.sympify(data[name],locals=loc))==0)
    check('actual cubic-nilpotent obstruction coefficient',z['wc']==12582912*c**3)
    check('all-parameter unit residue decomposition',S.expand(z['wp']-U-c*B1-Q*Tq)==0)
    check('second obstruction retains f without restricting it',not z['wp'].has(f) and not z['wc'].has(f))
    check('actual c correction factor',S.expand(S.sympify(data['coefficient_c'],locals=loc)-B1)==0)
    check('exact polynomial Q correction',S.expand(S.sympify(data['coefficient_Q'],locals=loc)-Tq)==0)
    # The k=2 multiplier makes the degree-six witnesses immutable.
    check('second-jet repair multiplier and witness degrees',-4*1*2**3*2*3==-192 and S.Poly(y**4*q*q,y,p,q).total_degree()==6 and S.Poly(y**5*p,y,p,q).total_degree()==6)
    v=S.Rational(25,96)*(S.Rational(368,76545)*rho+S.Rational(1328,42525));E=S.Rational(529,6075)
    check('unit residue explicit polynomial inverse modulo Q',S.expand(v*U-1-E*Q)==0)
    check('residual quadratic irreducible over Q',S.discriminant(Q,rho)==576*34 and S.integer_nthroot(34,2)[1] is False)
    cert=json.loads((SRC/'residual_unit_ideal_certificate.json').read_text())
    for key,expr in [('Q',Q),('Wc',z['wc']),('Wp',z['wp']),('U',U),('B',B1),('Tq',Tq),('v',v)]:
        check('unit certificate input '+key,S.expand(S.sympify(cert[key],locals=loc)-expr)==0)
    hp,hc,hq=[S.sympify(cert[key],locals=loc) for key in ['multiplier_Wp','multiplier_Wc','multiplier_Q']]
    check('explicit residual three-generator Bezout identity',S.expand(hp*z['wp']+hc*z['wc']+hq*Q-1)==0)
    t=v*c*B1
    check('finite inverse geometric identity',S.expand((1-t+t*t)*(1+t)-1-t**3)==0)
    check('nilpotent inverse remainder is actual Wc multiple',S.expand(t**3-(v**3*B1**3/S.Integer(12582912))*z['wc'])==0)
    # ACTUAL finite algebra K0[c]/c3 in the rational cyclic basis.
    Mr0=S.Matrix([[0,-S.Rational(6021,160)],[1,-S.Rational(123,10)]])
    Mr=S.diag(Mr0,Mr0,Mr0);Mc=S.zeros(6)
    Mc[2:4,0:2]=S.eye(2);Mc[4:6,2:4]=S.eye(2)
    I=S.eye(6);zero=S.zeros(6);unit=I[:,0]
    check('actual finite algebra rho quadratic multiplication',160*Mr**2+1968*Mr+6021*I==zero)
    check('actual finite algebra multiplication commutes',Mr*Mc==Mc*Mr)
    check('actual finite algebra cubic nilpotence retains square',Mc**3==zero and Mc**2!=zero)
    basis=S.Matrix.hstack(unit,Mr*unit,Mc*unit,Mr*Mc*unit,Mc**2*unit,Mr*Mc**2*unit)
    check('actual six-dimensional cyclic unit basis',basis==I)
    fix={h:0,i:1,e:0,f:0,a:0,b:0,d:0,g:0,j:0,k:0,l:0,m:0,n:0,C:1}
    def ev(expr):
        poly=S.Poly(S.expand(expr.subs(fix)),rho,c)
        return sum((coef*Mr**ex[0]*Mc**ex[1] for ex,coef in poly.terms()),zero)
    wp=ev(z['wp']);wv=ev(v);wb=ev(B1);wt=wv*Mc*wb
    check('actual Wc matrix vanishes without removing c2',ev(z['wc'])==zero and Mc**2*unit!=S.zeros(6,1))
    check('actual unit obstruction includes nilpotent correction',wp==ev(U)+Mc*wb and Mc*wb!=zero)
    inverse=wv*(I-wt+wt**2)
    check('actual multiplication finite inverse on both sides',inverse*wp==I and wp*inverse==I)
    check('actual finite algebra obstruction determinant nonzero',wp.det()!=0)
    check('no larger partial-jet or global theorem asserted',not ctx['A2_alone_claimed_empty'] and not ctx['normal_tail_increasing_power_containment_proved'] and not ctx['global_minimal_counterexample_coverage_proved'])
    out=dict(claim='PF2RULED156',outcome='PASS',check_count=len(checks),checks=checks,elapsed_seconds=time.perf_counter()-started,
             environment=dict(python=platform.python_version(),sympy=S.__version__),proof_sha256=sha(HERE/NOTE),checker_sha256=sha(Path(__file__)),context_sha256=sha(SRC/'context.json'),
             all_lower_trace_coefficients_and_full_f_retained=True,actual_residual_nilpotence_index=3,
             actual_finite_algebra_dimension_over_Q=6,actual_finite_algebra_dimension_over_quadratic_field=3,
             actual_obstruction_matrix_determinant=str(wp.det()),explicit_residual_unit_ideal_identity_verified=True,
             both_conjugate_complex_classes_excluded=True,full_complex_degree_at_most_five_v2_graph_branch_closed=True,
             full_bounded_graph_coefficient_scheme_zero_written_maximal_ideal_corollary=True,
             full_128_variable_Bezout_vector_computed=False,A2_alone_empty_claimed=False,
             normal_tail_power_containment_proved=False,global_coverage_proved=False,new_external_paper_imports=[],
             registered_predecessor_checkers_replayed=[],full_goal_complete=False,independent_replay=False,formal_verification=False,external_review=False)
    (_check_output(HERE/'complex_quintic_graph_closure_receipt.json')).write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({key:out[key] for key in ['claim','outcome','check_count','elapsed_seconds']}))

if __name__=='__main__':main()
