"""Exact producer controls; certificate truth is proved separately, not fitted."""
from pathlib import Path
import copy
import hashlib
import json
import subprocess
import sys
import sympy as S
import bernstein_slope_producer as B

p=Path(__file__).resolve().parent
proof=p/'EFFECTIVE_BERNSTEIN_SLOPE_SELECTION_AND_CUTOFFS.md'
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
X,Y=B.X,B.Y

def summary(polynomial,J):
    return dict(polynomial_terms=B.terms(S.Poly(polynomial,X,Y,domain=S.QQ)),validity_start=1,step_bound=J)

def run(name,certificate,phase,symbol,expected,pure=True):
    plan=B.make_plan(certificate)
    n=len(phase[0][0]);D=max(sum(e) for e,c in phase)
    jet=max(D*plan['symbol_jet_order_required_per_phase_degree'],max((sum(e) for e,c in symbol),default=0))
    evaluated=B.evaluate_rational(plan,phase,symbol,jet)
    out=B.resolve(plan,evaluated['evaluations'],all_positive_pure=pure)
    for key,value in expected.items():assert out[key]==value,(name,key,out)
    # Verify the coefficient-norm collar certificate algebra independently.
    for s in plan['strips']:
        q=B.polynomial(s['transformed_polynomial_terms'],(X,B.Z))
        deg=q.degree(X);f=B.polynomial(s['leading_x_coefficient_terms'],(B.Z,))
        assert S.expand(q.as_expr().coeff(X,deg)-f.as_expr())==0
        assert s['normal_radius']>S.Rational(s['leading_coefficient_root_bound'])
        assert s['collar_outer_radius']-s['normal_radius']>s['normal_step_bound']
        assert S.Rational(s['error_numerator_bound'])/s['collar_validity_start']<S.Rational(s['collar_lower_bound'])
    return dict(name=name,plan=plan,evaluated=evaluated,result=out),plan,evaluated

controls=[]
# Every monomial certificate below follows from tau_k=ceil(rk/D),not this prefix.
case,plan,evaluated=run('unary positive gap',summary((2*Y-3*X)*(2*Y-3*X-1),2),
    [[[2],'1']],[[[3],'1']],dict(selected_rational_slope='3/2',all_marker_eventual_safety=True))
controls.append(case)
case,_,_=run('lacunary formal unit relevant jet',summary((2*Y-3*X)*(2*Y-3*X-1),2),
    [[[2],'1']],[[[3],'1'],[[4],'1']],dict(selected_rational_slope='3/2',all_marker_eventual_safety=True))
controls.append(case)
case,_,_=run('unary unsafe without purity',summary((3*Y-X)*(3*Y-X-1)*(3*Y-X-2),1),
    [[[3],'1']],[[[1],'1']],dict(selected_rational_slope='1/3',all_marker_eventual_safety=False),pure=False)
controls.append(case)
case,_,_=run('inactive direction nilpotent alternative',summary(S.Integer(1),0),
    [[[1,0],'1']],[[[0,1],'1']],dict(alternative='eventually infinite first-survival orders',nilpotence_index_bound=1,all_marker_eventual_safety=True))
controls.append(case)
case,_,_=run('multiple distinct candidate slopes',summary((Y-X)*(Y-2*X),2),
    [[[1],'1']],[[[2],'1']],dict(selected_rational_slope='2',all_marker_eventual_safety=True))
controls.append(case)
case,_,_=run('unbounded nonlinear branches share the selected slope',summary((Y-2*X)*((Y-2*X)**2-X),2),
    [[[1],'1']],[[[2],'1']],dict(selected_rational_slope='2',all_marker_eventual_safety=True))
controls.append(case)
case,_,_=run('bounded nonlinear branch beside exact line',summary((Y-2*X)*(X*(Y-2*X-1)-1),2),
    [[[1],'1']],[[[2],'1']],dict(selected_rational_slope='2',all_marker_eventual_safety=True))
controls.append(case)

# Genuine all-order DVG seed: its retained identities prove purity and tau=k+1.
a,b,c,d,t=S.symbols('a b c d t')
phase=S.Poly((t+c)*(a*d+b*t),a,b,c,d,t)
symbol=[[[1,0,0,1,1],'1'],[[0,1,1,0,1],'-1']]
case,_,_=run('retained critical pure unsafe seed',summary(Y-X-1,1),B.terms(phase),symbol,
    dict(selected_rational_slope='1',all_marker_eventual_safety=False,critical_bad_marker_candidates_P_powers=[1]))
controls.append(case)

# Numerical cutoff transcription,checked against exact unary derivative degrees.
cutoffs=[]
for control in controls:
    result=control['result']
    if result.get('selected_rational_slope') not in ['3/2','2']:continue
    r,D=(3,2) if result['selected_rational_slope']=='3/2' else (2,1)
    C=result['affine_cutoff']['C'];E=result['affine_cutoff']['E']
    for q in [0,1,3]:
        m=C*q+E+1
        assert r*m>D*m+q
        cutoffs.append(dict(control=control['name'],marker_degree=q,moment=m))

rejections=[]
def rejects(name,fn):
    try:fn()
    except ValueError:rejections.append(name)
    else:raise AssertionError(name)
rejects('insufficient complete jet',lambda:B.evaluate_rational(plan,[[[2],'1']],[[[3],'1']],0))
rejects('nonzero symbol constant',lambda:B.evaluate_rational(plan,[[[2],'1']],[[[0],'1']],100))
rejects('missing original-power evaluation',lambda:B.resolve(plan,[]))
altered=copy.deepcopy(plan);altered['query_operator_power']+=1
rejects('altered derived plan',lambda:B.resolve(altered,evaluated['evaluations']))
rejects('zero order polynomial',lambda:B.make_plan(summary(S.Integer(0),1)))

# This deliberately FALSE certificate passes its one finite point. The program
# MUST NOT be described as a verifier of the actual all-index module evidence.
false_plan=B.make_plan(summary(Y-X-1,2))
false_evals=B.evaluate_rational(false_plan,[[[1],'1']],[[[2],'1']],2)
conditional=B.resolve(false_plan,false_evals['evaluations'])
assert conditional['selected_rational_slope']=='1'
assert conditional['conditional_on_supplied_actual_certificates'] is True
false_contract_control=dict(
    statement='A guessed p can pass the finite query and give a wrong conditional slope. Underlying actual Bernstein certificates and their validity bounds remain REQUIRED.',
    guessed_polynomial='Y-X-1',actual_pair='P=x,T=partial_x^2',actual_slope='2',
    conditional_slope_from_FALSE_contract=conditional['selected_rational_slope'],
    falsely_supplied_contract_accepted_as_verified=False)

# A finite adaptive query transcript cannot distinguish an inactive symbol
# from an unseen high active jet. The universal argument is in the proof.
queried={(0,1),(1,0),(2,0),(4,0),(11,0)}
M=max(e[0] for e in queried)+1
assert M>=2 and (M,0) not in queried
inactive={(0,1):1};perturbed={(0,1):1,(M,0):1}
assert all(inactive.get(e,0)==perturbed.get(e,0) for e in queried)
assert B._operator({(M,0):B.Fraction(1)},{e:B.Fraction(c) for e,c in inactive.items()})=={}
assert B._operator({(M,0):B.Fraction(1)},{e:B.Fraction(c) for e,c in perturbed.items()})!={}
oracle_control=dict(queried_exponents=[list(e) for e in sorted(queried)],unseen_active_order=M,
    statement='Identical finite coefficient transcript;inactive tau1=infinity versus perturbed tau1=M. Universal impossibility is proved by adversarial query extension,not by this example.')

# Exercise all documented CLI operations on one small exact input bundle.
first=controls[0]
example=dict(first['plan']['certificate_summary'],plan=first['plan'],
    phase_terms=[[[2],'1']],symbol_jet_terms=[[[3],'1']],complete_through=4,
    evaluations=first['evaluated']['evaluations'],all_positive_pure=True,
    expected_result=first['result'])
examplefile=p/'bernstein_slope_producer_example.json'
examplefile.write_text(json.dumps(example,indent=2)+'\n')
tmp=p.parents[2]/'research/.scratch/bernstein_slope_cli_controls_20261116';tmp.mkdir(parents=True,exist_ok=True)
cli=[]
for mode in ['plan','evaluate','resolve']:
    output=tmp/(mode+'.json')
    result=subprocess.run([sys.executable,str(p/'bernstein_slope_producer.py'),mode,str(examplefile),str(output)],capture_output=True,text=True)
    assert result.returncode==0,result.stderr
    got=json.loads(output.read_text())
    want={'plan':first['plan'],'evaluate':first['evaluated'],'resolve':first['result']}[mode]
    assert got==want,mode
    cli.append(dict(operation=mode,output_sha256=sha(output)))

receipt=dict(completed=True,kind='exact finite postprocessor and original-power controls,not universal verification',
    canonical_proof_sha256=sha(proof),producer_sha256=sha(p/'bernstein_slope_producer.py'),checker_sha256=sha(Path(__file__)),
    actual_input_controls=controls,cutoff_transcription_controls=cutoffs,rejection_controls=rejections,
    invalid_contract_countercontrol=false_contract_control,
    finite_oracle_transcript_control=oracle_control,
    cli_controls=cli,example_sha256=sha(examplefile),
    universal_strip_and_selection_theorem_verified_by='written proof only',
    underlying_general_Bernstein_certificate_generator_implemented=False,
    independent_replay=False,formal_verification=False,external_review=False)
(p/'bernstein_slope_producer_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps(dict(completed=True,actual_inputs=len(controls),cutoffs=len(cutoffs),rejections=len(rejections),invalid_contract_controls=1,
    computed_query_indices=[c['plan']['query_operator_power'] for c in controls])))
