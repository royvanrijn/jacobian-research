"""Focused exact controls; the universal argument is in the companion TeX note."""
import copy
import hashlib
import importlib.util
import json
from pathlib import Path
import subprocess
import sys
import tempfile
import sympy as sp

HERE=Path(__file__).resolve().parent
SCRIPT=HERE/"gaussian_polynomial_hom_sidecheck_20261116.py"
spec=importlib.util.spec_from_file_location("hom_sidecheck", SCRIPT)
G=importlib.util.module_from_spec(spec)
spec.loader.exec_module(G)
T=sp.Symbol("t")
ONE=[[0,0,"1"]]
DER=[[0,1,"1"]]
TIME=[[1,0,"1"]]
MINUS_ONE=[[0,0,"-1"]]
ZERO=[]

def scalar(op, anchor=None):
    data=dict(generator_count=1, relations=[[op]],
              coordinate_annihilators=[dict(operator=op,left_witnesses=[ONE])])
    if anchor is not None:data['unit_anchor']=[anchor]
    return data

checks=[]
def check(label, predicate):
    assert predicate, label
    checks.append(label)

# Exact polynomial solutions of elementary differential modules.
for label,op,dim in [
    ('Gaussian zero phase',DER,1),
    ('Gaussian constant nonzero phase',[[0,1,'1'],[0,0,'-1']],0),
    ('Gaussian nonisotropic linear phase',[[0,1,'1'],[1,0,'-1']],0),
    ('Gaussian quadratic unit period',[[0,1,'1'],[1,1,'-2'],[0,0,'-1']],0),
    ('origin point module',TIME,0),
    ('nonorigin point module',[[1,0,'1'],[0,0,'-1']],0),
    ('two scalar formal solutions with spurious polynomial quotient',[[0,2,'1'],[0,1,'-1']],1),
]:
    r=G.solve(scalar(op));check(label,r['polynomial_hom_dimension']==dim)

# An unbounded all-N family: the written leading-term proof is universal.
for n in [0,1,2,5,11,23]:
    op=[[1,1,'1'],[0,0,str(-n)]]
    anchor=[[0,n,str(sp.Rational(1,sp.factorial(n)))]]
    r=G.solve(scalar(op,anchor));check(f'Euler degree {n}',r['polynomial_hom_dimension']==1 and r['coordinate_degree_certificates'][0]['degree_bound']==n)
    f=sum(sp.Rational(c)*T**i for i,c in enumerate(r['polynomial_hom_basis'][0][0]['coefficients']))
    check(f'Euler independent derivative check {n}',sp.expand(T*sp.diff(f,T)-n*f)==0 and sp.degree(f,T)==n)
    check(f'Euler anchor {n}',r['unit_values_at_zero']==['1'])

# Two generators and a cross relation: (f0',f1')=(f1,0).
rows=[[DER,MINUS_ONE],[ZERO,DER]]
data=dict(generator_count=2,relations=rows,
 coordinate_annihilators=[
   dict(operator=[[0,2,'1']],left_witnesses=[DER,ONE]),
   dict(operator=DER,left_witnesses=[ZERO,ONE])])
r=G.solve(data);check('cross-relation Hom dimension',r['polynomial_hom_dimension']==2)
for i,basis in enumerate(r['polynomial_hom_basis']):
    f=[sum(sp.Rational(c)*T**j for j,c in enumerate(p['coefficients'])) for p in basis]
    check(f'cross-relation direct replay {i}',sp.expand(sp.diff(f[0],T)-f[1])==0 and sp.diff(f[1],T)==0)

# Polynomial module plus time-surjective point torsion retains just one map.
mixed=dict(generator_count=2, relations=[[DER,ZERO],[ZERO,TIME]],
 coordinate_annihilators=[dict(operator=DER,left_witnesses=[ONE,ZERO]),dict(operator=TIME,left_witnesses=[ZERO,ONE])])
r=G.solve(mixed);check('polynomial plus punctual module',r['polynomial_hom_dimension']==1 and r['unit_values_at_zero']==['1'])

# Left witness orientation matters: D*t - t*D = 1 presents the zero module.
zero_module=dict(generator_count=1,relations=[[TIME],[DER]],coordinate_annihilators=[dict(operator=ONE,left_witnesses=[DER,[[1,0,'-1']]])])
r=G.solve(zero_module);check('left noncommutative witness',r['polynomial_hom_dimension']==0)

# Pairwise disjoint annihilator degree candidates force zero only after coupling.
no_map=dict(generator_count=2,relations=[[[[1,1,'1']],ZERO],[ZERO,[[1,1,'1'],[0,0,'-1']]],[ONE,MINUS_ONE]],
 coordinate_annihilators=[dict(operator=[[1,1,'1']],left_witnesses=[ONE,ZERO,ZERO]),dict(operator=[[1,1,'1'],[0,0,'-1']],left_witnesses=[ZERO,ONE,ZERO])])
r=G.solve(no_map);check('cross relation kills separate candidates',r['polynomial_hom_dimension']==0)

# Automatic input-derived relation witnesses are produced, then replayed.
for label,data0,dim in [
 ('automatic zero phase',scalar(DER),1),
 ('automatic linear Gaussian',scalar([[0,1,'1'],[1,0,'-1']]),0),
 ('automatic Euler',scalar([[1,1,'1'],[0,0,'-11']]),1),
 ('automatic mixed punctual module',mixed,1),
 ('automatic coupled equations',data,2),
 ('automatic noncommutative zero module',zero_module,0),
 ('automatic cross-relation zero',no_map,0)]:
    d=copy.deepcopy(data0);d.pop('coordinate_annihilators')
    r=G.solve_presentation(d)
    check(label,r['outcome']=='COMPLETE' and r['polynomial_hom_dimension']==dim)
    enriched=dict(d,coordinate_annihilators=r['generated_coordinate_annihilators'])
    replay=G.solve(enriched)
    check(label+' witness-only replay',replay['polynomial_hom_dimension']==dim)
r=G.solve_presentation(dict(generator_count=1,relations=[]),max_weyl_degree=2)
check('nonholonomic free module is UNKNOWN',r['outcome']=='UNKNOWN' and not r['gaussian_safety_claimed'])
r=G.solve_presentation(dict(generator_count=1,relations=[[[[0,1,'1']]]]),max_witness_columns=1)
check('witness column cap is UNKNOWN',r['outcome']=='UNKNOWN')

# Exact rejection and capped computation are not a safety decision.
for label,alter in [
 ('false annihilator witness',lambda d:d['coordinate_annihilators'][0].update(left_witnesses=[ZERO])),
 ('false annihilator operator',lambda d:d['coordinate_annihilators'][0].update(operator=ONE)),
 ('zero annihilator',lambda d:d['coordinate_annihilators'][0].update(operator=ZERO))]:
    d=scalar(DER);alter(d)
    try:G.solve(d)
    except ValueError:check(label,True)
    else:raise AssertionError(label)
r=G.solve(scalar([[1,1,'1'],[0,0,'-23']]),max_unknowns=3)
check('cap returns UNKNOWN without safety',r['outcome']=='UNKNOWN' and not r['gaussian_safety_claimed'])

# The general Taylor/coinduction formula on O and a rank-one exponential.
# These finite checks exercise the commutator transcription, not the universal lemma.
for n in [0,1,3,7]:
    f=T**n+T+1
    l=lambda h:sp.expand(h).subs(T,0)
    image=sum(l(sp.diff(f,T,m))*T**m/sp.factorial(m) for m in range(n+3))
    check(f'Taylor reconstruction O degree {n}',sp.expand(image-f)==0)
    for m in [1,2,5]:
        check(f'Taylor t commutator {n},{m}',sp.expand(sp.diff(T*f,T,m)-T*sp.diff(f,T,m)-m*sp.diff(f,T,m-1))==0)

# Exercise JSON CLI with real newlines and round-trip the evidence.
with tempfile.TemporaryDirectory(prefix='gaussian-hom-sidecheck-') as temp:
    source=Path(temp)/'input.json';dest=Path(temp)/'output.json'
    cli_data=copy.deepcopy(mixed);cli_data.pop("coordinate_annihilators")
    source.write_text(json.dumps(cli_data))
    cp=subprocess.run([sys.executable,str(SCRIPT),str(source),str(dest)],capture_output=True,text=True)
    check('CLI exit',cp.returncode==0)
    check('CLI semantic result',json.loads(dest.read_text())['polynomial_hom_dimension']==1)

# Actual-map jet selection: stronger than a one-dimensional formal space.
SELECTED=HERE/'gaussian_selected_jets_sidecheck_20261116.py'
spec=importlib.util.spec_from_file_location('selected_side',SELECTED)
S=importlib.util.module_from_spec(spec);spec.loader.exec_module(S)
spurious=scalar([[0,2,'1'],[0,1,'-1']])
actual=S.gaussian_amplitude_jets([[[0],'1']],[[[[0,0],'1']]],[1],[[1]])
check('original exp(t) jets',actual==[['1','1']])
r=S.match_actual_jets(spurious,actual)
check('actual exp(t) rejects spurious polynomial quotient',r['outcome']=='NO_MATCH' and sp.Rational(r['nonzero_original_jet_value'])!=0)
check('finite original-jet obstruction is explicit',len(r['polynomial_map_jet_obstruction'])==2)

# Split Gaussian P=Z has exact source O, with W=t in that quotient.
# The supplied two-generator presentation explicitly retains W as well as1.
linear_cover=dict(generator_count=2,relations=[[DER,ZERO],[[[1,0,'-1']],ONE]])
hr=G.solve_presentation(linear_cover)
check('split Gaussian all-amplitude cover witnesses',hr['outcome']=='COMPLETE' and hr['polynomial_hom_dimension']==1)
certs=hr['generated_coordinate_annihilators']
degrees=[S.origin_certificate(G.parse_op(c['operator']))['jet_degree'] for c in certs]
phase=[[[1,0],'1']]
amps=[[[[0,0,0],'1']],[[[0,0,1],'1']]]
jets=S.gaussian_amplitude_jets(phase,amps,degrees,[[0,1],[1,0]])
r=S.match_actual_jets(linear_cover,jets,hr)
check('nonzero pure safe split Gaussian selected map',r['outcome']=='MATCH')
check('split W original marker has image t',r['selected_polynomial_images'][1]['expression']=='t')

# A genuine overpresentation with two formal directions may still work.
overcover=dict(generator_count=2,relations=[[DER,ZERO],[ZERO,[[0,1,'1'],[0,0,'-1']]]])
hr=G.solve_presentation(overcover)
cs=hr['generated_coordinate_annihilators'];degrees=[S.origin_certificate(G.parse_op(c['operator']))['jet_degree'] for c in cs]
jets=S.gaussian_amplitude_jets([],[[[[0,0],'1']],[]],degrees,[[1]])
r=S.match_actual_jets(overcover,jets,hr)
check('larger formal solution space with correct actual map',r['outcome']=='MATCH' and r['selected_polynomial_images'][0]['expression']=='1' and r['selected_polynomial_images'][1]['expression']=='0')

# A coordinate's resonance may require an arbitrarily late original jet.
for n in [0,3,12]:
    d=scalar([[1,1,'1'],[0,0,str(-n)]])
    jets=S.gaussian_amplitude_jets([],[[[[n,0],'1']]],[n],[[1]])
    check(f'Wick time-amplitude jet {n}',jets[0]==['0']*n+['1'])
    r=S.match_actual_jets(d,jets)
    check(f'origin resonance selected polynomial {n}',r['outcome']=='MATCH' and r['origin_certificates'][0]['jet_degree']==n)
    if n:
        try:S.match_actual_jets(d,[jets[0][:-1]])
        except ValueError:check(f'origin bound insufficient jet rejected {n}',True)
        else:raise AssertionError('Insufficient jets were accepted')

# Gaussian Wick engine checked against independent closed formulas.
standard=S.gaussian_amplitude_jets([[[1],'1']],[[[[0,0],'1']]],[6],[[1]])[0]
check('standard Gaussian linear exponential jets',standard==['1','0','1/2','0','1/8','0','1/48'])
split=S.gaussian_amplitude_jets([[[1,1],'1']],[[[[0,0,0],'1']]],[5],[[0,1],[1,0]])[0]
check('split ZW has geometric series coefficients',split==['1']*6)
geom=scalar([[0,1,'1'],[1,1,'-1'],[0,0,'-1']])
r=S.match_actual_jets(geom,[split])
check('nonpolynomial split moment map rejected',r['outcome']=='NO_MATCH')
try:S.gaussian_amplitude_jets([],[],[],[[1,2],[0,1]])
except ValueError:check('invalid covariance rejected',True)
else:raise AssertionError('Nonsymmetric covariance accepted')

# End-to-end ORIGINAL polynomial IBP, not a trusted integration presentation.
IBP_SCRIPT=HERE/'gaussian_ibp_annihilator_sidecheck_20261116.py'
spec=importlib.util.spec_from_file_location('original_ibp',IBP_SCRIPT)
IBP=importlib.util.module_from_spec(spec);spec.loader.exec_module(IBP)
ibp_examples=[]
for label,phase,alpha,C,expected in [
 ('zero Gaussian phase',[],[0],[[1]],'MATCH'),
 ('standard linear phase',[[[1],'1']],[0],[[1]],'NO_MATCH'),
 ('standard quadratic phase',[[[2],'1']],[0],[[1]],'NO_MATCH'),
 ('split safe nonzero marker',[[[1,0],'1']],[0,1],[[0,1],[1,0]],'MATCH'),
 ('split safe zero marker',[[[1,0],'1']],[1,0],[[0,1],[1,0]],'MATCH'),
 ('split ZW nonpolynomial unit',[[[1,1],'1']],[0,0],[[0,1],[1,0]],'NO_MATCH'),
]:
    result=IBP.find(phase,alpha,C,max_degree=3)
    check(label+' original-IBP result',result['outcome']=='COMPLETE' and result['actual_period_test']['outcome']==expected)
    # Reconstruct identity from emitted JSON with an independent expression path.
    xs=sp.symbols('u0:'+str(len(alpha)));t=sp.Symbol('v');variables=(t,*xs)
    P=sum(sp.Rational(c)*sp.prod(x**e for x,e in zip(xs,ex)) for ex,c in phase)
    Q=sp.prod(x**e for x,e in zip(xs,alpha));A=result['scalar_annihilator']
    lhs=sum(sp.Rational(c)*t**a*P**b*Q for a,b,c in A)
    inv=sp.Matrix(C).inv();rhs=0
    for i,terms in enumerate(result['polynomial_ibp_witnesses']):
        B=sum(sp.Rational(c)*sp.prod(x**e for x,e in zip(variables,ex)) for ex,c in terms)
        rhs+=sp.diff(B,xs[i])+(-sum(inv[i,j]*xs[j] for j in range(len(xs)))+t*sp.diff(P,xs[i]))*B
    check(label+' independent original identity replay',sp.expand(lhs-rhs)==0)
    ibp_examples.append(dict(name=label,phase_terms=phase,marker_exponent=alpha,covariance=C,certificate=result))
result=IBP.find([[[3],'1']],[0],[[1]],max_degree=0)
check('IBP exhaustion returns UNKNOWN',result['outcome']=='UNKNOWN')
(HERE/'gaussian-ibp-sidecheck-examples-2026-11-16.json').write_text(json.dumps(ibp_examples,indent=2)+'\n')

# Complete recurrence/low-marker pipeline, with independent exponential replay.
COVER_SCRIPT=HERE/'gaussian_euler_cover_sidecheck_20261116.py'
spec=importlib.util.spec_from_file_location('euler_cover',COVER_SCRIPT)
COVER=importlib.util.module_from_spec(spec);spec.loader.exec_module(COVER)
cover_examples=[]
for label,phase,C,expected in [
 ('zero',[],[[1]],'SAFE'),
 ('linear',[[[1],'1']],[[1]],'UNSAFE'),
 ('quadratic',[[[2],'1']],[[1]],'UNSAFE'),
 ('nonzero pure safe split phase',[[[1,0],'1']],[[0,1],[1,0]],'SAFE')]:
    result=COVER.decide(phase,C,cover_degree_cap=3,ibp_degree_cap=3)
    check(label+' full certificate decision',result['outcome']==expected)
    cert=result['cover'];n=len(C);xs=sp.symbols('q0:'+str(n));t=sp.Symbol('s')
    P=sum(sp.Rational(c)*sp.prod(x**e for x,e in zip(xs,ex)) for ex,c in phase)
    vector=sp.Matrix(xs);F=t*P-(vector.T*sp.Matrix(C).inv()*vector)[0]/2
    exponential=sp.exp(F)
    total=0;power=exponential
    for coefficient in cert['b_coefficients']:
        total+=sp.Rational(coefficient)*power
        power=-n*power-sum(x*sp.diff(power,x) for x in xs)
    for term in cert['remainder_terms']:
        ex=term['powers'];alpha,a,beta,b=ex[:n],ex[n],ex[n+1:2*n+1],ex[-1]
        value=sp.diff(exponential,t,b)
        for x,count in zip(xs,beta):value=sp.diff(value,x,count)
        total+=sp.Rational(term['coefficient'])*t**a*sp.prod(x**e for x,e in zip(xs,alpha))*value
        check(label+' positive spatial drop '+str(ex),sum(beta)>sum(alpha))
    check(label+' direct exponential recurrence replay',sp.simplify(total/exponential)==0)
    if expected=='UNSAFE':
        check(label+' actual bad marker retained','fixed_bad_marker_exponent' in result)
    cover_examples.append(dict(name=label,certificate=result))
check('split phase numerical all-marker cutoff',cover_examples[-1]['certificate']['affine_cutoff']['slope']=='1' and cover_examples[-1]['certificate']['affine_cutoff']['intercept']=='0')
result=COVER.find_cover([],[[1]],max_degree=0)
check('cover exhaustion UNKNOWN',result['outcome']=='UNKNOWN')
(HERE/'gaussian-euler-cover-examples-2026-11-16.json').write_text(json.dumps(cover_examples,indent=2)+'\n')

receipt=dict(outcome='PASS',check_count=len(checks),checks=checks,
 code_sha256=hashlib.sha256(SCRIPT.read_bytes()).hexdigest(),
 selected_jets_sha256=hashlib.sha256(SELECTED.read_bytes()).hexdigest(),
 original_ibp_sha256=hashlib.sha256(IBP_SCRIPT.read_bytes()).hexdigest(),
 euler_cover_sha256=hashlib.sha256(COVER_SCRIPT.read_bytes()).hexdigest(),
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 mathematical_scope='Exact polynomial Hom with generated left witnesses; selected-map original Gaussian jets and rational symmetric covariance. Both infinity bounds and origin injectivity are written arguments. Actual amplitude-map validity and full coverage remain external obligations',
 full_Gaussian_closure=False,certificate_only_pointwise_decision_controls=True,independent_full_proof_review=False)
output=HERE/'gaussian-polynomial-hom-sidecheck-2026-11-16.json'
output.write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','full_Gaussian_closure']}))
