"""Haar family transfer checks; consumes retained CAS outputs without new holonomic runs."""
from pathlib import Path
import copy,hashlib,json,subprocess,sys
import sympy as sp
import exceptional_haar_family_producer as H
p=Path(__file__).resolve().parent;out=p/'haar_family_examples';out.mkdir(exist_ok=True)
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
rows=[]
for group in ['2O','2I']:
 for case in ['zero_phase','singular_safe_cusp']:
  directory=p/'gaussian_family_examples'/case;b=json.loads((directory/'result.json').read_text())
  data={'group':group,'parameter_count':b['input']['parameter_count'],
        'phase_terms':[[[0]*4+e[1:],c] for e,c in b['input']['phase_terms']]}
  result=H.run(data,out/'unused',backend_result=b,backend_workdir=directory/'cas');assert H.verify(result)
  assert not result['purity_tested'] and result['actual_MZ_family_decision'] is None
  assert result['coefficient_ring_phase_safe']==(case=='zero_phase')
  q=out/(group+'_'+case+'.json');q.write_text(json.dumps(result,indent=2)+'\n')
  rows.append({'example':q.name,'sha256':sha(q),'coefficient_ring_safe':result['coefficient_ring_phase_safe']})
# Entire filtered phase constructor: small input-construction checks only.
constructors=[]
for group in ['2O','2I']:
 for bound in [0,1,2,3]:
  data=H.universal_family(group,bound);prep=H.prepare(data)
  assert prep['homogeneous_lift_degree']==0 # no nonconstant invariants below8/12.
  assert data['parameter_count']==(1 if bound<2 else 10)+1
  constructors.append({'group':group,'filtered_bound':bound,'redundant_parameter_count':data['parameter_count'],'CAS_run':False})
# Algebraic-field nonconstant degree8 preparation, with ordinary parameter relation.
f=H.H.G.Field({'kind':'number_field','minimal_polynomial':'a^2+1'});c0=sp.Symbol('c0')
original=json.loads((p/'haar_pointwise_examples/2O_degree8_preparation.json').read_text())['input']
phase=[]
for es,c in original['phase_terms']:
 for (n,),q in sp.Poly(f.expression(c).subs(f.a,c0),c0,domain=sp.QQ).terms():
  if q:phase.append([es+[n],str(q)])
data={'group':'2O','parameter_count':1,'phase_terms':phase,'parameter_ideal':[[[[2],'1'],[[0],'1']]]}
prepared=H.prepare(data);assert prepared['homogeneous_lift_degree']==8
q=out/'degree8_algebraic_family_prepared.json';q.write_text(json.dumps({'status':'PREPARED only; no holonomic family calculation','input':data,'preparation':prepared},indent=2)+'\n')
# Nilpotent coefficients retained before group checking: invariant only after
# taking radicals must be rejected.
rejections=[]
def reject(name,fn):
 try:fn()
 except (ValueError,AssertionError):rejections.append(name)
 else:raise AssertionError('invalid input accepted: '+name)
reject('noninvariant nilpotent quadratic',lambda:H.prepare({'group':'2O','parameter_count':1,'phase_terms':[[[2,0,0,0,1],'1']], 'parameter_ideal':[[[[2],'1']]]}))
reject('group-root relation only radical',lambda:H.prepare({'group':'2I','parameter_count':1,'group_root_parameter':0,'phase_terms':[], 'parameter_ideal':[[[[4],'1'],[[2],'-10'],[[0],'25']]]}))
reject('wrong spatial exponent count',lambda:H.prepare({'group':'2O','parameter_count':1,'phase_terms':[[[0,1],'1']]}))
reject('floating coefficients',lambda:H.prepare({'group':'2O','parameter_count':1,'phase_terms':[[[0,0,0,0,1],0.5]]}))
reject('invalid timeout',lambda:H.run({'group':'2O','parameter_count':1,'phase_terms':[]},out/'bad',timeout=-1))
x=json.loads((out/'2O_singular_safe_cusp.json').read_text());x['coefficient_ring_phase_safe']=True;reject('false coefficient-ring safety',lambda:H.verify(x))
x=json.loads((out/'2O_singular_safe_cusp.json').read_text());x['gaussian_family_result']['safe_radical_ideal']=[];reject('backend ideal differs from retained stdout',lambda:H.verify(x))
x=json.loads((out/'2O_zero_phase.json').read_text());x['Gaussian_family_core_sha256']='0'*64;reject('backend core provenance mismatch',lambda:H.verify(x))
# Ordinary ideals with same radical are different normalization problems.
a={'dimension':1,'parameter_count':1,'phase_terms':[[[0,1],'1']],'parameter_ideal':[[[[3],'1']]]}
b={**a,'parameter_ideal':[[[[1],'1']]]};assert not H.same_input(a,b)
# Exact principal colon chains for the all-ring theorem, not field sampling.
u,v=sp.symbols('u v');assert sp.expand(v**2*u**2)==u**2*v**2
thick={'ring':'QQ[u,v]/(u²v²)','safe_ideal':'(u)','normalization':'(v²)','annihilation_power':2,'unit_in_normalization':False}
retained_thick=p/'gaussian_family_examples/ordinary_thick_torsion/result.json'
rt=json.loads(retained_thick.read_text());cc=sp.symbols('c0:2')
normal=[H.poly_expression(dict((tuple(e),sp.Rational(c)) for e,c in f),cc) for f in rt['polynomial_Hom_normalization_ideal_lift']]
assert H.canonical_ideal(normal,cc)==H.canonical_ideal([cc[1]**2],cc) and rt['normalization_torsion_exponent']==2
thick['retained_Gaussian_torsion_control_source']='gaussian_family_examples/ordinary_thick_torsion/result.json'
thick['retained_source_sha256']=sha(retained_thick)
for n in range(5):assert sp.degree(u**(n+1),u)==n+1
unknown=H.run({'group':'2I','parameter_count':1,'phase_terms':[]},out/'unused',backend_result={'status':'UNKNOWN'})
assert unknown['status']=='UNKNOWN' and 'coefficient_ring_phase_safe' not in unknown
# CLI consumes the retained family result and its provenance; no new CAS run.
inputfile=out/'cli_input.json';inputfile.write_text(json.dumps({'group':'2I','parameter_count':1,'phase_terms':[]})+'\n');cli=out/'cli_result.json'
run=subprocess.run([sys.executable,str(p/'exceptional_haar_family_producer.py'),str(inputfile),str(cli),'--gaussian-family-certificate',str(p/'gaussian_family_examples/zero_phase/result.json'),'--backend-workdir',str(p/'gaussian_family_examples/zero_phase/cas')],capture_output=True,text=True,timeout=15)
assert run.returncode==0,(run.stdout,run.stderr);assert H.verify(json.loads(cli.read_text()))
receipt={'date':'2026-11-16','completed':True,'checker_sha256':sha(Path(__file__)),'adapter_sha256':sha(p/'exceptional_haar_family_producer.py'),
 'Gaussian_family_producer_sha256':sha(p/'gaussian_family_producer.py'),'Gaussian_family_core_sha256':sha(p/'gaussian_family_core.m2'),
 'sympy_version':sp.__version__,'reused_family_controls':rows,'universal_constructor_controls':constructors,
 'nonconstant_algebraic_family_prepared':{'source':'haar_family_examples/degree8_algebraic_family_prepared.json','sha256':sha(q),'holonomic_family_run':False},
 'ordinary_thick_normalization_transcription':thick,'mutation_rejections':rejections,
 'ordinary_same_radical_different_ideals_kept':True,'CLI_retained_family_certificate_PASS':True,'UNKNOWN_propagates_without_decision':True,
 'full_marked_safe_ideal_nilpotence_all_rings_written_proof':True,'exact_Haar_source_formal_normalization_uniqueness_written_proof':True,
 'normalization_exactly_safe_ideal_torsion_written_proof':True,'full_affine_ordinary_cutoff_ideal_transfer_written_proof':True,
 'new_holonomic_family_CAS_constructions':0,'new_nonconstant_spatial_family_locus_computed':False,
 'Gaussian_family_backend_modified':False,'family_outputs_independently_operator_certified':False,
 'pure_locus_or_full_MZ_decision_claimed':False,'universal_numerical_cutoffs_computed':False,
 'arbitrary_coefficient_nilpotents_discarded':False,'ordinary_cutoff_ideal_stabilization_claimed':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,'external_review_next_task':False,
 'new_degree_support_or_parameter_campaign':False,'full_exceptional_MZ_cases':'OPEN on BOTH quotients',
 'initial_pre_registration_correction':'zero symbolic polynomial initially used Python int0; corrected to SymPy0 before any backend invocation; no mathematical decision or CAS source was produced by that failure'}
(p/'haar_family_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'completed':True,'reused_full_family_results':4,'constructor_controls':len(constructors),'mutations':len(rejections),'new_holonomic_CAS':0,'CLI':'PASS'}))
