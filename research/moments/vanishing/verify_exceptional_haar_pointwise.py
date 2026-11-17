"""Finite certificate/transfer controls. Reuses retained backend certificates; no CAS campaign."""
from pathlib import Path
import copy,hashlib,json,subprocess,sys
import sympy as sp
import exceptional_haar_pointwise_producer as H
p=Path(__file__).resolve().parent;out=p/'haar_pointwise_examples';out.mkdir(exist_ok=True)
sha=lambda q:hashlib.sha256(q.read_bytes()).hexdigest()
backend={k:json.loads((p/'gaussian_source_hom_examples'/k/'result.json').read_text()) for k in ['zero','constant']}
rows=[]
for name in ['2O','2I']:
 for tag,terms in [('zero',[]),('constant',[[[0,0,0,0],'1']]),('sphere_zero',[[[2,0,0,0],'1'],[[0,2,0,0],'1'],[[0,0,2,0],'1'],[[0,0,0,2],'1'],[[0,0,0,0],'-1']])]:
  data={'group':name,'phase_terms':terms};b=backend['constant' if tag=='constant' else 'zero']
  result=H.decide(data,out/(name+'_'+tag),backend_result=b);assert H.verify(result)
  assert result['all_original_marker_eventual_safety']==(tag!='constant')
  assert result['eventual_purity']==(tag!='constant') and not result['actual_MZ_failure']
  if tag=='constant':
   assert result['original_fixed_bad_marker']['ambient_quaternion_monomial_terms']==[[[0,0,0,0],'1']]
  q=out/(name+'_'+tag+'.json');q.write_text(json.dumps(result,indent=2)+'\n');rows.append({'example':q.name,'sha256':sha(q),'status':'DECIDED','safe':result['all_original_marker_eventual_safety']})
# Actual nonconstant invariant preparation; no backend CAS invocation.
f=H.G.Field({'kind':'number_field','minimal_polynomial':'a^2+1'});a=f.a;v=sp.symbols('q0:4');X=v[0]-a*v[1];Y=v[2]-a*v[3]
P=sp.Poly(sp.expand(X**8+14*X**4*Y**4+Y**8),*v)
data={'group':'2O','field':f.spec,'phase_terms':[[list(e),f.encode(f.parse(str(c)))] for e,c in P.terms()]}
prepared=H.prepare(data);assert prepared['homogeneous_lift_degree']==8
_,_,raw=H.G.normalize({'dimension':4,'field':f.spec,'phase_terms':data['phase_terms']})
assert prepared['homogeneous_gaussian_phase_terms']==H.encode(raw,f)
q=out/'2O_degree8_preparation.json';q.write_text(json.dumps({'input':data,'preparation':prepared,'status':'PREPARED; Gaussian CAS NOT invoked'},indent=2)+'\n');rows.append({'example':q.name,'sha256':sha(q),'status':'PREPARED only'})
# Homogeneous marker convention: Reynolds(X0^2)=r^2/4 on BOTH groups.
marker_controls=[]
for name in ['2O','2I']:
 d,values,_=H.group_data(name);K=sp.QQ.algebraic_field(sp.sqrt(d));avg={}
 rat=lambda c:sp.Rational(c.numerator,c.denominator)
 for q0 in values:
  qs=[K.from_sympy(rat(x)+rat(y)*sp.sqrt(d)) for x,y in q0]
  avg=H.add(avg,H.substitute_linear({(2,0,0,0):K.one},H.right_matrix(qs),K),K)
 avg={e:c/K.convert(len(values)) for e,c in avg.items()}
 expected={tuple(2*int(i==j) for i in range(4)):K.convert(sp.Rational(1,4)) for j in range(4)}
 assert avg==expected
 assert 2**1*sp.factorial(2)*sp.Rational(1,4)==1 # Gaussian E[X0²]=1.
 marker_controls.append({'group':name,'ambient_marker_degree':2,'Reynolds_marker':'r²/4','original_sphere_marker':'1/4','radial_factor':4,'Gaussian_moment':1})
rejections=[]
def rejected(label,fn):
 try:fn()
 except (ValueError,AssertionError,StopIteration):rejections.append(label)
 else:raise AssertionError('mutation accepted: '+label)
rejected('odd noninvariant original phase',lambda:H.prepare({'group':'2O','phase_terms':[[[1,0,0,0],'1']]}))
rejected('even but noninvariant original phase',lambda:H.prepare({'group':'2I','phase_terms':[[[2,0,0,0],'1']]}))
rejected('wrong phase Gaussian certificate',lambda:H.decide({'group':'2O','phase_terms':[]},out/'bad',backend_result=backend['constant']))
rejected('wrong coefficient field Gaussian certificate',lambda:H.decide({'group':'2I','field':{'kind':'number_field','minimal_polynomial':'a^2-5'},'phase_terms':[]},out/'bad',backend_result=backend['zero']))
rejected('invalid timeout',lambda:H.decide({'group':'2O','phase_terms':[]},out/'bad',timeout=-1,backend_result=backend['zero']))
base=json.loads((out/'2O_constant.json').read_text())
x=copy.deepcopy(base);x['actual_MZ_failure']=True;rejected('nonpure unsafe phase mislabeled MZ failure',lambda:H.verify(x))
x=copy.deepcopy(base);x['original_fixed_bad_marker']['group_quaternion_pairs'].pop();rejected('incomplete fixed Reynolds marker',lambda:H.verify(x))
x=copy.deepcopy(base);x['preparation']['homogeneous_lift_degree']=2;rejected('wrong radial homogenization degree',lambda:H.verify(x))
x=copy.deepcopy(base);x['gaussian_result']['all_original_marker_eventual_safety']=True;rejected('mutated backend safety certificate',lambda:H.verify(x))
x=json.loads((out/'2O_zero.json').read_text());x['all_original_marker_cutoff']['slope']='1';rejected('wrong original cutoff',lambda:H.verify(x))
unknown=H.decide({'group':'2I','phase_terms':[]},out/'unused',backend_result={'status':'UNKNOWN','reason':'retained backend cap'})
assert unknown['status']=='UNKNOWN' and 'actual_MZ_failure' not in unknown
# CLI certificate consumption, no CAS invocation.
inputfile=out/'cli_input.json';inputfile.write_text(json.dumps({'group':'2I','phase_terms':[]})+'\n')
cli=out/'cli_result.json';r=subprocess.run([sys.executable,str(p/'exceptional_haar_pointwise_producer.py'),str(inputfile),str(cli),'--gaussian-certificate',str(p/'gaussian_source_hom_examples/zero/result.json')],capture_output=True,text=True,timeout=15)
assert r.returncode==0,(r.stdout,r.stderr);assert H.verify(json.loads(cli.read_text()))
receipt={'date':'2026-11-16','completed':True,'checker_sha256':sha(Path(__file__)),'adapter_sha256':sha(p/'exceptional_haar_pointwise_producer.py'),'gaussian_backend_sha256':sha(p/'gaussian_source_hom_producer.py'),'sympy_version':sp.__version__,
 'decided_reused_certificate_controls':rows[:6],'nonconstant_original_invariant_preparation':rows[6],
 'homogeneous_Reynolds_marker_radial_factor_controls':marker_controls,'mutation_rejections':rejections,
 'CLI_retained_certificate_consumption_pass':True,'UNKNOWN_propagates_without_decision':True,
 'new_CAS_constructions':0,'new_Gaussian_source_algorithm_implemented':False,
 'pointwise_ALL_degree_input_and_ALL_marker_certificate_transfer_written_proof':True,
 'fixed_original_marker_kept_homogeneous':True,'original_filtration_cutoff_transferred':True,
 'unit_period_selected_independently_of_zero_Hom':True,
 'new_nonconstant_phase_safety_decisions_run':False,'full_phase_family_decision_implemented':False,
 'actual_period_kernel_computed':False,'arbitrary_transcendental_coefficients_implemented':False,
 'new_pure_exceptional_failure_produced':False,'full_exceptional_purity_implies_extension_zero_proved':False,
 'independent_replay':False,'formal_verification':False,'external_review':False,'external_review_next_task':False,
 'new_degree_support_or_parameter_campaign':False,'full_exceptional_MZ_cases':'OPEN on BOTH quotients'}
(p/'haar_pointwise_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'completed':True,'decided_certificate_controls':6,'nonconstant_preparation':1,'mutation_rejections':len(rejections),'new_CAS_constructions':0,'CLI':'PASS'}))
