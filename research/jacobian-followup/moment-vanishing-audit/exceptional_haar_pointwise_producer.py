"""Certified pointwise Haar adapter; consumes the existing whole Gaussian producer.

Input coordinates are (a,b,c,d) on a^2+b^2+c^2+d^2=1, with the explicit
standard quaternion 2O/2I groups retained in verify_exceptional_sl2.py.
No degree search or Gaussian backend mutation is performed.
"""
from pathlib import Path
from fractions import Fraction
from functools import lru_cache
import argparse,hashlib,json,math
import sympy as sp
import gaussian_source_hom_producer as G
from verify_exceptional_sl2 import group as quaternion_group

ROOT=Path(__file__).resolve().parent
ZERO=(0,0,0,0)

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def add(p,q,K):
 out=dict(p)
 for e,c in q.items():out[e]=out.get(e,K.zero)+c
 return {e:c for e,c in out.items() if c}

def mul(p,q,K):
 out={}
 for e,c in p.items():
  for f,d in q.items():
   es=tuple(x+y for x,y in zip(e,f));out[es]=out.get(es,K.zero)+c*d
 return {e:c for e,c in out.items() if c}

def power(p,n,K):
 out={ZERO:K.one}
 while n:
  if n&1:out=mul(out,p,K)
  n//=2
  if n:p=mul(p,p,K)
 return out

def encode(poly,field):return [[list(e),field.encode(c)] for e,c in sorted(poly.items()) if c]

def sphere_reduce(poly,K):
 base={ZERO:K.one, (2,0,0,0):-K.one,(0,2,0,0):-K.one,(0,0,2,0):-K.one}
 out={};cache={}
 for e,c in poly.items():
  n=e[3]//2
  if n not in cache:cache[n]=power(base,n,K)
  mon=(e[0],e[1],e[2],e[3]%2)
  out=add(out,mul({mon:c},cache[n],K),K)
 return out

@lru_cache(None)
def group_data(name):
 half=Fraction(1,2);quarter=Fraction(1,4)
 if name=='2O':d=2;extra=[[(0,half),(0,half),(0,0),(0,0)]];order=48
 elif name=='2I':d=5;extra=[[(0,0),(half,0),(quarter,-quarter),(-quarter,-quarter)]];order=120
 else:raise ValueError('group must be 2O or 2I')
 values,_,_=quaternion_group(d,extra)
 if len(values)!=order:raise ValueError('wrong generated group order')
 scalar=lambda a,b=0:(Fraction(a),Fraction(b))
 generators=[(scalar(0),scalar(1),scalar(0),scalar(0)),(scalar(0),scalar(0),scalar(1),scalar(0)),(scalar(half),)*4]
 generators += [tuple(scalar(*c) for c in q) for q in extra]
 return d,tuple(sorted(values)),tuple(generators)

def right_matrix(q):
 a,b,c,d=q
 return [[a,-b,-c,-d],[b,a,d,-c],[c,-d,a,b],[d,c,-b,a]]

def substitute_linear(poly,matrix,K):
 linear=[]
 for row in matrix:
  linear.append({tuple(int(i==j) for i in range(4)):c for j,c in enumerate(row) if c})
 cache={};out={}
 for es,c in poly.items():
  term={ZERO:c}
  for i,n in enumerate(es):
   if (i,n) not in cache:cache[i,n]=power(linear[i],n,K)
   term=mul(term,cache[i,n],K)
  out=add(out,term,K)
 return out

def prepare(data):
 if set(data)-{'group','field','phase_terms'}:raise ValueError('unknown Haar input keys')
 name=data['group'];d,values,generators=group_data(name)
 field,_,raw=G.normalize({'dimension':4,'phase_terms':data['phase_terms'],**({'field':data['field']} if 'field' in data else {})})
 K=field.K;canonical=sphere_reduce(raw,K)
 if any(sum(e)%2 for e in canonical):raise ValueError('phase is not invariant under the binary centre on the sphere quotient')
 D=max((sum(e) for e in canonical),default=0)
 r2={tuple(2*int(i==j) for i in range(4)):K.one for j in range(4)}
 lifted={}
 for e,c in canonical.items():lifted=add(lifted,mul({e:c},power(r2,(D-sum(e))//2,K),K),K)
 # Group coefficients may require an extension unrelated to the phase field.
 # Construct its exact compositum; no assumption of linear disjointness.
 L=sp.QQ.algebraic_field(*(([field.root] if field.root is not None else [])+[sp.sqrt(d)]))
 converted={e:L.from_sympy(K.to_sympy(c)) for e,c in lifted.items()}
 rational=lambda x:sp.Rational(x.numerator,x.denominator)
 for q in generators:
  qs=[L.from_sympy(rational(a)+rational(b)*sp.sqrt(d)) for a,b in q]
  if substitute_linear(converted,right_matrix(qs),L)!=converted:
   raise ValueError('phase is not invariant under the supplied exceptional group; averaging phases is not allowed')
 active=sorted({i for e in lifted for i,n in enumerate(e) if n}) or [0]
 reduced={tuple(e[i] for i in active):c for e,c in lifted.items()}
 gaussian={'dimension':len(active),'phase_terms':encode(reduced,field),'field':field.spec}
 group_rows=[[[str(a),str(b)] for a,b in q] for q in values]
 group_hash=hashlib.sha256(json.dumps(group_rows,sort_keys=True).encode()).hexdigest()
 return {'coordinate_system':'quaternion (a,b,c,d); matrix [[a+ib,c+id],[-c+id,a-ib]]',
         'sphere_relation':'a^2+b^2+c^2+d^2-1', 'group':name,'group_order':len(values),'group_quadratic_radicand':d,
         'group_quaternion_pairs':group_rows,'group_data_sha256':group_hash,
         'coefficient_field':field.spec,'canonical_original_phase_terms':encode(canonical,field),
         'homogeneous_lift_degree':D,'homogeneous_gaussian_phase_terms':encode(lifted,field),
         'active_gaussian_coordinates':active,'gaussian_input':gaussian,
         'invariance_checked_on_exact_generators':True,'quotient_representative_reduced_before_homogenization':True,
         'inactive_coordinate_all_marker_transfer':'independent Wick factorization',
         'radial_factor':'2^((D*m+e)/2)*((D*m+e)/2+1)! for a homogeneous marker of even degree e'}

def same_gaussian_input(a,b):
 fa,na,pa=G.normalize(a);fb,nb,pb=G.normalize(b)
 return fa.spec==fb.spec and na==nb and pa==pb

def summarize(prepared,backend):
 if not same_gaussian_input(prepared['gaussian_input'],backend['input']):raise ValueError('Gaussian certificate belongs to a different phase or field')
 G.verify_algebraic_parts(backend)
 f,n,_=G.normalize(backend['input']);K=f.K;source=backend['exact_unlocalized_source']
 unit=next(i for i,b in enumerate(source['amplitude_basis']) if len(b)==1 and not any(b[0][0]) and f.parse(b[0][1])==K.one)
 period=backend['actual_period_selection']['selected_original_generator_periods'][unit]
 eventual=period is not None
 unitpoly={} if period is None else {e[0]:f.parse(c) for e,c in period}
 allpositive=eventual and unitpoly=={0:K.one}
 D=prepared['homogeneous_lift_degree'];moments=None
 if eventual:
  moments=[[m,f.encode(c*K.convert(math.factorial(m))/K.convert(2**(D*m//2)*math.factorial(D*m//2+1)))] for m,c in sorted(unitpoly.items())]
 safe=backend['all_original_marker_eventual_safety'];marker=None
 if not safe:
  terms=backend['original_bad_marker_terms']
  if not terms or len(terms)!=1:raise ValueError('backend bad amplitude must be one original monomial')
  exps,c=terms[0];degree=sum(exps)
  if degree%2:raise ValueError('odd marker cannot have a nonpolynomial period for an even homogeneous phase')
  full=[0]*4
  for i,j in enumerate(prepared['active_gaussian_coordinates']):full[j]=exps[i]
  marker={'kind':'fixed regular invariant polynomial given by exact finite Reynolds orbit average',
          'formula':'Q(q)=1/|Gamma| sum_gamma R(q*gamma), restricted to r^2=1',
          'ambient_quaternion_monomial_terms':[[full,c]],'homogeneous_marker_degree':degree,
          'group':prepared['group'],'group_data_sha256':prepared['group_data_sha256'],'group_quadratic_radicand':prepared['group_quadratic_radicand'],'original_quotient_marker_filtration_degree_upper_bound':degree,
          'coefficient_field':f.spec,'group_quaternion_pairs':prepared['group_quaternion_pairs'],
          'regular_and_invariant_by_group_permutation':True,
          'infinite_original_marked_tail':'Gaussian actual nonpolynomial certificate plus invariance averaging and nonzero radial factors in EVERY order'}
 cutoff=None if not safe else {**backend['all_marker_cutoff'],'degree_convention':'original quotient matrix-entry filtration degree (equivalently quaternion degree)',
                              'transfer':'same slope/intercept; homogenize each original invariant marker at its filtration degree'}
 return {'all_original_marker_eventual_safety':safe,'eventual_purity':eventual,'all_positive_pure':allpositive,
         'actual_MZ_failure':eventual and not safe,'original_unit_moments_if_eventually_zero':moments,
         'original_fixed_bad_marker':marker,'all_original_marker_cutoff':cutoff,
         'pure_test_independent_of_zero_Hom':True,'actual_period_image_kernel_computed':False,
         'finite_bounded_prefix_used_without_indicial_certificate':False,'phase_averaged_to_force_invariance':False}

def decide(data,workdir,timeout=30,backend_result=None,prepare_only=False):
 if type(timeout) is not int or timeout<0:raise ValueError('timeout must be a nonnegative integer; zero is uncapped')
 prep=prepare(data)
 if prepare_only:return {'status':'PREPARED','input':data,'preparation':prep}
 if backend_result is None:backend=G.decide(prep['gaussian_input'],workdir,timeout)
 else:backend=backend_result
 if backend.get('status')=='UNKNOWN':return {'status':'UNKNOWN','input':data,'preparation':prep,'gaussian_result':backend,'reason':backend.get('reason','backend UNKNOWN')}
 if backend.get('status')!='DECIDED':raise ValueError('backend did not supply a decided certificate')
 result={'status':'DECIDED','input':data,'preparation':prep,'gaussian_result':backend,**summarize(prep,backend),
         'adapter_sha256':sha(__file__),'gaussian_backend_sha256':sha(G.__file__),'group_generator_source_sha256':sha(ROOT/'verify_exceptional_sl2.py'),
         'scope':'pointwise exact QQ/simple-numberfield phases in ANY degree;no whole exceptional MZ or coefficient-family decision',
         'independent_replay':False,'formal_verification':False,'external_review':False}
 return result

def verify(result):
 if result.get('status')!='DECIDED':raise ValueError('only DECIDED carries a decision certificate')
 if result['preparation']!=prepare(result['input']):raise ValueError('incorrect quotient lift or group certificate')
 if result['adapter_sha256']!=sha(__file__) or result['gaussian_backend_sha256']!=sha(G.__file__) or result['group_generator_source_sha256']!=sha(ROOT/'verify_exceptional_sl2.py'):raise ValueError('source hash drift')
 for k,v in summarize(result['preparation'],result['gaussian_result']).items():
  if result.get(k)!=v:raise ValueError('incorrect Haar decision field '+k)
 return True

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
 ap.add_argument('--workdir',type=Path);ap.add_argument('--timeout',type=int,default=30);ap.add_argument('--gaussian-certificate',type=Path);ap.add_argument('--prepare-only',action='store_true')
 a=ap.parse_args();data=json.loads(a.input.read_text());backend=None if a.gaussian_certificate is None else json.loads(a.gaussian_certificate.read_text())
 result=decide(data,a.workdir or a.output.with_suffix('.cas'),a.timeout,backend,a.prepare_only)
 a.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
