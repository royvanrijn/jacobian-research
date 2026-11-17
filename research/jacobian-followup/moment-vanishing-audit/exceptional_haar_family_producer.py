"""Full polynomial coefficient-family Haar adapter; Gaussian backend unchanged.

All coefficients are rational polynomials in c0,...,c_(r-1). Algebraic
coefficient fields can be given by parameter_ideal. If group_root_parameter
is supplied, that parameter is the SAME group-field root; its square must
be 2 (2O) or5 (2I) modulo the ORDINARY ideal. Otherwise the group field is
used solely in the exact invariant-family check.
"""
from pathlib import Path
import argparse,hashlib,json,itertools
import sympy as sp
import gaussian_family_producer as G
import exceptional_haar_pointwise_producer as H

ROOT=Path(__file__).resolve().parent
sha=lambda p:hashlib.sha256(Path(p).read_bytes()).hexdigest()

def poly_expression(terms,variables):
 return sum(c*sp.prod(x**n for x,n in zip(variables,e)) for e,c in terms.items())

def encoded_parameter(expr,cs):
 return [[list(e),str(c)] for e,c in sp.Poly(sp.expand(expr),*cs,domain=sp.QQ).terms() if c]

def canonical_ideal(ideal,cs):
 gb=sp.groebner(ideal,*cs,domain=sp.QQ)
 return [p.as_expr() for p in gb.polys]

def prepare(data):
 allowed={'group','parameter_count','phase_terms','parameter_ideal','group_root_parameter'}
 if set(data)-allowed:raise ValueError('unknown Haar family input keys')
 r=data['parameter_count']
 if type(r) is not int or r<1:raise ValueError('positive parameter_count required')
 name=data['group'];radicand,values,generators=H.group_data(name)
 _,_,raw,ideal,_=G.normalize({'dimension':4,'parameter_count':r,'phase_terms':data['phase_terms'],'parameter_ideal':data.get('parameter_ideal',[])})
 cs=sp.symbols('c:'+str(r));xs=sp.symbols('q:4');K=sp.QQ.poly_ring(*cs)
 expressions=canonical_ideal([poly_expression(f,cs) for f in ideal],cs)
 gb=sp.groebner(expressions,*cs,domain=sp.QQ)
 reduce=lambda f:gb.reduce(sp.expand(f))[1]
 root=data.get('group_root_parameter')
 if root is not None:
  if type(root) is not int or not 0<=root<r:raise ValueError('invalid group_root_parameter')
  if reduce(cs[root]**2-radicand)!=0:raise ValueError('group root relation must hold in the ORDINARY coefficient algebra')
  domain=sp.QQ;sqrt=cs[root]
 else:domain=sp.QQ.algebraic_field(sp.sqrt(radicand));sqrt=sp.sqrt(radicand)
 rawpoly={}
 for e,c in raw.items():
  spatial=e[:4];coefficient=K.from_sympy(c*sp.prod(v**n for v,n in zip(cs,e[4:])))
  rawpoly[spatial]=rawpoly.get(spatial,K.zero)+coefficient
 canonical=H.sphere_reduce(rawpoly,K)
 canonical={e:K.from_sympy(reduce(K.to_sympy(c))) for e,c in canonical.items()}
 canonical={e:c for e,c in canonical.items() if c}
 if any(sum(e)%2 for e in canonical):raise ValueError('family is not centre-invariant modulo the ORDINARY parameter ideal')
 D=max((sum(e) for e in canonical),default=0)
 r2={tuple(2*int(i==j) for i in range(4)):K.one for j in range(4)}
 lifted={}
 for e,c in canonical.items():lifted=H.add(lifted,H.mul({e:c},H.power(r2,(D-sum(e))//2,K),K),K)
 lifted={e:K.from_sympy(reduce(K.to_sympy(c))) for e,c in lifted.items()}
 lifted={e:c for e,c in lifted.items() if c}
 expr=sum((K.to_sympy(c)*sp.prod(x**n for x,n in zip(xs,e)) for e,c in lifted.items()),sp.S.Zero)
 invgb=sp.groebner(expressions,*cs,domain=domain)
 rat=lambda q:sp.Rational(q.numerator,q.denominator)
 for gen in generators:
  matrix=H.right_matrix([rat(a)+rat(b)*sqrt for a,b in gen])
  change={x:sum(c*y for c,y in zip(row,xs)) for x,row in zip(xs,matrix)}
  difference=sp.Poly(sp.expand(expr.subs(change,simultaneous=True)-expr),*xs)
  if any(invgb.reduce(sp.expand(c))[1]!=0 for c in difference.coeffs()):
   raise ValueError('family is not invariant modulo the ORDINARY coefficient ideal; reduced invariance is insufficient')
 active=sorted({i for e in lifted for i,n in enumerate(e) if n}) or [0]
 phase=[]
 for e,c in sorted(lifted.items()):
  for pe,pc in sp.Poly(K.to_sympy(c),*cs,domain=sp.QQ).terms():phase.append([list(e[i] for i in active)+list(pe),str(pc)])
 backend={'dimension':len(active),'parameter_count':r,'phase_terms':phase,
          'parameter_ideal':[encoded_parameter(f,cs) for f in expressions]}
 return {'group':name,'group_order':len(values),'group_root_parameter':root,'group_radicand':radicand,
         'parameter_count':r,'ordinary_parameter_ideal':backend['parameter_ideal'],
         'homogeneous_lift_degree':D,'active_gaussian_coordinates':active,'gaussian_family_input':backend,
         'canonical_original_phase_terms':[[list(e),str(K.to_sympy(c))] for e,c in sorted(canonical.items())],
         'ordinary_nilpotents_retained':True,'invariance_checked_mod_ordinary_ideal':True,
         'fixed_homogeneous_degree_used_on_every_fibre':True,'numerical_cutoffs_computed':False}

def same_input(a,b):
 na,ra,pa,ia,ha=G.normalize(a);nb,rb,pb,ib,hb=G.normalize(b)
 if (na,ra,pa,ha)!=(nb,rb,pb,hb):return False
 cs=sp.symbols('c:'+str(ra))
 return canonical_ideal([poly_expression(f,cs) for f in ia],cs)==canonical_ideal([poly_expression(f,cs) for f in ib],cs)

def check_trusted_backend(backend,backend_workdir,prepared):
 if backend.get('status')!='DECIDED':raise ValueError('no completed family result')
 if not same_input(prepared['gaussian_family_input'],backend['input']):raise ValueError('family certificate has different phase, parameters or ORDINARY ideal')
 if backend['core_sha256']!=sha(G.CORE):raise ValueError('Gaussian family core hash drift')
 directory=Path(backend_workdir)
 for f,h in backend['cas_files_sha256'].items():
  if sha(directory/f)!=h:raise ValueError('family CAS artifact hash mismatch')
 if (directory/'family.m2').read_text()!=G.make_script(backend['input']):raise ValueError('family CAS script differs from supplied input')
 parsed=G.parse_output((directory/'family.stdout').read_text())
 if any(backend.get(k)!=v for k,v in parsed.items()):raise ValueError('family result differs from retained CAS stdout')
 if not backend.get('ordinary_coefficient_nilpotents_retained') or backend.get('relative_Hom_base_change_assumed'):raise ValueError('wrong family source contract')
 return True

def summarize(backend):
 return {'original_field_safe_radical_ideal':backend['safe_radical_ideal'],
         'original_polynomial_Hom_normalization_ideal_lift':backend['polynomial_Hom_normalization_ideal_lift'],
         'normalization_torsion_exponent':backend['normalization_torsion_exponent'],
         'normalization_colon_stable':backend['normalization_colon_stable'],
         'coefficient_ring_phase_safe':backend['coefficient_ring_phase_safe'],
         'purity_tested':False,'actual_MZ_family_decision':None,
         'normalization_formula':'Gamma_I(A)=union_L Ann_A(I^L); safety iff 1 belongs',
         'all_original_markers_preserved':True,'ordinary_nilpotents_retained':True,
         'universal_numerical_cutoff_computed':False,
         'assurance_boundary':'Gaussian family holonomic/characteristic-cycle and commutative-algebra operations trusted as CAS; NOT independent original-operator certificates'}

def run(data,workdir,timeout=30,backend_result=None,backend_workdir=None,prepare_only=False):
 if type(timeout) is not int or timeout<0:raise ValueError('timeout must be nonnegative; zero is uncapped')
 prep=prepare(data)
 if prepare_only:return {'status':'PREPARED','input':data,'preparation':prep}
 directory=Path(backend_workdir or workdir)
 backend=G.run(prep['gaussian_family_input'],directory,timeout) if backend_result is None else backend_result
 if backend.get('status')=='UNKNOWN':return {'status':'UNKNOWN','input':data,'preparation':prep,'gaussian_family_result':backend}
 check_trusted_backend(backend,directory,prep)
 return {'status':'DECIDED','input':data,'preparation':prep,'gaussian_family_result':backend,
         'backend_workdir':str(directory.resolve()),**summarize(backend),
         'adapter_sha256':sha(__file__),'group_generator_source_sha256':sha(ROOT/'verify_exceptional_sl2.py'),'pointwise_transfer_adapter_sha256':sha(H.__file__),
         'Gaussian_family_producer_sha256':sha(G.__file__),'Gaussian_family_core_sha256':sha(G.CORE),
         'independent_replay':False,'formal_verification':False,'external_review':False}

def verify(result):
 if result.get('status')!='DECIDED':raise ValueError('only completed results carry family decisions')
 if result['preparation']!=prepare(result['input']):raise ValueError('incorrect invariant family lift')
 for key,path in [('adapter_sha256',__file__),('group_generator_source_sha256',ROOT/'verify_exceptional_sl2.py'),('pointwise_transfer_adapter_sha256',H.__file__),('Gaussian_family_producer_sha256',G.__file__),('Gaussian_family_core_sha256',G.CORE)]:
  if result[key]!=sha(path):raise ValueError('adapter/backend source hash drift')
 check_trusted_backend(result['gaussian_family_result'],result['backend_workdir'],result['preparation'])
 if any(result.get(k)!=v for k,v in summarize(result['gaussian_family_result']).items()):raise ValueError('changed original family result')
 return True

def universal_family(group,degree):
 """Surjective linear Reynolds parametrization of the ENTIRE filtered phase space.

Parameters are redundant; the final parameter is the group-field root.
This constructs inputs only, never runs a family calculation.
"""
 if type(degree) is not int or degree<0:raise ValueError('nonnegative matrix degree required')
 D=2*(degree//2);d,values,_=H.group_data(group)
 f=H.G.Field({'kind':'number_field','minimal_polynomial':'a^2-'+str(d)});K=f.K
 # Backend root index0 is -sqrt(d); express averages using root parameter +sqrt(d).
 sqrt=-f.parse('a');rational=lambda q:K.convert(sp.Rational(q.numerator,q.denominator))
 monomials=sorted(e for e in itertools.product(range(D+1),repeat=4) if sum(e)==D)
 parameter_count=len(monomials)+1;phase=[];g=sp.Symbol('g')
 for j,e in enumerate(monomials):
  average={}
  for q in values:
   matrix=H.right_matrix([rational(a)+rational(b)*sqrt for a,b in q])
   average=H.add(average,H.substitute_linear({e:K.one},matrix,K),K)
  for ex,c in average.items():
   expr=f.expression(f.encode(c/K.convert(len(values)))).subs(f.a,-g)
   for (root_power,),coefficient in sp.Poly(expr,g,domain=sp.QQ).terms():
    pe=[0]*parameter_count;pe[j]=1;pe[-1]=root_power
    if coefficient:phase.append([list(ex)+pe,str(coefficient)])
 root_relation=[[[0]*(parameter_count-1)+[2],'1'],[[0]*parameter_count,str(-d)]]
 return {'group':group,'parameter_count':parameter_count,'group_root_parameter':parameter_count-1,
         'phase_terms':phase,'parameter_ideal':[root_relation]}

def main():
 ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('input',type=Path);ap.add_argument('output',type=Path)
 ap.add_argument('--workdir',type=Path);ap.add_argument('--timeout',type=int,default=30);ap.add_argument('--prepare-only',action='store_true')
 ap.add_argument('--gaussian-family-certificate',type=Path);ap.add_argument('--backend-workdir',type=Path)
 a=ap.parse_args();backend=None if a.gaussian_family_certificate is None else json.loads(a.gaussian_family_certificate.read_text())
 result=run(json.loads(a.input.read_text()),a.workdir or a.output.with_suffix('.cas'),a.timeout,backend,a.backend_workdir,a.prepare_only)
 a.output.write_text(json.dumps(result,indent=2)+'\n')

if __name__=='__main__':main()
