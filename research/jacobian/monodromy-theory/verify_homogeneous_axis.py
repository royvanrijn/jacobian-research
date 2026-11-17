"""Exact homogeneous-axis and exceptional11 decomposition controls KMFI118--120."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
import sympy as S
import json,hashlib,time,platform
from collections import Counter
from homogeneous_axis import norm_control,quotient_rows
from sweep_differential import gradient,coframe
from psl11_residues import exceptional_actions,inverse
from psl11_curve_spectrum import integer_group
from local_image import orbit_lengths
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,truth):
    if not truth:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
normmodels=[]
for r in [2,3,4]:
    d=norm_control(r);ps=d['parameters'];M=d['M'];F=d['F'];f=d['f'];g=d['g'];X=d['X'];n=r
    T=d['T'];A,Z,G=S.symbols('A Z gamma')
    check('rank'+str(r)+' separable rational companion polynomial',
       S.degree(d['d'],T)==r and S.gcd(S.Poly(d['d'],T),S.Poly(S.diff(d['d'],T),T)).degree()==0)
    check('rank'+str(r)+' roots avoid0 and1',d['d'].subs(T,0)!=0 and d['d'].subs(T,1)!=0)
    check('rank'+str(r)+' constant rank-one nilpotency identity',
       S.expand((r+1)*d['d']-T*S.diff(d['d'],T)-(T-1)**r)==0)
    check('rank'+str(r)+' norm form homogeneous',all(sum(mon)==r for mon,co in S.Poly(F,*ps).terms()))
    # These whole polynomial moment identities verify the diagonal norm reduction
    # used to turn the constant characteristic identity into global nilpotency.
    for j in range(r+1):
        check('rank'+str(r)+' entire norm gradient moment'+str(j),
          zero((gradient(F,ps)*(M**j)*S.Matrix(ps))[0]-S.trace(M**j)*F))
    P=(g.jacobian(ps)-S.Matrix(ps)*gradient(f,ps)).applyfunc(S.expand)
    H=g-f*S.Matrix(ps)
    check('rank'+str(r)+' whole homogeneous Euler identity',
          mz(P*S.Matrix(ps)-(n+1)*H-f*S.Matrix(ps)))
    check('rank'+str(r)+' whole matrix homogeneous',all(sum(mon)==n for v in P for mon,co in S.Poly(v,*ps).terms() if co!=0))
    check('rank'+str(r)+' H exact norm multiplication',mz(H-F*M*S.Matrix(ps)))
    pt={p:(1 if i==0 else 0) for i,p in enumerate(ps)};P0=P.subs(pt);df0=gradient(f,ps).subs(pt)
    C0=S.Matrix.vstack(*[df0*P0**j for j in range(r)])
    check('rank'+str(r)+' actual full-rank coframe atunit',C0.det()!=0)
    check('rank'+str(r)+' actual nilpotent regular Jordan type atunit',
        P0**r==S.zeros(r) and P0.rank()==r-1)
    source=X+G*S.Matrix([1,*ps])
    check('rank'+str(r)+' whole gamma-polynomial Jacobian atunit',
       zero(source.jacobian([G,*ps]).subs(pt).det()-G**r))
    if r<=3:
        check('rank'+str(r)+' entire polynomial nilpotency',mz((P**r).applyfunc(S.expand)))
    V=S.Matrix(ps)+F*M*S.Matrix(ps)
    cert=(d['Q'].adjugate()*V).applyfunc(S.expand)
    want=S.zeros(r,1);want[0]=F;want[1]=F**2
    check('rank'+str(r)+' entire axis ideal adjugate certificate',mz(cert-want))
    check('rank'+str(r)+' exact individual axis ideal reductions',
       all(zero(V[i]-F*(M*S.Matrix(ps))[i]-ps[i]) for i in range(r)))
    det=S.expand((A*S.eye(r)+Z*M).det())
    check('rank'+str(r)+' exact binomial norm eliminant',
       zero(Z*det-((A+Z)**(r+1)-A**(r+1))/S.Integer(r+1)))
    YY=(A*S.eye(r)+Z*M)*S.Matrix(ps)
    QY=sum((YY[j]*M**j for j in range(r)),S.zeros(r))
    check('rank'+str(r)+' whole multiplication-algebra target identity',
       mz(QY-(A*S.eye(r)+Z*M)*d['Q']))
    normmodels.append({'r':r,'matrix_degree':r,'norm_terms':len(S.Poly(F,*ps).terms()),
      'companion_polynomial':str(d['d']),'M':[[str(e) for e in row] for row in M.tolist()],
      'F':str(F),'X':[str(e) for e in X],'coframe_atunit':str(C0.det()),
      'actual_axis_ideal':'(p1,...,pr)','actual_axis_points':1,'generic_degree':r+1,
      'geometric_group':'regular C'+str(r+1),'arithmetic_constants':'Q(zeta_'+str(r+1)+')',
      'Keller_completion_allowed':False})
# Positive axis orbits occur at rank-drop points. Exact triangular controls.
positive=[]
for r,n in [(2,2),(2,3),(4,2),(4,5)]:
    ps=list(S.symbols('x1:'+str(r+1)));z=ps[-1];f=z**n
    g=S.zeros(r,1);g[-1]=S.Rational(n,n+1)*z**(n+1)
    P=(g.jacobian(ps)-S.Matrix(ps)*gradient(f,ps)).applyfunc(S.expand)
    H=g-f*S.Matrix(ps);V=H+S.Matrix(ps)
    basis=S.groebner(list(V),*ps,order='lex',domain=S.QQ)
    expected=[*ps[:-1],z**(n+1)-(n+1)*z]
    check('positive rank'+str(r)+' degree'+str(n)+' entire axis ideal',
      all(basis.reduce(v)[1]==0 for v in expected) and
      all(S.groebner(expected,*ps,order='lex',domain=S.QQ).reduce(v)[1]==0 for v in V))
    check('positive rank'+str(r)+' degree'+str(n)+' finite reduced axis algebra',
      S.gcd(S.Poly(expected[-1],z),S.Poly(S.diff(expected[-1],z),z)).degree()==0)
    nonzero={p:0 for p in ps[:-1]}
    tail=P.subs(nonzero)
    check('positive rank'+str(r)+' degree'+str(n)+' nonzero axis rank collapse',
      tail==S.zeros(r) and S.gcd(S.Poly(z**n-(n+1),z),S.Poly(S.diff(f,z),z)).degree()==0)
    check('positive rank'+str(r)+' degree'+str(n)+' normalized matrix nilpotency',mz(P**r))
    positive.append({'rank':r,'matrix_degree':n,'k':n+1,'generic_rank_DX':1,
                     'nonzero_axis_rank_P':0,'axis_ideal':[str(v) for v in expected]})
# Check the residue contradiction and its essential rank-one exception.
residues=[]
for r in range(1,7):
 for n in [1,2,5]:
    residual=1+S.Rational((n+1)**r-1,(-n)**r)
    check('Euler Jacobi hypothetical finite rank'+str(r)+' n'+str(n),
          residual==0 if r==1 else residual!=0)
    residues.append({'r':r,'n':n,'hypothetical_finite_residue_sum':str(residual),
                     'contradiction':r>=2})

data,orders,triangles,a5classes,actions=exceptional_actions()
elements,one,table=integer_group(data);idx={g:i for i,g in enumerate(elements)}
inverses=[idx[data['canon'](inverse(g))] for g in elements]
groupmodels=[];expected_counts={(2,1,2,11):15,(4,2,2,3):15,(3,1,3,11):10,(5,1,5,11):6}
expected_k={(2,1):[3,5,7,9],(4,2):[3],(3,1):[4,7,10],(5,1):[6,11]}
for j,action in enumerate(actions):
    H=frozenset(idx[g] for g in action['H']);perms=[action['permutations'][g] for g in elements]
    groups,pairs,accepted=quotient_rows(H,one,table,inverses,perms)
    label='exceptional_action'+str(j)
    check(label+' A5 order60 fixes distinguished sheet',len(H)==60 and
       sum(all(perms[g][i]==i for g in H) for i in range(11))==1)
    check(label+' all59 A5 subgroups',len(groups)==59 and
       Counter(map(len,groups))=={1:1,2:15,3:10,4:5,5:6,6:10,10:6,12:5,60:1})
    check(label+' all67 nontrivial cyclic quotient pairs',len(pairs)==67)
    rows=Counter((a['D_order'],a['I_order'],a['n'],len(a['I_fixed_sheets'])) for a in accepted)
    check(label+' complete46 accepted pairs and four rows',rows==expected_counts and len(accepted)==46)
    reject=[a for a in pairs if a not in accepted and len(a['I_fixed_sheets'])<a['n']+1]
    check(label+' all21 exclusions exactly bad inertia marks',len(reject)==21 and
      Counter((a['D_order'],a['I_order'],a['n']) for a in reject)=={(6,3,2):10,(10,5,2):6,(12,4,3):5})
    check(label+' allowed actual fiber counts all rows',all(
       a['allowed_actual_k']==expected_k[a['D_order'],a['I_order']] for a in accepted))
    check(label+' V4 inertia has four omitted two-sheet local factors',all(
       orbit_lengths(frozenset(a['I']),perms)==(1,1,1,2,2,2,2)
       for a in accepted if a['D_order']==4))
    check(label+' finite-near-axis degree11 only n5',[(a['n'],a['D_order'],a['I_order']) for a in accepted
       if 11 in a['allowed_actual_k']]==[(5,5,1)]*6)
    groupmodels.append({'action':j,'subgroups':[
      {'elements':sorted(D),'generators':list(gens)} for D,gens in sorted(groups.items(),key=lambda v:(len(v[0]),tuple(sorted(v[0]))))],
      'all_cyclic_quotient_pairs':pairs,'accepted_pairs':accepted,'accepted_pair_count':len(accepted)})
manifest=json.loads((ROOT/'sources/homogeneous-axis-manifest.json').read_text())
for src in manifest['sources']:
    raw=(ROOT/src['path']).read_bytes()
    check('primary bytes '+src['path'],len(raw)==src['bytes'])
    check('primary hash '+src['path'],hashlib.sha256(raw).hexdigest()==src['sha256'])
model={'claims':['KMFI118','KMFI119','KMFI120'],
 'homogeneous_matrix_degree':'n>=1','all_actual_axis_fibers':'finite reduced etale; k=1+n*m',
 'nonzero_axis_points':'f=(n+1)A,gamma=-nA; rankP<=r-2 and rankDX=rankP+1 forr>=2',
 'all_r_ge2_homogeneous_sweeps_globally_nonfinite':True,
 'rank1_finiteness_excluded':False,
 'hypothetical_finite_residue_controls':residues,
 'positive_axis_orbit_controls':positive,'rational_full_rank_norm_controls':normmodels,
 'PSL11_complete_positive_axis_rows':[
 {'n':2,'D':'C2','I':'1','actual_k':[3,5,7,9]},
 {'n':2,'D':'V4','I':'C2','actual_k':[3],'omitted_geometric_normalization_points':4,'omitted_local_generic_rank':8},
 {'n':3,'D':'C3','I':'1','actual_k':[4,7,10]},
 {'n':5,'D':'C5','I':'1','actual_k':[6,11]}],
 'full_A5_normal_quotient_transcripts':groupmodels,
 'finite_axis_criterion':'k=N by local generic ranks; no flatness assumption on normalization',
 'PSL11_finite_near_axis_required_degree':5,
 'PSL11_positive_axis_orbits_forced_to_exist':False,
 'PSL11_other_homogeneous_degrees':'only distinguished axis point, omitted local generic rank10',
 'local_generic_rank_deficit_promoted_to_point_or_divisor_count':False,
 'rows_claimed_as_actual_PSL_realizations':False,
 'actual_PSL11_input_supplied':False,'absolute_exceptional11_action_decided':False,
 'full_action_spectrum_complete':False,'goal_complete':False}
(ROOT/'homogeneous_axis_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['HOMOGENEOUS_AXIS_OBSTRUCTION.md','homogeneous_axis.py','verify_homogeneous_axis.py','homogeneous_axis_model.json',
 'sources/homogeneous-axis-manifest.json','sweep_differential.py','SWEEP_DIFFERENTIAL_RECONSTRUCTION.md',
 'local_image.py','psl11_residues.py','psl11_curve_spectrum.py','HIGHER_RANK_COMPLETION.md']
paths += [a['path'] for a in manifest['sources']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-started,3),'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 'assurance':'Written all-degree homogeneous fiber and nonfiniteness proofs, torus residue and exhaustive A5 quotient argument, local-rank finite-axis criterion and rational all-rank full-rank cyclic norm control. Exact checks supplement proofs; no formal, independent whole-proof or external review claimed.',
 'goal_complete':False}
(ROOT/'homogeneous_axis_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'runtime_seconds':receipt['runtime_seconds'],
 'A5_subgroups_per_action':59,'cyclic_pairs_per_action':67,'accepted_pairs_per_action':46,
 'finite_near_axis_PSL11_matrix_degree':5,'actual_PSL11_input_supplied':False,'goal_complete':False}))
