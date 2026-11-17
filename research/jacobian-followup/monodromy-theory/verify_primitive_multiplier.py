"""Exact KMFI124--126 cyclicity, degree-free multiplier and quadratic shear controls."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
import sympy as S
import json,time,hashlib,platform
from integrated_sweeps import family
from polynomial_wave import polynomial_chart
from sweep_differential import divergence_columns,gradient
from homogeneous_axis import norm_control
from primitive_multiplier import feedback_form,logarithmic_residues,quadratic_shear,direction_shape
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,value):
    if not value:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def mz(M):return all(zero(e) for e in M)
t=S.Symbol('t');controls=[]
for r in [2,3,4]:
 d=family(r,t*t,t);controls.append(('integrated_rank'+str(r),d['parameters'],d['X']))
v,w=S.symbols('v w')
controls.append(('binary_full_rank',[v,w],polynomial_chart(v*v,[(1,0)],v,w)['X']))
gauss=[]
for name,ps,X in controls:
 r=len(ps);P=(X[1:,:].jacobian(ps)-S.Matrix(ps)*gradient(X[0],ps)).applyfunc(S.expand)
 top=(P**(r-1)).applyfunc(S.expand)
 row=next(i for i in range(r) if any(top[i,j]!=0 for j in range(r)));left=top[row,:]
 normal=S.Matrix([-(left*S.Matrix(ps))[0],*list(left)])
 monoms=set()
 for e in normal:monoms.update(S.Poly(e,*ps).monoms())
 coeff=S.Matrix([[S.Poly(e,*ps).coeff_monomial(mon) for e in normal] for mon in sorted(monoms)])
 check(name+' whole Gauss normal has no constant relation',coeff.rank()==r+1)
 q0=S.Symbol('q0');qs=list(S.symbols('q1:'+str(r+1)));vec=S.Matrix(qs)-q0*S.Matrix(ps)
 for m in range(1,r):
  value=(P**m*vec).applyfunc(S.expand)
  actual=S.expand(sum(S.diff(value[j],ps[j]) for j in range(r)))
  target=-(r-m)*(gradient(X[0],ps)*P**(m-1)*vec)[0]
  check(name+' entire affine-Euler divergence identity power'+str(m),zero(actual-target))
 gauss.append({'control':name,'dual_Gauss_coefficient_rank':r+1,'source_rank':r})
# Cyclicity is generic, not a claim at every point.
d=family(4,t*t,t);ps=d['parameters'];P=(d['X'][1:,:].jacobian(ps)-S.Matrix(ps)*gradient(d['X'][0],ps)).applyfunc(S.expand)
point={p:(1 if i==0 else 0) for i,p in enumerate(ps)}
vec=S.Matrix(ps)-S.Matrix([1,0,0,0]);cyclic=(P**3*vec).applyfunc(S.expand)
check('affine-Euler cyclicity generic but vector vanishes at center',
 not mz(cyclic) and cyclic.subs(point)==S.zeros(4,1))
positive=[]
for r in [2,3,4]:
 d=norm_control(r);ps=d['parameters'];V=d['M']*S.Matrix(ps)
 shape=direction_shape(V,ps);f=feedback_form(V,ps,S.Integer(1));logs=logarithmic_residues(f,ps)
 check('norm rank'+str(r)+' Jacobian power shape',shape['power_compatible'])
 check('norm rank'+str(r)+' exact cyclic feedback solve',mz(f['eta']*f['C']-f['moments']))
 T=f['T'];A=f['A']
 characteristic=S.expand((T-1)**r)-f['chi']+(f['eta']*(T*S.eye(r)-A).adjugate()*V)[0]
 check('norm rank'+str(r)+' entire feedback characteristic identity',zero(characteristic))
 check('norm rank'+str(r)+' logarithmic coefficient matching accepted',logs['compatible'])
 check('norm rank'+str(r)+' all residues positive integers',
 all(v.is_Integer and v>0 for v in logs['residues']))
 check('norm rank'+str(r)+' entire polynomial multiplier exact',
 all(zero(f['eta'][i]-S.diff(d['F'],ps[i])/d['F']) for i in range(r)))
 check('norm rank'+str(r)+' Euler degree condition',
 zero((f['eta']*S.Matrix(ps))[0]-r) and zero(1/shape['delta']-1-r))
 at={p:(1 if i==0 else 0) for i,p in enumerate(ps)}
 N=f['N'].subs(at);rho=f['rho'].subs(at)
 O=S.Matrix.vstack(*[rho*N**j for j in range(r)])
 check('norm rank'+str(r)+' full-rank rational coframe',O.det()!=0)
 positive.append({'r':r,'det_C':str(f['detC']),'eta':[str(e) for e in f['eta']],
                  'logarithmic_factors':[str(e) for e in logs['factors']],
                  'residues':[str(e) for e in logs['residues']],'coframe_atunit':str(O.det()),
                  'actual_geometric_group':'regular C'+str(r+1),'new_Keller_action_claimed':False})
# The fourth-power condition itself is a real obstruction.
u,v,w,t,lam=S.symbols('u v w t lam');ps=[u,v,w,t]
square=direction_shape(S.Matrix([p*p for p in ps]),ps)
check('primitive quadratic square-map fails perfect-fourth-power gate',
 not square['power_compatible'] and square['determinant']!=0)
V=S.Matrix([u*u,3*u*v,4*u*w,5*u*t+v*w])
shape=direction_shape(V,ps);f=feedback_form(V,ps,lam*u);logs=logarithmic_residues(f,ps,lam)
check('quadratic shear passes determinant-power shape',shape['power_compatible'] and shape['h0']==u and shape['delta']==120)
check('whole degree4 parameter feedback identity',mz(f['eta']*f['C']-f['moments']))
T=f['T'];char=(T-lam*u)**4-f['chi']+(f['eta']*(T*S.eye(4)-f['A']).adjugate()*V)[0]
check('quadratic shear entire characteristic solve',zero(char))
check('isolated parameter not discarded by generic linear solving',
 zero(logs['parameter_condition'].as_expr()-(lam-5)**4))
check('all residue and consistency degrees at most4',
 all(S.degree(e,lam)<=4 for e in logs['residues']) and logs['parameter_condition'].degree()==4)
eta5=[S.factor(e.subs(lam,5)) for e in f['eta']]
check('exact closed lambda5 form',eta5==[S.Rational(29,3)/u,-S.Rational(4,3)/v,S.Rational(1,12)/w,0])
check('closedness alone cannot accept the multiplier',
 all(zero(S.diff(eta5[j],ps[i])-S.diff(eta5[i],ps[j])) for i in range(4) for j in range(i+1,4)) and
 any(not e.subs(lam,5).is_Integer or e.subs(lam,5)<0 for e in logs['residues']))
check('exact degree-free Euler constraint',
 zero((f['eta']*S.Matrix(ps))[0]-(lam**4/S.Integer(60)-2)))
negative={'V':[str(e) for e in V],'det_C':str(f['detC']),'eta':[str(e) for e in f['eta']],
 'logarithmic_factors':[str(e) for e in logs['factors']],'residue_polynomials':[str(e) for e in logs['residues']],
 'parameter_condition':str(logs['parameter_condition'].as_expr()),'only_candidate':5,'closed_candidate_form':[str(e) for e in eta5],
 'candidate_rejected_by_nonnegative_integer_residues':True}

# All-parameter identities, not a parameter sample proving the family theorem.
shear=quadratic_shear();pp=shear['parameters'];u,v,w,t=pp
a,b,c=shear['a'],shear['b'],shear['c'];V=shear['V'];A=shear['A'];Q=shear['Q'];KK=shear['K'];RR=shear['R']
C=S.Matrix.hstack(*[A**j*V for j in range(4)])
expected=-a*b*(a-1)*(a-b)*(b-1)*u**10*v*w*(KK*u*t+RR)
check('whole quadratic shear Krylov determinant all coefficients',zero(C.det(method='domain-ge')-expected))
adj=(c*u*S.eye(4)-A).adjugate()*V
check('whole quadratic shear c-eigenvalue adjugate all coefficients',
 all(zero(adj[i]-(u**3*(KK*u*t+RR) if i==3 else 0)) for i in range(4)))
A1=A.subs(u,1);V1=V.subs(u,1)
adja=(a*S.eye(4)-A1).adjugate()*V1
expecta=S.Matrix([0,a*v*(a-1)*(a-b)*(a-c),0,a*v*(a-1)*(a-b)*S.diff(Q,v)])
adjb=(b*S.eye(4)-A1).adjugate()*V1
expectb=S.Matrix([0,0,b*w*(b-1)*(b-a)*(b-c),b*w*(b-1)*(b-a)*S.diff(Q,w)])
check('whole a-eigenvalue residue identity including resonance2',mz(adja-expecta))
check('whole b-eigenvalue residue identity including resonance2',mz(adjb-expectb))
alpha,beta,gamma=shear['q']
Sq=alpha*v*v/(c-a)+gamma*w*w/(c-b)
trace_res=(alpha*v*v+gamma*w*w)+(a-1)*v*S.diff(Sq,v)+(b-1)*w*S.diff(Sq,w)-(c-1)*Sq
check('whole final Darboux pole contradiction',
 zero(trace_res-alpha*(a-1)*v*v/(c-a)-gamma*(b-1)*w*w/(c-b)))
check('mixed derivative forces cross coefficient zero',
 zero(S.diff(Q,v,w)*(1/(c-a)-1/(c-b))-beta*(a-b)/((c-a)*(c-b))))
# Fixed resonant controls test that the implementation never divides by a-2.
resonant=[]
for aval,bval,cval in [(2,3,4),(3,4,2),(3,4,3),(3,4,1)]:
    ps=pp;vv=V.subs({a:aval,b:bval,c:cval,alpha:0,beta:1,gamma:0})
    rdata=feedback_form(vv,ps,lam*u);rlogs=logarithmic_residues(rdata,ps,lam)
    check('resonant shear'+str((aval,bval,cval))+' isolated condition retained',
       zero(rlogs['parameter_condition'].as_expr()-(lam-cval)**4))
    check('resonant shear'+str((aval,bval,cval))+' candidate independent of t',
       zero(rdata['eta'][3].subs(lam,cval)))
    resonant.append({'weights':[aval,bval,cval],'condition':str(rlogs['parameter_condition'].as_expr())})
# Cyclicity must reject these before multiplier reconstruction.
for aval,bval,cval in [(1,4,5),(3,3,5)]:
    vv=V.subs({a:aval,b:bval,c:cval,alpha:0,beta:1,gamma:0})
    try:feedback_form(vv,pp,lam*u)
    except ValueError:checks.append('singular primitive Krylov direction rejected '+str((aval,bval,cval)))
    else:raise AssertionError('singular cyclic direction accepted')
# A real polynomial multiplier can pass the log test but still have rank3.
vv=V.subs({a:1+S.I,b:1-S.I,c:2,alpha:0,beta:1,gamma:0});ps=pp
sharp=feedback_form(vv,ps,2*u);ll=logarithmic_residues(sharp,ps)
check('sharp quadratic direction has polynomial multiplier',
 ll['compatible'] and mz(sharp['eta']-S.Matrix([[0,1/v,1/w,0]])) and
 all(co.is_Integer and co>=0 for co in ll['residues']))
F=v*w;X=S.Matrix([-2*u*F,*list(F*(vv-2*u*S.Matrix(ps)))])
P=(X[1:,:].jacobian(ps)-S.Matrix(ps)*S.Matrix([[S.diff(X[0],p) for p in ps]])).applyfunc(S.expand)
check('sharp quadratic normalized whole matrix nilpotent',mz((P**4).applyfunc(S.expand)))
gam=S.Symbol('gam')
check('sharp quadratic whole normalized sweep Jacobian',
 zero((X+gam*S.Matrix([1,*ps])).jacobian([gam,*ps]).det(method='domain-ge')-gam**4))
at={p:1 for p in ps};O=S.Matrix.vstack(*[sharp['rho']*sharp['N']**j for j in range(4)])
check('sharp log-positive example fails full-rank test',
 X.jacobian(ps).subs(at).rank()==3 and S.simplify(O.subs(at).det())==0 and
 all(zero(S.diff(e,t)) for e in X))
model={'claims':['KMFI124','KMFI125','KMFI126'],
 'universal_generic_cyclic_vectors':'every nonzero q_vec-q0*p',
 'simultaneous_pointwise_cyclicity_claimed':False,
 'Gauss_image_spans_full_dual_projective_space':True,
 'translation_cylinders_for_full_rank_branch_image':'none',
 'varying_horizontal_rulings_excluded':False,
 'Gauss_controls':gauss,
 'complete_fixed_primitive_direction_criterion':{
 'necessary_shape':'detDV=delta*h0^r nonzero; h=lambda*h0; C=(V,DV*V,...,(DV)^(r-1)*V) invertible',
 'eta':'first r moments of 1-(T-h)^r/chi_DV times C^-1',
 'polynomial_multiplier':'eta=sum nonnegative integer b_j*dlog(D_j), D_j factors ofdetC',
 'full_rank':'rho=dh+h*eta, N=DV+V*eta-hI; cyclic rho coframe invertible',
 'parameter_compatibility_degree_upper':'r',
 'all_lambda_values_claimed_finite':False,
 'generic_group_determined_by_passing_criterion':False},
 'accepted_norm_multipliers':positive,'nonlogarithmic_quadratic_countercontrol':negative,
 'quadratic_shear_entire_parameter_identities':{'V':[str(e) for e in V],'Q':str(Q),'K':str(KK),'R':str(RR),'det_C':str(expected)},
 'resonant_numeric_controls':resonant,
 'sharp_polynomial_multiplier_rank_deficient':{'weights':['1+i','1-i','2'],'Q':'v*w','F':'v*w','h':'2*u','X':[str(e) for e in X],'rankDX':3},
 'all_quadratic_shear_directions_closed_for_all_multiplier_degrees':True,
 'all_primitive_quadratic_directions_classified':False,
 'actual_PSL11_input_supplied':False,'absolute_exceptional11_action_decided':False,
 'full_action_spectrum_complete':False,'goal_complete':False}
(ROOT/'primitive_multiplier_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['PRIMITIVE_MULTIPLIER_CRITERION.md','primitive_multiplier.py','verify_primitive_multiplier.py','primitive_multiplier_model.json',
 'SWEEP_DIFFERENTIAL_RECONSTRUCTION.md','sweep_differential.py','integrated_sweeps.py',
 'polynomial_wave.py','homogeneous_axis.py','FACTORIZED_HOMOGENEOUS_FIELDS.md']
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-started,3),'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 'assurance':'Unbounded Euler/Gauss cyclicity, complete scalar multiplier criterion for a fixed primitive direction, and all-parameter quadratic shear full-rank obstruction. Exact whole identities and positive/negative/resonant controls supplement written proofs; no formal, independent whole-proof or external review claimed.',
 'goal_complete':False}
(ROOT/'primitive_multiplier_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'runtime_seconds':receipt['runtime_seconds'],
 'parameter_degree_bound':'r','all_quadratic_shear_multiplier_degrees_closed':True,
 'actual_PSL11_input_supplied':False,'goal_complete':False}))
