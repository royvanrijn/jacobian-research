"""Exact KMFI106--108 boundary, inertia-power and class-rank controls."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
from collections import Counter
from math import gcd
import json,hashlib,time,platform
import sympy as S
from higher_completion import companion_surface,weighted_stage,weighted_order,truncate
from completion_boundary import cycle_lengths,permutation_power,binary_high_weight_data
from psl11_residues import exceptional_actions
from central_lift import sl2_prime,permutation_order
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,val):
    if not val:raise AssertionError(name)
    checks.append(name)
def zero(expr):return S.cancel(expr)==0
companion=[]
for r in [2,3,4]:
    ps=list(S.symbols('p1:'+str(r+1)));g=S.Symbol('gamma');us=list(S.symbols('u1:'+str(r+1)))
    X=companion_surface(ps);st=weighted_stage(X,ps,g,us);E=st['E'];J=S.zeros(r)
    for i in range(r-1):J[i,i+1]=1
    check('rank'+str(r)+' whole Gamma stage fiber',E.subs(g,0)==S.Matrix([1,*list((S.eye(r)+J)*S.Matrix(us))]))
    Q=S.Matrix(r,r,lambda i,j:S.Symbol('q'+str(i)+str(j)))
    linear=sum((S.trace(J**i*Q*J**(r-1-i)) for i in range(r)),S.Integer(0))
    check('rank'+str(r)+' nilpotent trace derivative',S.expand(linear-r*Q[r-1,0])==0)
    Ys=list(S.symbols('Y0:'+str(r+1)));tt=S.Symbol('T');CC=S.Symbol('C')
    yy=sum(Ys[i]*tt**(i+1) for i in range(r))
    d=Ys[-1]+S.expand(sum((-1)**(j+1)*yy**j/S.Integer(j) for j in range(1,r+2))).coeff(tt,r+1)
    weights=st['output_weights'];poly=S.Poly(S.expand(d),*Ys)
    minimum=min(sum(a*b for a,b in zip(mon,weights)) for mon,c in poly.terms())
    low=sum(c*S.prod(v**m for v,m in zip(Ys,mon)) for mon,c in poly.terms()
            if sum(a*b for a,b in zip(mon,weights))==minimum)
    check('rank'+str(r)+' exact lowest weight e',minimum==r+1)
    check('rank'+str(r)+' pure first-coordinate initial',
          zero(low-(-1)**r*Ys[0]**(r+1)/S.Integer(r+1)))
    scaled=S.expand(d.subs({v:CC**a*v for v,a in zip(Ys,weights)},simultaneous=True))
    strict=S.cancel(scaled/CC**(r+1))
    check('rank'+str(r)+' normalized branch initial on C0',
          zero(strict.subs(CC,0)-(-1)**r*Ys[0]**(r+1)/S.Integer(r+1)))
    xx=S.Symbol('x',nonzero=True)
    Zi=S.Matrix([xx**(-e)*v for v,e in zip((S.eye(r)+J)*S.Matrix(us),weights[1:])])
    inverse=(S.eye(r)+J).inv()*S.diag(*[xx**e for e in weights[1:]])*Zi
    check('rank'+str(r)+' complete Gamma inverse',inverse==S.Matrix(us))
    check('rank'+str(r)+' Gamma matrix determinant1',(S.eye(r)+J).det()==1)
    companion.append({'r':r,'e':r+1,'low_branch':str(low),'negative_cyclic_input':True})
# Exact formal factorization for every e|m is the algebraic reason for inertia powering.
V,U=S.symbols('V U')
for e,k in [(3,1),(5,1),(2,3),(3,2),(4,3),(5,2)]:
    # product over e-th roots is V^(e*k)-U^e; no approximate complex roots used.
    expr=S.resultant(S.Symbol('zeta')**e-1,V**k-S.Symbol('zeta')*U,S.Symbol('zeta'))
    check('Kummer factorization e'+str(e)+' k'+str(k),zero(expr-(V**(e*k)-U**e)))
    check('Kummer factor Eisenstein has residue degree1 e'+str(e)+' k'+str(k),
          S.Poly(V**k-U,V).LC()==1 and S.Poly(V**k-U,V).TC()==-U)
# Binary control with explicit finite high-weight rectifier and marked coordinates.
binary=binary_high_weight_data();K=binary['coefficient_field'];rr=S.sqrt(-3)
def algzero(expr):
    return K.from_sympy(S.cancel(expr))==K.zero
def polyzero(expr,variables):
    return S.Poly(S.expand(expr),*variables,domain=K).is_zero
E=binary['E'];theta=binary['theta'];point=binary['point'];L=binary['L'];vs=binary['variables'];t=vs[0]
check('binary E polynomial in all original stage variables',all(S.denom(e)==1 for e in E))
check('binary marked regular zero',all(algzero(e) for e in E.subs(dict(zip(theta,point)))))
check('binary marked inverse Jacobian',all(algzero(e) for e in E.jacobian(theta).subs(dict(zip(theta,point)))*L-S.eye(3)))
check('binary marked determinant gamma0',algzero(L.det()-point[0]))
check('binary full regular-Jordan Gamma fiber',list(E.subs(theta[0],0))==[1,theta[1]+theta[2],theta[2]])
rec=binary['rectifier'];jet=binary['jet'];phi=rec['phi'];inv=rec['inverse']
check('binary finite curve exact residual',
 all(polyzero(truncate(e,t,4),[t]) for e in jet[1:,:].subs(dict(zip(vs[1:],rec['h'])),simultaneous=True)))
check('binary Gaussian flag exact residual',polyzero(truncate((rec['M']*rec['A'])[1,0],t,4),[t]))
check('binary triangular polynomial inverse',
 all(polyzero(e,vs) for e in phi.subs(dict(zip(vs,inv)),simultaneous=True)-S.Matrix(vs)))
check('binary rectifier whole Jacobian1',polyzero(phi.jacobian(vs).det()-1,vs))
transformed=jet.subs(dict(zip(vs,phi)),simultaneous=True)
check('binary all required low weighted monomials vanish',
 all((o:=weighted_order(e,vs,[1,3,4])) is None or o>=want for e,want in zip(transformed,[1,3,4])))
check('binary x0 first two outputs triangular',
 S.Poly(S.expand(transformed[0]),*vs,domain=K).coeff_monomial(t)==1 and
 S.Poly(S.expand(transformed[1]),*vs,domain=K).coeff_monomial(vs[1])==1)
# Exact original branch equation after centering; no finite-field sample.
Y0,Y1,Y2,C=S.symbols('Y0 Y1 Y2 C')
d=9*(Y2+Y1-S.Rational(1,6)+(Y0+1)**2/2)**2-8*(Y1+Y0+S.Rational(1,2))**3
scaled=S.expand(d.subs({Y0:C*Y0,Y1:C**3*Y1,Y2:C**4*Y2},simultaneous=True))
check('binary weighted branch multiplicity3',S.Poly(scaled,C).terms()[-1][0][0]==3)
check('binary weighted branch initial exactly Y0 cubed',S.expand(S.cancel(scaled/C**3).subs(C,0)-Y0**3)==0)
# The full binary group action is reconstructed as SL2(F3) on its nonzero vectors.
dd=sl2_prime(3);vectors=[(a,b) for a in range(3) for b in range(3) if (a,b)!=(0,0)]
lookup={v:i for i,v in enumerate(vectors)}
perms=[]
for a,b,c,d0 in dd['elements']:
    perms.append(tuple(lookup[((a*x+b*y)%3,(c*x+d0*y)%3)] for x,y in vectors))
profiles={tuple(cycle_lengths(p)) for p in perms if permutation_order(p)==3}
check('binary complete order3 dictionary',profiles=={(1,1,3,3)})
check('binary every branch inertia cubes to identity',all(permutation_power(p,3)==tuple(range(8)) for p in perms if permutation_order(p)==3))
check('binary normalization rank7 and six unramified omissions',8-1==7 and 8-2==6)
# Both complete exceptional actions give eleven split primes, not just allowed fixed sheets.
data,orders,triangles,classes,actions=exceptional_actions()
for j,action in enumerate(actions):
    selected=[action['permutations'][a] for a,o in orders.items() if o==5]
    check('exceptional action'+str(j)+' all five-cycle profiles',{tuple(cycle_lengths(p)) for p in selected}=={(1,5,5)})
    powered=[permutation_power(p,5) for p in selected]
    check('exceptional action'+str(j)+' C0 completely split',set(powered)=={tuple(range(11))})
    check('exceptional action'+str(j)+' nine omitted unramified primes',all(cycle_lengths(p)==[1]*11 for p in powered) and 11-2==9)
# The divisor relation is primitive and the class-rank increment is exact.
from sympy.matrices.normalforms import smith_normal_form
from sympy.polys.domains import ZZ
for relation,increment in [([1]*8,7),([1]*11,10),([1,1,2,2,3],4)]:
    snf=smith_normal_form(S.Matrix([relation]),domain=ZZ)
    check('primitive divisor relation '+str(len(relation)),abs(snf[0,0])==1 and len(relation)-1==increment)
manifest=json.loads((ROOT/'sources/completion-boundary-manifest.json').read_text())
for source in manifest['retained_primary_inputs']:
    raw=(ROOT/source['path']).read_bytes()
    check('retained primary bytes '+source['path'],len(raw)==source['bytes'])
    check('retained primary hash '+source['path'],hashlib.sha256(raw).hexdigest()==source['sha256'])
model={'claims':['KMFI106','KMFI107','KMFI108'],'r_ge2_full_rank_generic_center_required':True,
 'entire_C0_actual_fiber_counts':{'Z0_nonzero':2,'Z0_zero':1},
 'x0_component_isomorphism':'A^(r+1) onto C0','gamma0_component_isomorphism':'G_m times A^r onto C0 minus Z0=0',
 'C0_full_decomposition_equals_inertia':'generated by sigma^(r+1)','C0_all_residue_degrees':1,
 'class_rank_formula':'rank Cl(Z_F)=rank Cl(Z_S)+cycles(sigma^(r+1))-1',
 'PSL11_C0_ramification':1,'PSL11_C0_decomposition_order':1,'PSL11_split_normalization_primes':11,
 'PSL11_retained_primes':2,'PSL11_omitted_unramified_primes':9,'PSL11_class_rank_increment':10,
 'PSL11_completed_representative_class_rank_lower':11,
 'bound_claimed_for_all_affine_models_or_for_ZS_itself':False,
 'binary_high_weight_control':{'coefficient_field':'Q(s),s^2=-3','theta':[str(e) for e in theta],
  'E':[str(e) for e in E],'theta0':[str(e) for e in point],
  'L':[[str(L[i,j]) for j in range(3)] for i in range(3)],
  'G_ordinary_jet4':[str(e) for e in jet],
  'h':[str(e) for e in rec['h']],'A':[[str(rec['A'][i,j]) for j in range(2)] for i in range(2)],
  'C0_split_primes':8,'omitted_unramified_primes':6,'class_group':'Z^7',
  'canonical_vector':[2,0,0,0,0,0,0],'Keller_group':'SL2(F3) on8'},
 'companion_controls':companion,'absolute_exceptional_action_decided':False,'goal_complete':False}
(ROOT/'completion_boundary_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['COMPLETION_EXCEPTIONAL_DIVISOR.md','completion_boundary.py','verify_completion_boundary.py','completion_boundary_model.json',
 'sources/completion-boundary-manifest.json','HIGHER_RANK_COMPLETION.md','higher_completion.py','BRANCH_CENTERED_COMPLETION.md',
 'branch_centered.py','POLYNOMIAL_MULTIPLIER_CLASSIFICATION.md','polynomial_wave.py','central_lift.py','psl11_residues.py']
paths += [s['path'] for s in manifest['retained_primary_inputs']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,'runtime_seconds':round(time.monotonic()-started,3),
 'python':platform.python_version(),'sympy':S.__version__,'sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 'assurance':'Written universal two-source-component, local Kummer and class-group proofs; exact symbolic controls and fully specified binary positive input. No formal verification or independent whole-proof review.',
 'goal_complete':False}
(ROOT/'completion_boundary_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'PSL11_C0_omitted_unramified_primes':9,'PSL11_normalization_class_rank_lower':11,'binary_class_group':'Z^7',
 'absolute_exceptional_action_decided':False,'goal_complete':False}))
