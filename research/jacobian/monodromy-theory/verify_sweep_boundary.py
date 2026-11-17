"""Exact controls for KMFI100--102; universal compactification proof is separate."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
from collections import Counter
import json,hashlib,time,platform
import sympy as S
from sweep_boundary import line_plucker,ruled_data,nonfinite_control,binary_control
from polynomial_wave import polynomial_chart
from psl11_residues import exceptional_actions
from central_lift import sl2_prime,permutation_order
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,val):
    if not val:raise AssertionError(name)
    checks.append(name)
def zero(x):return S.cancel(x)==0
w,z,g,t,v,k=S.symbols('w z gamma t v kappa')
# The compactification retains direction points, rather than just source poles.
x1,x2,x3=S.symbols('X1 X2 X3');X=S.Matrix([x1,x2,x3])
pl=line_plucker(X,w,z)
check('all six line Plucker coordinates',pl=={(0,1):1,(0,2):w,(0,3):z,
 (1,2):x1*w-x2,(1,3):x1*z-x3,(2,3):x2*z-x3*w})
check('Plucker quadric',zero(pl[0,1]*pl[2,3]-pl[0,2]*pl[1,3]+pl[0,3]*pl[1,2]))
q0,q1,q2,q3,a,b,c,d=S.symbols('q0 q1 q2 q3 a b c d')
point=S.Matrix([q0,q1,q2,q3]);dir_inf=S.Matrix([0,0,a,b])
linepoint=point*c+dir_inf*d
check('boundary affine line is horizontal',S.diff(linepoint[1]/linepoint[0],d)==0)
check('boundary marked direction has first two coordinates zero',list(dir_inf[:2])==[0,0])
# Full differential identity for arbitrary rational ruled coefficients.
ff=S.Function('f')(t);rr=S.Function('r')(t);ss=S.Function('s')(t)
rd=ruled_data(ff,rr,ss,t,v,w);normal=rd['normal'];surface=rd['X']
check('ruled normal tangent t',zero(normal.dot(surface.diff(t))))
check('ruled normal tangent v',zero(normal.dot(surface.diff(v))))
II=S.Matrix([[normal.dot(surface.diff(t,2)),normal.dot(surface.diff(t,v))],
             [normal.dot(surface.diff(v,t)),normal.dot(surface.diff(v,2))]])
vec=S.Matrix([1/S.diff(ff,t),w])
check('full asymptotic differential identity',
 zero((vec.T*II*vec)[0]*S.diff(ff,t)**3-rd['asymptotic_numerator']))
check('horizontal ruling is asymptotic',II[1,1]==0)
vv=-(rd['B']+2*S.diff(rr,t)*S.diff(ff,t)**2*w)/rd['A']
ztan=rr*w+(S.diff(rr,t)*vv+S.diff(ss,t))/S.diff(ff,t)
check('complete v elimination',zero(ztan-rd['p']*w-rd['q']))
check('A zero is constant derivative ratio',zero(S.diff(S.diff(rr,t)/S.diff(ff,t),t)
   -rd['A']/S.diff(ff,t)**2))
# Affine covariance of the normal direction and its one-variable property.
rho=S.symbols('rho');T=S.Matrix([[a,b],[c,d]])
normal2=T.T*S.Matrix([-rho,1])
check('normal slope affine Mobius transform',
 zero(-normal2[0]/normal2[1]-(a*rho-c)/(d-b*rho)))
check('nonconstant slope stays nonconstant',
 zero(S.diff((a*rho-c)/(d-b*rho),rho)-(a*d-b*c)/(d-b*rho)**2))
P=S.Matrix([[-rho,1],[-rho**2,rho]])
check('nilpotent affine covariance',(T.inv()*P*T)**2==S.zeros(2))
# Polynomial inverse to a generic degree-one rational pencil.
p0,p1,q0,q1,r0,r1,Z,W,U=S.symbols('p0 p1 q0 q1 r0 r1 Z W U')
den=p1*W+q1-r1*Z;num=r0*Z-p0*W-q0
inv=num/den
check('Mobius-pencil exact inverse',
 zero((Z*(r1*U+r0)-W*(p1*U+p0)-(q1*U+q0)).subs(U,inv)))
check('nonconstant denominator coefficients exactly p1 and minus r1',
 S.Poly(den,W,Z).coeff_monomial(W)==p1 and S.Poly(den,W,Z).coeff_monomial(Z)==-r1)
check('constant-denominator inverse is affine',
 S.Poly(S.expand(num.subs({p1:0,r1:0})),W,Z).total_degree()==1)
# A=0 horizontal ruling, including two geometric lines above a generic Y1.
simple=polynomial_chart(w,[(0,1)],w,z)
check('simple horizontal control exact surface',simple['X']==S.Matrix([w,z,w*z-w**3/3]))
for h in [w,w**2,w**3+w+1]:
    f=polynomial_chart(h,[(0,1)],w,z)
    Xf=f['X'];rt=h;st=S.cancel((Xf[2]-h*Xf[1]).subs(w,t))
    data=ruled_data(h.subs(w,t),h.subs(w,t),st,t,v,w)
    check('A0 ruled graph h='+str(h),data['A']==0 and zero(data['asymptotic_numerator'].subs(w,t)))
    check('A0 polynomial normalization h='+str(h),
      zero((Xf+g*S.Matrix([1,w,z])).jacobian([g,w,z]).det()-g*g))
    check('finite chart control h='+str(h),f['N'] in [3,5,7] and f['delta']==1)
# A!=0: swap the same finite surface. The multiplier is initially1/z.
swapped=S.Matrix([z,z*w-z**3/3,w])
rd2=ruled_data(t,1/t,t*t/3,t,v,w)
check('A nonzero swap control',zero(rd2['A']-2/t**3))
check('A nonzero elimination returns affine t=z',rd2['p']==0 and zero(rd2['q']-t))
check('swapped surface graph',zero(swapped[2]-swapped[1]/z-z*z/3))
NN=S.Matrix([[-w,1,0],[-z,0,1]]);PP=NN*swapped.jacobian([w,z])
check('swapped normalization rational one-coordinate slope',zero(PP[1,1]/PP[0,1]-1/z))
check('swapped sweep determinant',zero((swapped+g*S.Matrix([1,w,z])).jacobian([g,w,z]).det()-g*g))
# A genuinely nonfinite sweep has a distinct omitted plane, not a missing D prime.
minimal=k*k+k+1
def red(expr):
    num,den=S.fraction(S.cancel(expr))
    return S.rem(S.rem(num,minimal,k)*S.invert(S.rem(den,minimal,k),minimal,k),minimal,k).expand()
nf=nonfinite_control(w,z,g,k);a0=nf['a'];XX=nf['X'];YY=nf['S']
check('nonfinite control normalized determinant',red(YY.jacobian([g,w,z]).det()-g*g)==0)
check('nonfinite control collapsed line',all(e.subs({w:0,g:0})==0 for e in YY))
check('nonfinite control rank two',red(XX.diff(w).cross(XX.diff(z))[2])!=0)
Y1,Y2,Y3=S.symbols('Y1 Y2 Y3');targetD=Y2*Y3-Y1**3/(3*a0)
check('nonfinite exact ramification surface',red(targetD.subs({Y1:XX[0],Y2:XX[1],Y3:XX[2]},simultaneous=True))==0)
check('nonfinite complete ramification pullback',
 red(targetD.subs({Y1:YY[0],Y2:YY[1],Y3:YY[2]},simultaneous=True)+g**3/(3*a0))==0)
check('nonfinite D not horizontal-line ruled',
 set(S.Poly(targetD.subs(Y1,1),Y2,Y3).monoms())=={(1,1),(0,0)} and
 red(targetD.subs({Y1:1,Y2:0,Y3:0}))!=0)
check('nonfinite control birational branch inverse',
 red(a0*XX[1]/((a0+1)*XX[0])-w)==0 and red(-a0*XX[2]/XX[0]-z)==0)
check('nonfinite inverse root fields',
 red(YY[1]/(g+(a0+1)*(YY[0]-g)/a0)-w)==0 and
 red(YY[2]/(g-(YY[0]-g)/a0)-z)==0)
ep,V0,BB=S.symbols('epsilon V0 BB',nonzero=True)
arc={w:ep,z:V0/ep,g:V0+ep*BB/V0}
arcY=YY.subs(arc,simultaneous=True)
check('nonfinite exact escaping arc',all(red(a-b)==0 for a,b in zip(arcY,
 [(a0+1)*V0+ep*BB/V0,(a0+2)*V0*ep+ep**2*BB/V0,BB])))
check('nonfinite boundary target plane and horizontal lines',
 all(red(S.limit(val,ep,0)-want)==0 for val,want in zip(arcY,[(a0+1)*V0,0,BB])))
check('nonfinite boundary plane differs from D',red(targetD.subs(Y2,0))!=0)
# Sharp delta2 finite binary control.
binary=binary_control(w,z,g);XB,YB=binary['X'],binary['S']
check('sharp binary determinant',zero(YB.jacobian([g,w,z]).det()-g*g))
disc=9*(Y3+Y1**2/2)**2-8*Y2**3
check('sharp binary branch pullback',
 zero(disc.subs({Y1:YB[0],Y2:YB[1],Y3:YB[2]},simultaneous=True)-g**3*(9*g+4*w**3)/4))
check('sharp binary generic two-to-one branch',
 all(zero(e) for e in XB-XB.subs({w:-w,z:z-2*w**3},simultaneous=True)))
check('sharp binary unramified tail',S.degree((9*g+4*w**3)/4,g)==1)
# The actual exceptional actions have no possible complete visible inertia type.
data,orders,triangles,classes,actions=exceptional_actions()
def cycle_lengths(p):
    seen=set();out=[]
    for i in range(len(p)):
        if i in seen:continue
        j=i;n=0
        while j not in seen:seen.add(j);n+=1;j=p[j]
        out.append(n)
    return sorted(out)
for j,action in enumerate(actions):
    profiles={tuple(cycle_lengths(action['permutations'][g0])) for g0,o in orders.items() if o==3}
    check('all exceptional order3 profiles action'+str(j),profiles=={(1,1,3,3,3)})
    check('no complete sweep inertia action'+str(j),not profiles&{(1,)*8+(3,), (1,)*5+(3,3)})
# Small projective controls support the universal order3 fixed-point proof.
natural=[]
for q in [3,5,7,11,13]:
    dg=sl2_prime(q);pp=[dg['act'](mat) for mat in dg['quotient']]
    profiles={tuple(cycle_lengths(perm)) for perm in pp if permutation_order(perm)==3}
    fixed={profile.count(1) for profile in profiles}
    check('natural PSL2 '+str(q)+' order3 fixed sheets',fixed=={1 if q==3 else 2 if q%3==1 else 0})
    natural.append({'q':q,'order3_profiles':[list(a) for a in sorted(profiles)]})
# Positive small-q field witnesses.
rankone=S.Matrix([z*z,0,2*z**3/3]);sr=rankone+g*S.Matrix([1,w,z])
check('q2 rank-one determinant',zero(sr.jacobian([g,w,z]).det()-g*g))
check('q2 exact cubic pencil',zero(sr[2]-sr[0]*z+z**3/3))
a5=polynomial_chart(w*w,[(0,1)],w,z)
check('q4 natural A5 finite witness',a5['N']==5 and a5['geometric_group']=='A5')
manifest=json.loads((ROOT/'sources/sweep-boundary-manifest.json').read_text())
for source in manifest['sources']:
    raw=(ROOT/source['path']).read_bytes()
    check('primary bytes '+source['path'],len(raw)==source['bytes'])
    check('primary hash '+source['path'],hashlib.sha256(raw).hexdigest()==source['sha256'])
model={'claims':['KMFI100','KMFI101','KMFI102'],
 'horizontal_ruling_implies_finite_polynomial_chart':True,
 'all_boundary_graph_curves_included':True,'generic_D_normalization_complete':True,
 'global_finiteness_claimed_for_every_sweep':False,
 'full_ramification_surface_inertia':['3*1^(N-3)','3^2*1^(N-6)'],
 'both_exceptional_PSL11_sweep_field_actions_excluded':True,
 'all_primitive_pair_and_scalar_degrees_covered':True,
 'natural_projective_PSL2_sweep_occurrence_q':[2,4],
 'small_natural_group_controls':natural,
 'nonfinite_control':{'coefficient_equation':'kappa^2+kappa+1=0',
 'degree':3,'delta':1,'ramification_surface':str(targetD),
 'complete_pullback':'-gamma^3/(3*(-kappa-2))',
 'source_line_collapsed':True,'globally_finite':False,'Keller_completion':False},
 'sharp_delta2_control':'binary SL2(F3) degree8 full inertia3^2*1^2',
 'absolute_exceptional_Keller_action_decided':False,'full_action_spectrum_complete':False,'goal_complete':False}
(ROOT/'sweep_boundary_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['SWEEP_BOUNDARY_COMPLETENESS.md','sweep_boundary.py','verify_sweep_boundary.py',
 'sweep_boundary_model.json','sources/sweep-boundary-manifest.json',
 'SWEEP_ASYMPTOTICS.md','POLYNOMIAL_MULTIPLIER_CLASSIFICATION.md',
 'BRANCH_CENTERED_COMPLETION.md','polynomial_wave.py','psl11_residues.py','central_lift.py']
paths += [s['path'] for s in manifest['sources']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-started,3),'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{name:hashlib.sha256((ROOT/name).read_bytes()).hexdigest() for name in paths},
 'assurance':'Written universal horizontal-ruling and line-family compactification proof; symbolic identities and finite-group controls are supplemental, not formal or independent whole-proof verification.',
 'goal_complete':False}
(ROOT/'sweep_boundary_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],
 'all_normalized_sweep_PSL11_routes_excluded':True,'absolute_exceptional_action_decided':False,'goal_complete':False}))
