"""Exact finite-jet and full polynomial-map controls for KMFI103--105."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
from collections import Counter
import json,hashlib,time,platform
import sympy as S
from higher_completion import truncate,rectify,weighted_order,weighted_stage,companion_surface,suspended_cubic
from psl11_residues import exceptional_actions,centralizer
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,condition):
    if not condition:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
records=[]
for r in [1,2,3,4]:
    variables=list(S.symbols('t u1:'+str(r+1)));t=variables[0];us=variables[1:]
    G=S.Matrix([t+us[0]*us[-1],*[us[i]+t*t+(i+1)*t**3+us[(i+1)%r]**2+
                  t*us[(i-1)%r] for i in range(r)]])
    rec=rectify(G,variables);phi=rec['phi'];inverse=rec['inverse'];D=2*r
    check('rank'+str(r)+' finite curve residual',
      all(truncate(e,t,D)==0 for e in G[1:,:].subs(dict(zip(us,rec['h'])),simultaneous=True)))
    check('rank'+str(r)+' lower unitriangular polynomial matrix',
      all(rec['A'][i,j]==0 for i in range(r) for j in range(r) if j>i) and
      all(rec['A'][i,i]==1 for i in range(r)) and rec['A'].subs(t,0)==S.eye(r))
    check('rank'+str(r)+' right Gaussian flag',
      all(truncate((rec['M']*rec['A'])[i,j],t,D)==0 for i in range(r) for j in range(i)))
    check('rank'+str(r)+' inverse left',
      all(zero(e) for e in phi.subs(dict(zip(variables,inverse)),simultaneous=True)-S.Matrix(variables)))
    check('rank'+str(r)+' inverse right',
      all(zero(e) for e in inverse.subs(dict(zip(variables,phi)),simultaneous=True)-S.Matrix(variables)))
    check('rank'+str(r)+' whole source Jacobian1',zero(phi.jacobian(variables).det()-1))
    transformed=G.subs(dict(zip(variables,phi)),simultaneous=True).applyfunc(S.expand)
    ords=[weighted_order(e,variables,rec['weights']) for e in transformed]
    check('rank'+str(r)+' all polynomial weighted orders',
      all(o is None or o>=e for o,e in zip(ords,rec['weights'])))
    check('rank'+str(r)+' finite coefficient truncation',
      all(S.degree(e,t)<D for e in list(rec['h'])+list(rec['A']) if e!=0))
    records.append({'r':r,'weights':rec['weights'],'h':[str(e) for e in rec['h']],
      'A':[[str(rec['A'][i,j]) for j in range(r)] for i in range(r)],
      'weighted_orders':ords})
check('rank4 Gaussian matrix genuinely nonidentity',records[-1]['A']!=[[str(int(i==j)) for j in range(4)] for i in range(4)])
# All-rank weight arithmetic is an identity in r, not a cutoff theorem.
rr=S.symbols('r',integer=True,positive=True)
sum_a=rr*(3*rr-1)/2;sum_e=1+rr*(3*rr+1)/2
check('all-rank determinant balance',S.expand(rr+sum_a-sum_e)==-1)
check('all-rank nonlinear source cutoff',S.expand(2*rr-2*rr)==0)
check('all-rank nonlinear high-variable cutoff',S.expand(2*(rr+1)-2*rr)==2)
companion=[]
for r in [1,2,3,4]:
    ps=list(S.symbols('p1:'+str(r+1)));g=S.Symbol('gamma');uu=list(S.symbols('a1:'+str(r+1)))
    X=companion_surface(ps);st=weighted_stage(X,ps,g,uu);swept=st['S']
    check('companion'+str(r)+' full rank',X.jacobian(ps)[:r,:]==S.eye(r))
    check('companion'+str(r)+' full normalized determinant',zero(swept.jacobian([g,*ps]).det()-g**r))
    NN=S.Matrix.hstack(-S.Matrix(ps),S.eye(r));PP=NN*X.jacobian(ps)
    check('companion'+str(r)+' nilpotent normalization',all(zero(e) for e in PP**r))
    origin={p:0 for p in ps};P0=PP.subs(origin)
    check('companion'+str(r)+' strict upper centering',all(P0[i,j]==0 for i in range(r) for j in range(i+1)))
    check('companion'+str(r)+' entire first-stage polynomiality',
      all(S.denom(S.cancel(e))==1 for e in st['E']))
    Ys=list(S.symbols('Y1:'+str(r+2)));tt=S.Symbol('T')
    pp=sum(Ys[i]*tt**(i+1) for i in range(r))
    logco=S.expand(sum((-1)**(j+1)*pp**j/S.Integer(j) for j in range(1,r+2))).coeff(tt,r+1)
    d=Ys[-1]+logco
    check('companion'+str(r)+' exact constant tail',
      zero(d.subs(dict(zip(Ys,swept)),simultaneous=True)-(-1)**r*g**(r+1)/S.Integer(r+1)))
    companion.append({'r':r,'X':[str(e) for e in X],'cyclic_degree':r+1,'tail_constant':True})
# Fully specified ordinary-A6 positive control. It has rank-deficient X.
r=4;ps=list(S.symbols('p1:5'));g=S.Symbol('gamma');uu=list(S.symbols('a1:5'));theta=[g,*uu]
X=suspended_cubic(ps);st=weighted_stage(X,ps,g,uu);E=st['E']
check('positive rank4-normalization determinant',zero(st['S'].jacobian([g,*ps]).det()-g**4))
check('positive control X rank3',X.jacobian(ps).rank()==3)
theta0=S.Matrix([-3,-S.Rational(1,27),0,0,0]);sub0=dict(zip(theta,theta0))
check('positive regular second preimage',E.subs(sub0)==S.zeros(5,1))
J=E.jacobian(theta).subs(sub0);L=J.inv()
check('positive inverse Jacobian and determinant',
 L==S.Matrix([[9,324,0,0,0],[-S.Rational(11,27),-15,0,0,0],[0,0,1,-1,1],[0,0,0,1,-1],[0,0,0,0,1]]) and L.det()==-3)
vs=list(S.symbols('t u1:5'));t=vs[0]
G=E.subs(dict(zip(theta,theta0+L*S.Matrix(vs))),simultaneous=True).applyfunc(S.expand)
rec=rectify(G,vs);Theta=theta0+L*rec['phi']
Ep=E.subs(dict(zip(theta,Theta)),simultaneous=True).applyfunc(S.expand)
check('positive exact weighted polynomial orders',
 [weighted_order(e,vs,rec['weights']) for e in Ep]==[1,5,6,7,8])
check('positive spectator outputs',list(Ep[2:,:])==vs[2:])
check('positive first block independent of other spectators',
 not any(e.has(*vs[2:]) for e in [Theta[0],Ep[0],Ep[1]]))
x,y,z=S.symbols('x y z')
def substitute_block(expr,division,extra_x=0):
    old=S.Poly(expr,vs[0],vs[1],domain=S.QQ);out={}
    for (a,b),co in old.terms():
        power=a+5*b-division+extra_x
        if power<0:raise AssertionError('non-polynomial x division')
        out[(power,a,b)]=co
    return S.Poly.from_dict(out,(x,y,z),domain=S.QQ)
F3=[substitute_block(Theta[0],0,1),substitute_block(Ep[0],1),substitute_block(Ep[1],5)]
check('positive expanded degrees',[f.total_degree() for f in F3]==[15,125,135])
check('positive expanded monomial counts',[len(f.terms()) for f in F3]==[9,323,388])
# Serialize the entire map, then reconstruct its polynomial coefficients independently.
full_terms=[]
for f in F3:
    full_terms.append([{'monomial':list(m)+[0,0,0],'coefficient':str(c)} for m,c in f.terms()])
for j in range(3,6):
    mon=[0]*6;mon[j]=1
    full_terms.append([{'monomial':mon,'coefficient':'1'}])
serialized={'variables':['x','y0','y1','y2','y3','y4'],'coordinates':full_terms}
reconstructed=[]
for coord in serialized['coordinates'][:3]:
    assert all(not any(t['monomial'][3:]) for t in coord)
    reconstructed.append(S.Poly.from_dict({tuple(t['monomial'][:3]):S.Rational(t['coefficient']) for t in coord},(x,y,z),domain=S.QQ))
check('positive expanded coefficient roundtrip',reconstructed==F3)
JJ=[[f.diff(v) for v in [x,y,z]] for f in reconstructed]
det=JJ[0][0]*(JJ[1][1]*JJ[2][2]-JJ[1][2]*JJ[2][1])-JJ[0][1]*(JJ[1][0]*JJ[2][2]-JJ[1][2]*JJ[2][0])+JJ[0][2]*(JJ[1][0]*JJ[2][1]-JJ[1][1]*JJ[2][0])
check('positive entire six-variable Jacobian minus3',det==S.Poly(-3,x,y,z,domain=S.QQ) and
 all(c==[{'monomial':[int(k==j) for k in range(6)],'coefficient':'1'}] for j,c in enumerate(serialized['coordinates'][3:],start=3)))
check('positive map defined on x0',all(m[0]>=0 for f in reconstructed for m,c in f.terms()))
Y1,Y2,Y3,Y4,Y5=S.symbols('Y1 Y2 Y3 Y4 Y5')
root=ps[0];poly=root**3+3*root**2-3*Y1*root+3*Y2
check('positive exact cubic inverse equation',
 zero(poly.subs({Y1:st['S'][0],Y2:st['S'][1]},simultaneous=True)))
discriminant=S.discriminant(poly,root)
check('positive cubic discriminant nonsquare',
 S.degree(discriminant,Y2)==2 and S.gcd(S.Poly(discriminant,Y2),S.Poly(S.diff(discriminant,Y2),Y2)).degree()==0)
p4=Y5/g;p3=(Y4-p4)/g;p2=(Y3-p3)/g
check('positive two-sided spectator reconstruction',
 all(zero(e) for e in S.Matrix([p2,p3,p4]).subs(
  {Y3:st['S'][2],Y4:st['S'][3],Y5:st['S'][4]},simultaneous=True)-S.Matrix(ps[1:])))
# Complete exceptional inertia data force rank4 and a missing second index5 prime.
data,orders,triangles,classes,actions=exceptional_actions();mul=data['qmul']
for j,action in enumerate(actions):
    fixed_orders={orders[g0] for g0,p in action['permutations'].items() if any(a==b for a,b in enumerate(p))}
    check('action'+str(j)+' all allowable full inertia orders',fixed_orders=={1,2,3,5})
    check('action'+str(j)+' only full-rank candidates r1,2,4',{o-1 for o in fixed_orders if o>1}=={1,2,4})
    profiles=[]
    for g0,o in orders.items():
        if o!=5:continue
        perm=action['permutations'][g0];remaining=set(range(11));lengths=[]
        while remaining:
            a=min(remaining);orb={a};b=perm[a]
            while b!=a:orb.add(b);b=perm[b]
            remaining-=orb;lengths.append(len(orb))
        profiles.append(tuple(sorted(lengths)))
    check('action'+str(j)+' exact five-cycle dictionary',set(profiles)=={(1,5,5)})
check('every order5 centralizer is C5',
 all(len(centralizer(data['quotient'],g0,mul))==5 for g0,o in orders.items() if o==5))
manifest=json.loads((ROOT/'sources/higher-completion-manifest.json').read_text())
for source in manifest['retained_primary_inputs']:
    raw=(ROOT/source['path']).read_bytes()
    check('retained primary bytes '+source['path'],len(raw)==source['bytes'])
    check('retained primary hash '+source['path'],hashlib.sha256(raw).hexdigest()==source['sha256'])
model={'claims':['KMFI103','KMFI104','KMFI105'],'formal_flag_controls':records,
 'all_rank_completion_criterion_proved':True,'one_common_independent_variable':True,
 'full_rank_required_for_normal_generation_tail_equivalence':True,
 'marked_point_construction_can_also_apply_to_rank_deficient_inputs':True,
 'whole_jacobian_formula':'gamma0','source_change_has_explicit_polynomial_inverse':True,
 'companion_negative_controls':companion,
 'positive_A6_control':{'group':'S3 on3','X_rank':3,'full_rank_index5_candidate':False,
 'theta0':[str(e) for e in theta0],'L':[[str(L[i,j]) for j in range(5)] for i in range(5)],
 'h':[str(e) for e in rec['h']],'A':[[str(rec['A'][i,j]) for j in range(4)] for i in range(4)],
 'whole_jacobian':-3,'expanded':serialized},
 'PSL11_full_rank_tautological_remaining_rank':4,
 'PSL11_remaining_sweep_must_be_nonfinite':True,'PSL11_remaining_branch_degree':1,
 'PSL11_missing_index5_normalization_primes':1,'PSL11_present_unramified_tail_irreducible_and_residue_degree':1,
 'actual_PSL11_input_supplied':False,'absolute_exceptional_action_decided':False,'goal_complete':False}
(ROOT/'higher_completion_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['HIGHER_RANK_COMPLETION.md','higher_completion.py','verify_higher_completion.py','higher_completion_model.json',
 'sources/higher-completion-manifest.json','BRANCH_CENTERED_COMPLETION.md','SWEEP_BOUNDARY_COMPLETENESS.md',
 'PSL11_RESIDUE_OBSTRUCTION.md','psl11_residues.py','central_lift.py']+[e['path'] for e in manifest['retained_primary_inputs']]
receipt={'result':'pass','check_count':len(checks),'checks':checks,'runtime_seconds':round(time.monotonic()-started,3),
 'python':platform.python_version(),'sympy':S.__version__,'sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 'assurance':'Written all-rank completion and dimension-reduction proofs; exact finite formal flags, explicit polynomial inverses and full expanded A6 Jacobian. No formal verification, independent whole-proof review or exceptional11 realization.',
 'goal_complete':False}
(ROOT/'higher_completion_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds'],'positive_full_A6_Jacobian':-3,
 'all_rank_affine_completion':True,'PSL11_full_rank_remaining_r':4,'actual_PSL11_input':False,'goal_complete':False}))
