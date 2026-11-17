"""Exact KMFI115--117 reconstruction, rank and moving-flag controls."""
import os
os.environ.setdefault('SYMPY_GROUND_TYPES','flint')
from pathlib import Path
import sympy as S
from sympy.polys.matrices import DomainMatrix
import time,json,hashlib,platform,math
from integrated_sweeps import family
from polynomial_wave import polynomial_chart
from sweep_differential import *
ROOT=Path(__file__).resolve().parent;started=time.monotonic();checks=[]
def check(name,truth):
    if not truth:raise AssertionError(name)
    checks.append(name)
def zero(e):return S.cancel(e)==0
def matrix_zero(M):return all(zero(e) for e in M)
t=S.Symbol('t');controls=[]
for r in [2,3,4]:
    f=family(r,t*t,t)
    controls.append(('integrated_rank'+str(r),f['parameters'],f['X']))
w,z=S.symbols('w z');binary=polynomial_chart(w*w,[(1,0)],w,z)
controls.append(('binary_rank2',[w,z],binary['X']))
models=[]
for label,ps,X in controls:
    r=len(ps);oldf=X[0];g=X[1:,:];P=(g.jacobian(ps)-S.Matrix(ps)*gradient(oldf,ps)).applyfunc(S.expand)
    d=from_g(g,ps,oldf.subs({p:0 for p in ps}))
    check(label+' Euler reconstruction all coordinates',matrix_zero(d['X']-X))
    check(label+' nilpotency all trace coefficients',all(zero(e) for e in nilpotent_residuals(P)))
    intr=from_matrix(P,ps)
    check(label+' all first-order curl conditions',all(zero(e) for e in intr['curl']))
    check(label+' reconstructed gradient exact',matrix_zero(gradient(intr['f'],ps)-gradient(oldf,ps)))
    check(label+' closedness of intrinsic differential',all(zero(e) for e in intr['closed']))
    check(label+' exact entire matrix reconstruction',matrix_zero(intr['rebuilt']-P))
    check(label+' reconstructed vector derivatives exact',matrix_zero(intr['g'].jacobian(ps)-g.jacobian(ps)))
    C=coframe(oldf,P,ps);detC=S.factor(C.det(method='domain-ge'))
    check(label+' full cyclic coframe',detC!=0 and X.jacobian(ps).rank()==r)
    power=S.eye(r)
    for m in range(1,r+1):
        power=(power*P).applyfunc(S.expand)
        check(label+' entire divergence power'+str(m),
          matrix_zero(divergence_columns(power,ps)+(r-m)*gradient(oldf,ps)*(P**(m-1))))
    # Exact polynomial identity at the generic point, with a frozen kernel direction.
    top=(P**(r-1)).applyfunc(S.expand)
    col=next(j for j in range(r) if any(e!=0 for e in top[:,j]));u=top[:,col]
    pivot=next(i for i in range(r) if u[i]!=0);v=(top[pivot,:]/u[pivot]).applyfunc(S.cancel)
    check(label+' exact rank-one last power',matrix_zero(top-u*v))
    check(label+' kernel and left kernel',matrix_zero(P*u) and matrix_zero(v*P))
    lam=(gradient(oldf,ps)*u)[0];D=directional_matrix(P,ps,u)
    check(label+' nonzero kernel derivative eigenvalue',lam!=0 and matrix_zero(v*D-lam*v))
    check(label+' nilpotent last-power trace pairing',zero((v*u)[0]))
    models.append({'label':label,'r':r,'coframe_determinant':str(detC),
                   'kernel_direction_eigenvalue':str(S.factor(lam)),
                   'all_Piola_and_eigenvector_residuals':'zero'})
# Reject a trace-only arbitrary coefficient input.
ps=list(S.symbols('p1:5'));g=S.Matrix([p*ps[0] for p in ps]);radial=from_g(g,ps)
check('Euler parameterization forces trace but not higher nilpotency',zero(S.trace(radial['P'])) and
      any(not zero(e) for e in nilpotent_residuals(radial['P'])[1:]))
try:from_g(S.Matrix(ps),ps)
except ValueError:checks.append('nonzero linear trace rejected')
else:raise AssertionError('nonzero linear trace accepted')
# Matrix reconstruction must not discard its differential compatibility test.
for r in [2,4]:
    ps=list(S.symbols('p1:'+str(r+1)));bad=S.zeros(r);bad[0,1]=ps[0]**2
    intr=from_matrix(bad,ps)
    check('rank'+str(r)+' nilpotent incompatible matrix rejected',
      all(zero(e) for e in nilpotent_residuals(bad)) and not matrix_zero(intr['rebuilt']-bad))
    if r==2:
        check('rank2 curl alone is tautological countercontrol',
           all(zero(e) for e in intr['curl']) and any(not zero(e) for e in intr['closed']))
    else:check('rank4 first-order curl detects incompatibility',any(not zero(e) for e in intr['curl']))
# The new affine-matrix obstruction is sharp, even with regular Jordan type.
sharp=[]
for r in range(2,7):
    ps=list(S.symbols('p1:'+str(r+1)));g=S.Matrix([ps[0]**2/2,*ps[2:],0]);d=from_g(g,ps);P=d['P'];X=d['X']
    check('rank'+str(r)+' affine sharp control entire X',X==S.Matrix([ps[0],ps[0]**2/2,*ps[2:],0]))
    check('rank'+str(r)+' affine sharp control nilpotent',all(zero(e) for e in nilpotent_residuals(P)))
    check('rank'+str(r)+' affine sharp control regular generically',P.rank()==r-1 and not matrix_zero(P**(r-1)))
    check('rank'+str(r)+' affine sharp control rank collapse',X.jacobian(ps).rank()==r-1 and
            zero(coframe(d['f'],P,ps).det(method='domain-ge')))
    gam=S.Symbol('gamma');swept=X+gam*S.Matrix([1,*ps])
    check('rank'+str(r)+' affine sharp control whole normalized Jacobian',
          zero(swept.jacobian([gam,*ps]).det(method='domain-ge')-gam**r))
    sharp.append({'r':r,'rank_DX':r-1,'generic_Jordan_block_size':r,'degree_P':1})
# The left-eigenvalue identity at the actual marked branch point.
f=family(4,t*t,t);ps=f['parameters'];X=f['X'];P=(X[1:,:].jacobian(ps)-S.Matrix(ps)*gradient(X[0],ps)).applyfunc(S.expand)
pt={p:(1 if i==0 else 0) for i,p in enumerate(ps)};P0=P.subs(pt);u=P0**3*S.eye(4)[:,3]
D=directional_matrix(P,ps,u).subs(pt);T=S.Symbol('T')
check('rank4 actual marked kernel vector',u==S.Matrix([1,1,S.Rational(4,3),S.Rational(2,15)]))
check('rank4 actual directional characteristic polynomial',
 S.expand(D.charpoly(T).as_expr()-(T-2)*(T-1)**2*(T+4))==0)
check('rank4 actual derivative nonnilpotent',(gradient(X[0],ps).subs(pt)*u)[0]==2 and S.trace(D**2)!=0)
H=haantjes_component(P,ps,0,0,2)
check('Haantjes replacement fails valid full-rank input',H!=0)

# Exact dimensions of the linear compatibility space, before nilpotency.
import itertools
def linear_dimensions(r,n):
    ps=list(S.symbols('x1:'+str(r+1)))
    mons=[S.prod(p**j for p,j in zip(ps,exps)) for exps in itertools.product(range(n+1),repeat=r) if sum(exps)==n]
    rowids={};entries={};closedrows=set();column=0
    for i in range(r):
     for j in range(r):
      for mon in mons:
       P=S.zeros(r);P[i,j]=mon;d=compatibility(P,ps)
       equations=[('trace',0,d['trace'])]+[('curl',k,e) for k,e in enumerate(d['curl'])]+[('closed',k,e) for k,e in enumerate(d['closed'])]
       for kind,k,e in equations:
        for powers,co in S.Poly(e,*ps).terms():
         if co==0:continue
         key=(kind,k,powers)
         if key not in rowids:rowids[key]=len(rowids)
         row=rowids[key];entries[row,column]=co
         if kind=='closed':closedrows.add(row)
       column+=1
    mat=S.MutableSparseMatrix(len(rowids),column,entries)
    rank=DomainMatrix.from_Matrix(mat).convert_to(S.QQ).rank()
    openrows=[i for i in range(len(rowids)) if i not in closedrows]
    no_closed=mat.extract(openrows,list(range(column)))
    rank_no=DomainMatrix.from_Matrix(no_closed).convert_to(S.QQ).rank()
    return dict(r=r,n=n,unknown_coefficients=column,rank_constraints=rank,
                dimension=column-rank,dimension_without_closedness=column-rank_no)
dimensions=[]
for r,n in [(2,0),(2,1),(2,2),(3,1),(3,2),(4,0),(4,1),(4,2)]:
    d=linear_dimensions(r,n);expected=r*r-1 if n==0 else r*math.comb(n+r,r-1)
    check('rank'+str(r)+' homogeneous degree'+str(n)+' exact compatibility-space dimension',d['dimension']==expected)
    if r>=3:check('rank'+str(r)+' degree'+str(n)+' first-order system already forces closedness',
                   d['dimension_without_closedness']==d['dimension'])
    dimensions.append(d)
model={'claims':['KMFI115','KMFI116','KMFI117'],
 'all_degree_coefficient_reconstruction':'f_n=div(g_(n+1))/n; impose traceless linear g and traces P^2,...,P^r',
 'rank4_remaining_nonlinear_equations':['tr(P^2)=0','tr(P^3)=0','tr(P^4)=0'],
 'intrinsic_matrix_reconstruction':'first-order curl plus trace; r2 also needs da=0; r>=3 closedness follows',
 'homogeneous_linear_space_dimensions':dimensions,'full_rank_controls':models,
 'sharp_affine_controls':sharp,
 'generic_full_rank_test':'det(df;dfP;...;dfP^(r-1))!=0',
 'kernel_direction_derivative_identity':'v DP[u]=df(u) v, with df(u)!=0 at every full-rank point',
 'actual_rank4_point_directional_derivative':[[str(e) for e in row] for row in D.tolist()],
 'actual_rank4_directional_characteristic_polynomial':str(S.factor(D.charpoly(T).as_expr())),
 'actual_Haantjes_nonzero_component':str(S.factor(H)),
 'new_obstruction':'full-rank inputs cannot be affine P, lie in a fixed nilpotent matrix space, or have all g components degree<=2',
 'fixed_point_Jordan_frame_excluded':False,
 'arbitrary_quadratic_Keller_maps_excluded':False,
 'all_nonlinear_nilpotent_inputs_classified':False,
 'actual_PSL11_input_supplied':False,'absolute_exceptional11_action_decided':False,
 'full_action_spectrum_classified':False,'goal_complete':False}
(ROOT/'sweep_differential_model.json').write_text(json.dumps(model,indent=2)+'\n')
paths=['SWEEP_DIFFERENTIAL_RECONSTRUCTION.md','sweep_differential.py','verify_sweep_differential.py',
 'sweep_differential_model.json','integrated_sweeps.py','polynomial_wave.py',
 'HIGHER_RANK_COMPLETION.md','ALL_RANK_INTEGRATED_ACTIONS.md']
receipt={'result':'pass','check_count':len(checks),'checks':checks,
 'runtime_seconds':round(time.monotonic()-started,3),'python':platform.python_version(),'sympy':S.__version__,
 'sha256':{f:hashlib.sha256((ROOT/f).read_bytes()).hexdigest() for f in paths},
 'assurance':'Unbounded written Euler/curl reconstruction, cyclic coframe, Piola and derivative-eigenvalue proof, with a universal affine/nilpotent-space rank-collapse obstruction. Exact controls supplement the proofs. No formal, independent whole-proof or external review claimed.',
 'goal_complete':False}
(ROOT/'sweep_differential_receipt.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({'result':'pass','checks':len(checks),'runtime_seconds':receipt['runtime_seconds'],
 'linear_space_dimensions':[(d['r'],d['n'],d['dimension']) for d in dimensions],
 'actual_rank4_directional_eigenvalues':[2,1,1,-4],
 'actual_PSL11_field_supplied':False,'goal_complete':False}))
