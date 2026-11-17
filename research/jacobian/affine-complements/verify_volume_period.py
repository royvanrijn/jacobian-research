"""Exact controls for the global period obstruction and marked Wright ring.

The affine cover, global exactness argument and conditional compatibility
comparison are written proofs, not formal verification by these controls.
"""
from pathlib import Path
from time import monotonic
import sympy as S,json,hashlib,platform
HERE=Path(__file__).resolve().parent
start=monotonic();checks=[]
def check(name,value):
    assert value,name
    checks.append(name)
def zero(name,value):check(name,S.cancel(value)==0)
x,y,u,v,a,b,s,t,alpha,beta,lam=S.symbols('x y u v a b s t alpha beta lam')
local={str(z):z for z in (x,y,u,v,a,b,s,t,alpha,beta,lam)}
def parse(value):return S.sympify(value,locals=local)
producer=HERE/'compute_volume_period.py'
source_hash=hashlib.sha256(producer.read_bytes()).hexdigest()
period=json.loads((HERE/'volume_period_certificates.json').read_text())
check('period producer seal',period['source_sha256']==source_hash)
check('two full boundary charts retained',len(period['cases'])==2)
for row in period['cases']:
    name=row['name'];marked=name=='contact_at_marked';r=1 if marked else 2
    tr={x:1/u,y:u**3*v+u**r}
    jac=S.det(S.Matrix([[S.diff(tr[x],u),S.diff(tr[x],v)],
                       [S.diff(tr[y],u),S.diff(tr[y],v)]]))
    zero('chart Jacobian '+name,jac+u)
    cox_t={s:x,t:1,alpha:y,beta:1-x*y} if marked else {s:x,t:1,alpha:1-x*x*y,beta:y}
    cox_s={s:1,t:u,alpha:1+u*u*v,beta:-v} if marked else {s:1,t:u,alpha:-v,beta:1+u*v}
    for symbol,weight in ((s,1),(t,1),(alpha,-1),(beta,-2)):
        zero('Cox chart transition '+name+' '+str(symbol),
             u**weight*S.sympify(cox_t[symbol]).subs(tr,simultaneous=True)-cox_s[symbol])
    B=alpha*s+beta*t*t if marked else alpha*t+beta*s*s
    E=(a*s+b*t)*B+(beta*s**3 if marked else -alpha*s*s)
    zero('B trivial on t chart '+name,B.subs(cox_t,simultaneous=True)-1)
    zero('B trivial on s chart '+name,B.subs(cox_s,simultaneous=True)-1)
    f0=parse(row['local_volume_f0']);f1=parse(row['local_volume_f1'])
    zero('volume is the E equation on t chart '+name,f0-E.subs(cox_t,simultaneous=True))
    zero('volume is minus E on s chart '+name,f1+E.subs(cox_s,simultaneous=True))
    zero('full two-form gluing '+name,jac*f0.subs(tr,simultaneous=True)-f1)
    F0=parse(row['local_antiderivative_F0']);F1=parse(row['local_primitive_on_u_chart'])
    zero('first local primitive differential '+name,S.diff(F0,y)-f0)
    zero('second local primitive differential '+name,-S.diff(F1,v)-f1)
    difference=S.expand(-F0.subs(tr,simultaneous=True)*S.diff(tr[x],u)-F1)
    zero('literal primitive overlap difference '+name,difference-parse(row['primitive_difference']))
    residue=S.expand(u*difference).coeff(u,0)
    zero('exact global period '+name,residue-(b if marked else a))
    zero('saved period '+name,residue-parse(row['period']))
    zero('primitive difference is horizontal '+name,S.diff(difference,v))

laurent=sum(S.Symbol('h'+str(j+8))*u**j for j in range(-8,9))
zero('Laurent derivatives have zero residue',S.expand(u*S.diff(laurent,u)).coeff(u,0))
check('marked period never zero in allowed parameter domain',period['cases'][0]['period']=='b')
check('away zero-period condition is a=0',period['cases'][1]['period']=='a')
eta0=parse(period['remaining_family']['global_primitive_x_coefficient'])
eta1=parse(period['remaining_family']['global_primitive_u_coefficient'])
zero('remaining primitive differential on x chart',-S.diff(eta0,y)-(b-x*x+x**4*y))
zero('remaining primitive differential on u chart',-S.diff(eta1,v)-(-b*u-v))
zero('remaining primitive glues regularly',
     -eta0.subs({x:1/u,y:u**3*v+u*u},simultaneous=True)/u**2-eta1)
B=alpha*t+beta*s*s;E=(a*s+b*t)*B-alpha*s*s
psi={s:lam*s,t:t,alpha:alpha,beta:beta/lam**2}
zero('normalizing automorphism preserves B',B.subs(psi,simultaneous=True)-B)
zero('full parameter normalization',E.subs(psi,simultaneous=True)
     -lam**2*E.subs({a:a/lam,b:b/lam**2},simultaneous=True))
zero('a=0,b=lambda^2 is the fixed boundary pair',
     E.subs({a:0,b:lam**2}).subs(psi,simultaneous=True)
     -lam**2*((t*t-s*s)*alpha+s*s*t*beta))
zero('remaining fixed determinant',t*(s*s*t)-s*s*(t*t-s*s)-s**4)

# Independent literal verification of the two subring inclusions.
X,Y,A,B,D,F0,F1,F2,F3=S.symbols('X Y A B D F0 F1 F2 F3')
local.update({str(z):z for z in (X,Y,A,B,D,F0,F1,F2,F3)})
marking=json.loads((HERE/'wright_marking_certificates.json').read_text())
check('marked-ring producer seal',marking['source_sha256']==source_hash)
al=1+X*X*Y;be=-al*al-Y;tt=1+X*X*al
old=[be*tt*tt,be*X*tt,al*tt,al*X]
standard=[Y,X*Y,X*X*Y,X**3*Y+X]
sub={A:standard[1],B:standard[2],D:standard[3]}
for j,(original,saved) in enumerate(zip(old,marking['old_generators_in_wright'])):
    zero('old generator in Wright ring '+str(j),original-parse(saved).subs(sub,simultaneous=True))
for j,(expected,saved) in enumerate(zip(standard,marking['wright_generators_in_old'])):
    zero('Wright generator in old ring '+str(j),expected-parse(saved).subs(
         dict(zip((F0,F1,F2,F3),old)),simultaneous=True))
relations=list(map(parse,marking['kernel_generators']))
for j,rel in enumerate(relations):zero('source relation '+str(j),rel.subs(sub,simultaneous=True))
gb=S.groebner(relations,D,B,A,Y,order='lex')
check('stored exact kernel basis',list(map(str,gb))==marking['kernel_groebner_basis'])
elim=S.groebner([A-X*Y,B-X*X*Y,D-X**3*Y-X],X,D,B,A,Y,order='lex')
kernel=[e.as_expr() for e in elim.polys if not e.as_expr().has(X)]
check('full elimination ideal, not just vanishing equations',
      list(S.groebner(kernel,D,B,A,Y,order='lex'))==list(gb))
def modzero(name,value):
    num=S.cancel(value).as_numer_denom()[0]
    zero(name,gb.reduce(S.expand(num))[1])
for chart,Xvalue in [('Y_nonzero',A/Y),('Bplus1_nonzero',D/(B+1))]:
    for symbol,expected in [(A,Xvalue*Y),(B,Xvalue**2*Y),(D,Xvalue**3*Y+Xvalue)]:
        modzero('source-open inverse '+chart+' '+str(symbol),symbol-expected)
near={Y:u**3*v-u*u,A:u*u*v-u,B:u*v-1,D:v}
for j,rel in enumerate(relations):zero('near-boundary relation '+str(j),rel.subs(near,simultaneous=True))
for symbol in (Y,A,B,D):
    modzero('near-boundary inverse '+str(symbol),
            symbol-near[symbol].subs({u:A/B,v:D},simultaneous=True))
zero('uniformizer inverse near R',(A/B).subs(near,simultaneous=True)-u)
check('R has the full affine coordinate v',near[D]==v and
      [near[z].subs(u,0) for z in (Y,A,B)]==[0,0,-1])
for symbol,expr,expected in [(X,1/u,-1),(Y,near[Y],2),(A,near[A],1),(B,near[B],0),(D,near[D],0)]:
    num,den=S.fraction(S.cancel(expr))
    def val(p):return min(m[0] for m,c in S.Poly(p,u).terms() if c!=0)
    check('actual omitted valuation '+str(symbol),val(num)-val(den)==expected)
    check('stored omitted valuation '+str(symbol),marking['omitted_valuations'][str(symbol)]==expected)
zero('Wright volume has simple zero along R',
     S.det(S.Matrix([[S.diff(1/u,u),0],[S.diff(near[Y],u),S.diff(near[Y],v)]]))+u)
# Check the period-normalized form in the actual global source coordinates.
r=X/tt;z=tt*tt*be
volume=1-r*r+r**4*z
zero('normalized E equation on source slice',volume-1/tt)
source_jac=S.diff(r,X)*S.diff(z,Y)-S.diff(r,Y)*S.diff(z,X)
zero('normalized period form is minus source volume',volume*source_jac+1)

report=dict(as_of='2026-11-14',outcome='PASS',check_count=len(checks),checks=checks,
    all_height_periods=dict(contact_at_marked='b != 0: excluded',contact_away_marked='a: forces a=0'),
    surviving_complex_boundary_orbits=1,
    surviving_boundary=dict(B=['t','s**2'],E=['t**2-s**2','s**2*t']),
    marked_ring='C[Y,XY,X^2Y,X^3Y+X] inside C[X,Y]',
    omitted_prime=['Y','A','B+1'],omitted_source_values=dict(X=-1,Y=2),
    subsequent_all_height_pole_three_exclusion="PF2RULED8 (separate intrinsic proof and receipt)",
            frontier_superseded_by='PF2RULED16',
    frontier_handoff_scope='historical pre-PF2RULED16; this receipt does not prove the current frontier',
    first_open_reduced_pairs=[[6,8],[7,8]],
    frontier_adoption='PF2RULED12; subsequent handoff, outside this checker intrinsic proof coverage',
    assurance=dict(exactness_obstruction='written all-height proof',
       marked_comparison='literal equality inside the actual polynomial source ring',
       actual_ruled_compatibility_proved=False,remaining_zero_period_model_excluded=False,
       finite_map_constructed=False,quartic_packets_excluded=False,
       formal_verification=False,independent_full_reproof=False,external_review=False),
    source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [
       Path(__file__),HERE/'compute_volume_period.py',HERE/'volume_period_certificates.json',
       HERE/'wright_marking_certificates.json',HERE/'VOLUME_PERIOD_OBSTRUCTION.md',HERE/'WRIGHT_MARKED_MODEL.md']},
    python=platform.python_version(),sympy=S.__version__,elapsed_seconds=round(monotonic()-start,6))
(HERE/'volume_period_receipt.json').write_text(json.dumps(report,indent=2)+'\n')
print(f'PASS {len(checks)} exact controls; one zero-period boundary orbit; marked Wright ring identified; finite maps still open')
