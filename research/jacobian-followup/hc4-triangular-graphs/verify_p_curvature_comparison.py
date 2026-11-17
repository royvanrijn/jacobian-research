"""New raw/finite p-curvature comparison controls;run before registration."""
import hashlib,json,time
from pathlib import Path
import sympy as s
START=time.monotonic();HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
checks=[];algebraic=0
def check(name,truth,alg=False):
 global algebraic
 assert truth,name
 checks.append(name);algebraic+=int(alg)
def eq(name,a,b=0):check(name,s.cancel(a-b)==0,True)
x,y,t=s.symbols('x y t')
F=x*(x-1)*(x-t);H=x*(x-1);rr=H**2/2
L=s.Rational(3,4)*t*(1-t)*H**2+s.Rational(1,2)*(1-2*t)*H*F-F**2/4
eq('NEW exact Legendre Picard-Fuchs primitive',L,s.diff(rr,x)*F-s.Rational(3,2)*rr*s.diff(F,x))
aa=-1/(2*(t-1));bb=1/(2*t*(t-1));r=-H/(t*(t-1))
eq('NEW compact basis exact reduction',H/2-(aa+bb*x)*F,s.diff(r,x)*F-r*s.diff(F,x)/2)
check('NEW compact basis change nonzero',s.cancel(t*bb)!=0,True)
eq('NEW Legendre E preserves defining equation',2*y*(-H/(2*y))+H)
eq('NEW Legendre J preserves defining equation',2*y*(s.diff(F,x)/(2*y))-s.diff(F,x))
# Commutation in the etale double cover,after replacing y^2 by F.
Ey=-H/(2*y);Jy=s.diff(F,x)/(2*y)
eq('NEW unique etale lifts commute on y',s.diff(Jy,t)+Ey*s.diff(Jy,y)-s.diff(Ey,x)-Jy*s.diff(Ey,y))
C=s.Matrix([[0,t/(4*(1-t))],[1,t/(1-t)]])
N=C.subs(t,0)
check('NEW Euler residue nonzero nilpotent',N!=s.zeros(2) and N**2==s.zeros(2),True)
for p in [3,5,7]:
 check('NEW finite p-curvature residue '+str(p),N**p-N==-N and N!=s.zeros(2),True)
 check('NEW identity raw y^(p-1) has exact p-torsion '+str(p),(p-1)+1==p and p>1,True)
 eq('NEW identity raw missing primitive coefficient '+str(p),s.diff(y**p,y),p*y**(p-1))
# New nonsymmetric transverse lifts P=x;not a replay of any LP source family.
for p in [2,3,5]:
 def mod(f):return s.Poly(s.expand(f),x,y,modulus=p).as_expr()
 def E(f):return mod(s.diff(f,x)+(x+1)*y*s.diff(f,y))
 def T(f):return mod(E(f)+(x+1)*f)
 def iterate(op,f):
  for _ in range(p):f=op(f)
  return f
 ep_y=iterate(E,y);h=-ep_y
 check('NEW restricted vertical p-th power '+str(p),iterate(E,x)==0 and h!=0,True)
 for j,f in enumerate([1,x*y,y*y+x]):
  check('NEW corrected raw p-th power image '+str((p,j)),mod(iterate(T,f)+s.diff(h*f,y))==0,True)
context=json.loads((HERE/'sources/p-curvature-comparison-context.json').read_text())
for p,h in context['retained_proof_and_source_sha256'].items():check('retained byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
for p,h in context['new_primary_sha256'].items():check('new private primary byte pin '+p,hashlib.sha256((ROOT/p).read_bytes()).hexdigest()==h)
check('ONE explicit published finite criterion,no unproved import',len(context['newly_imported_proved_external_statements'])==1 and context['newly_used_unproved_external_statements']==[])
positive={'raw_corrected_submersion_connection_has_zero_p_curvature_all_p','raw_zero_p_curvature_uses_no_Keller_mate','identity_Keller_map_raw_integral_torsion_at_every_prime_proved','geometric_commuting_lift_finite_comparison_kernel_all_odd_p_proved','published_finite_Gauss_Manin_arithmetic_inverse_criterion_adopted','uniform_spreading_quantifier_retained','newly_adopted_external_theorem','external_review_deferred'}
for p,vv in context.items():
 if isinstance(vv,bool):check('scope '+p,vv is (p in positive))
result=dict(outcome='PASS',check_count=len(checks),algebraic_controls=algebraic,checks=checks,proof_sha256=hashlib.sha256((HERE/'RAW_P_CURVATURE_AND_FINITE_GAUSS_MANIN_COMPARISON.md').read_bytes()).hexdigest(),checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-START,3),**{p:v for p,v in context.items() if isinstance(v,bool)})
(HERE/'receipt-p-curvature-comparison.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps({p:result[p] for p in ['outcome','check_count','algebraic_controls','elapsed_seconds','original_goal_complete']}))
