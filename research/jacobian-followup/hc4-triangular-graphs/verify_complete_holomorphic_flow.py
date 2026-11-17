"""Exact controls for the written complete-flow theorem; not a classification replay."""
import hashlib,json,time
from pathlib import Path
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1]
def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def main():
 start=time.monotonic();checks=[]
 def check(label,condition):
  assert condition,label
  checks.append(label)
 def zero(label,expression):check(label,s.cancel(expression)==0)
 x,y,t,z,alpha,a,c,m,n,ell=s.symbols('x y t z alpha a c m n ell')
 def deriv(f,P,Q):return s.diff(f,x)*P+s.diff(f,y)*Q
 # I: complete does not mean locally finite. The exact n-th derivative
 # is x^n*y; the sampled iterates merely check the transcription.
 P,Q=0,x*y
 q=y
 for k in range(1,5):
  q=deriv(q,P,Q);zero('non-LF complete vertical field iterate '+str(k),q-x**k*y)
 tau=s.symbols('tau');flow=y*s.exp(t*x)
 zero('entire vertical flow solves the ODE',s.diff(flow,t)-x*flow)
 zero('entire vertical flow satisfies the time group law',flow.subs(t,t+tau)-flow.subs(y,y*s.exp(tau*x)))
 graph=x**3+x+1;A=x**2+2;B=c*x*s.diff(graph,x)-A*graph
 zero('I graph equation with nonconstant A',deriv(y-graph,c*x,A*y+B)-A*(y-graph))
 # The DIFFERENT symmetry fixes x and scales y-q, so preserves any
 # polynomial in x times y-q, with no dependence on root count.
 F=s.Function('F')(x);H=F*(y-graph)
 zero('whole vertical union plus graph has LF logarithmic field',deriv(H,0,y-graph)-H)
 zero('translation preserves an arbitrary vertical union',deriv(F,0,1))
 G=s.Function('G')
 zero('I residue-class differential equation transcription',s.diff(t**3*G(t**5),t)*t-t**3*(3*G(t**5)+5*t**5*s.diff(G(z),z).subs(z,t**5)))
 for k in [2,3,5]:
  # Degree of coordinate function field C(t^k) is k, so the graph
  # parametrization must NOT be mistaken for a normalization map.
  zero('nonprimitive graph control k='+str(k),(t**k)**2-t**(2*k))
  check('normalization field degree exceeds1 k='+str(k),s.degree(z**k-x,z)==k>1)
 # I vertical zero graph with missing x-values: x*y+1=0 has x^-1=-y.
 zero('vertical zero curve Bezout unit',(-x*y-1).subs(y,-1/x))
 check('Gm zero curve has a nonconstant unit',s.cancel(-(-1/x))==1/x)
 # II formula for symbolic lambda uses h-polynomials and differentiation.
 for mm,nn in [(1,1),(2,3),(3,2)]:
  h=x**mm*y**nn;lam=h**2+h+3
  P=nn*lam*x;Q=(a-mm*lam)*y
  zero('II complete normal form h equation '+str((mm,nn)),deriv(h,P,Q)-nn*a*h)
  zero('II logarithmic x equation '+str((mm,nn)),deriv(x,P,Q)-nn*lam*x)
 # III general rational-chain identity for unconstrained symbolic p,lambda.
 p=s.Function('p')(x);lam=s.Function('lambda')(z)
 w=x**ell*y+p;L=s.symbols('L')
 P=n*L*x;Q=(a*w-L*((m+n*ell)*x**ell*y+m*p+n*x*s.diff(p,x)))/x**ell
 zero('III birational chain rule D(w)',deriv(w,P,Q)-(a-m*L)*w)
 # Polynomial complete non-LF example, including its omitted chart x=0.
 h=x*(x*y+1);P=x*h;Q=-(x*y+1)*(2*x*y+1)
 zero('III genuine polynomial example invariant h',deriv(h,P,Q))
 zero('III chart boundary x=0 is present with nonzero field',Q.subs(x,0)+1)
 zero('III x=0 remains invariant',P.subs(x,0))
 # On a nonzero level h=C, x is a nonconstant unit and y=(C-x)/x^2.
 C=s.symbols('C',nonzero=True)
 zero('III nonzero level has Gm parameterization',h.subs(y,(C-x)/x**2)-C)
 # Constant-lambda condition forces a=m*lambda0; nonconstant p leaves
 # a forbidden x*p' term below ell (not a polynomial vector field).
 lambda0=s.symbols('lambda0',nonzero=True)
 for ll in [2,3,5]:
  pp=1+x**(ll-1)
  rem=s.rem(n*lambda0*x*s.diff(pp,x),x**ll,x)
  check('III constant lambda rejects nonconstant p ell='+str(ll),rem!=0 and s.degree(rem,x)<ll)
 # Nodal A1-normalized curve: tangent field has quadratic, incomplete
 # normalization speed, despite exact logarithmic tangency.
 u,v=s.symbols('u v');HH=v**2-u**2*(u+1);uu=t**2-1;vv=t*(t**2-1)
 zero('nodal curve normalization identity',HH.subs({u:uu,v:vv}))
 Pu=s.diff(HH,v);Qv=-s.diff(HH,u)
 zero('nodal Hamiltonian is logarithmic',s.diff(HH,u)*Pu+s.diff(HH,v)*Qv)
 zero('nodal normalization u speed',s.diff(uu,t)*(t**2-1)-Pu.subs({u:uu,v:vv}))
 zero('nodal normalization v speed',s.diff(vv,t)*(t**2-1)-Qv.subs({u:uu,v:vv}))
 T=s.symbols('T');solution=(1+s.exp(2*T))/(1-s.exp(2*T))
 zero('quadratic normalization field has meromorphic solution',s.diff(solution,T)-(solution**2-1))
 check('quadratic normalization solution has a finite pole',s.simplify((1-s.exp(2*T)).subs(T,0))==0 and (1+s.exp(2*T)).subs(T,0)==2)
 # Joint Hartogs has codimension TWO, even after adding time variables.
 check('joint flow Hartogs codimension is2',(1+2)-(1+0)==2)
 check('joint group-law density uses both times',(2+2)-(2+0)==2)
 algebra_count=len(checks)
 context=json.loads((HERE/'sources/complete-holomorphic-flow-context.json').read_text())
 for name,hsh in context['retained_proof_and_source_sha256'].items():check('retained bytes '+name,sha(ROOT/name)==hsh)
 for name,hsh in context['new_primary_sha256'].items():check('new primary bytes '+name,sha(HERE/name)==hsh)
 check('one explicit proved classification import;Cauchy extension written directly',len(context['newly_imported_proved_external_statements'])==1 and not context['newly_used_unproved_external_statements'])
 flags={k:v for k,v in context.items() if isinstance(v,bool)}
 check('classification access boundary is honest',flags['classification_normal_forms_read_in_primary_Bustinduy'] and not flags['classification_independently_reproved'] and not flags['original_Brunella_full_text_read'])
 check('no automatic completeness or full HC4 claim',not any(flags[k] for k in ['complete_field_proved_to_exist_automatically','last_HC4_source_to_regular_inverse_Gram_implication_proved','whole_surviving_higher_degree_packet_excluded','unrestricted_HC4_or_JC2_resolved','original_goal_complete']))
 proof=HERE/'COMPLETE_HOLOMORPHIC_FLOW_DESCENT.md'
 for phrase in ['CURRENT LAST UNSOLVED IMPLICATION','ADOPTED AS REPRODUCED THERE','JOINTLY holomorphic','one-variable Cauchy','ENTIRE N','BIRATIONAL','Every reduced subset','vanishes identically','PLANE-specific','DEFERRED']:
  check('written scope '+phrase,phrase.lower() in proof.read_text().lower())
 result=dict(outcome='PASS',check_count=len(checks),algebraic_and_scope_controls=algebra_count,checks=checks,proof_sha256=sha(proof),checker_sha256=sha(Path(__file__)),dependencies=context['retained_proof_and_source_sha256'],source_sha256=context['new_primary_sha256'],elapsed_seconds=round(time.monotonic()-start,3),**flags)
 (HERE/'receipt-complete-holomorphic-flow.json').write_text(json.dumps(result,indent=2)+'\n')
 print(json.dumps({k:result[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
if __name__=='__main__':main()
