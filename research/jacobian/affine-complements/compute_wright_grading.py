"""Exact algebra for the Wright grading obstruction and trace countercontrol.

Universal weight and trace proofs are written in the associated notes;
finite graded controls do not establish all-height coverage by themselves.
"""
from pathlib import Path
import sympy as S,json,hashlib
HERE=Path(__file__).resolve().parent
X,Y,z,n,l,e=S.symbols('X Y z n l e')
f,g,df,dg=S.symbols('f g df dg')
formal=(n*f+2*z*df)*(l*g+z*dg)-z*df*(e*g+2*z*dg)
expected=n*l*f*g+n*z*f*dg+(2*l-e)*z*df*g
assert S.expand(formal-expected)==0
# The actual finite quartic countercontrol has an invariant smooth ramification
# line, a homogeneous volume and a regular semisimple hyperbolic grading.
x,w,p,q,t=S.symbols('x w p q t')
P=x+w+x**3;Q=-x*x/2-x*w-S.Rational(3,4)*x**4
J=S.expand(S.diff(P,x)*S.diff(Q,w)-S.diff(P,w)*S.diff(Q,x))
assert J==w
E=lambda F:S.expand((-x*S.diff(F,x)+2*w*S.diff(F,w))/3)
F=x**4+2*x*x-4*p*x-4*q
# Multiplication by x on the exact rank-four target-module basis.
M=S.zeros(4)
for j in range(4):
    remainder=S.rem(x**(j+1),F,x)
    for i in range(4):M[i,j]=S.expand(remainder).coeff(x,i)
traces={j:S.expand(S.trace(M**j)) for j in range(5)}
def traced(expr):
    expr=S.Poly(S.rem(S.expand(expr.subs(w,p-x-x**3)),F,x),x)
    return S.expand(sum(coef*traces[power[0]] for power,coef in expr.terms()))
TP=traced(E(P));TQ=traced(E(Q))
assert S.expand(TP+S.Rational(52,3)*p)==0
assert S.expand(TQ-S.Rational(64,3)*q-8)==0
assert S.diff(TP,p)+S.diff(TQ,q)==4
pc=t+t**3;qc=-t*t/2-S.Rational(3,4)*t**4
transverse=S.factor(TP.subs({p:pc,q:qc})*S.diff(qc,t)-TQ.subs({p:pc,q:qc})*S.diff(pc,t))
disc=S.factor(S.discriminant(F,x))
disc_derivative=S.expand(TP*S.diff(disc,p)+TQ*S.diff(disc,q))
assert disc.subs({p:0,q:0})==0
assert disc_derivative.subs({p:0,q:0})!=0
inverse=p*(1-12*q)/(1+9*p*p+4*q)
assert S.cancel(inverse.subs({p:pc,q:qc})-t)==0
# Record the pole-filtration degree of the real Wright grading on V1.
A=X*Y;B=X*X*Y;D=X**3*Y+X;C=1+B+D*D
old=[-Y-C*C-2*A*D-B*D*D,-A-D-2*B*D-D**3,C,D]
delta=lambda H:S.expand(-X*S.diff(H,X)+2*Y*S.diff(H,Y))
new=[2*old[0]+2*C*C-4*old[1]*D,old[1]-2*D+4*C*D,-2*D*D,-D]
for f0,g0 in zip(old,new):assert S.expand(delta(f0)-g0)==0
payload=dict(as_of='2026-11-14',outcome='EXACT_GRADING_IDENTITIES_COMPUTED',
 grading_weights=dict(X=-1,Y=2),
 weight_modules=dict(zero='C[B]',negative='A_-n=D^n*C[B] for n>=1',
       even_positive='A_2m=Y^m*C[B] for m>=0',odd_positive='A_(2m+1)=A*Y^m*C[B] for m>=0'),
 homogeneous_jacobian_formula='B^(l-1)*(n*l*f*g+n*B*f*gprime+(n+1)*B*fprime*g)',
 pole_generator_derivatives=['2*f0+2*f2**2-4*f1*f3','f1-2*f3+4*f2*f3','-2*f3**2','-f3'],
 trace_countercontrol=dict(source='A2_(x,w), not the Wright normalization or a Keller counterexample',
   P=str(P),Q=str(Q),jacobian=str(J),normalized_grading='(-x*d_dx+2*w*d_dw)/3',
   monic_quartic=str(F),multiplication_matrix_x=[[str(M[i,j]) for j in range(4)] for i in range(4)],
   power_traces={str(j):str(value) for j,value in traces.items()},
   target_trace_P=str(TP),target_trace_Q=str(TQ),target_divergence='4',
   branch_parameterization=[str(pc),str(qc)],branch_rational_inverse=str(inverse),
   branch_transversality=str(transverse),branch_transversality_at_zero=str(transverse.subs(t,0)),
   discriminant=str(disc),trace_derivative_of_discriminant_at_zero=str(disc_derivative.subs({p:0,q:0}))),
 source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 assurance=dict(unrestricted_wright_pairs_excluded=False,trace_tangency_to_actual_wright_branch_proved=False,
   any_keller_counterexample_constructed=False,homogeneous_obstruction='written all-height proof'))
(HERE/'wright_grading_certificates.json').write_text(json.dumps(payload,indent=2)+'\n')
print('Formal homogeneous identity verified; Wright grading raises B-poles by at most one')
print('Quartic control trace:',TP,TQ,'divergence4; transverse at branch origin:',transverse.subs(t,0))
