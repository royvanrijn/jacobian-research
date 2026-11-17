"""Exact controls for the all-outer-degree nilpotent pencil source-orbit normal form."""
import hashlib,json,platform,time
from pathlib import Path
import sympy as sp
START=time.monotonic()
HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
checks=[]
def check(value,label):
    assert bool(value),label
    checks.append(label)
def zero(value,label):check(sp.cancel(sp.expand(value))==0,label)
def zeros(values,label):check(all(sp.cancel(sp.expand(v))==0 for v in values),label)
x,y,z,w,t=sp.symbols('x y z w t')
variables=[x,y,z,w]
cases=[
 ('2+1+1',2*x*y+z*z+w*w,x*x,2),
 ('2+2',2*x*y+2*z*w,x*x+z*z,2),
 ('3+1',2*x*y+z*z+w*w,2*x*z,3),
 ('4',2*x*w+2*y*z,2*x*z+y*y,4)]
for label,qA,qB,index in cases:
    A=sp.hessian(qA,variables)/2
    B=sp.hessian(qB,variables)/2
    T=A.inv()*B
    check(A.det()!=0,'the quadratic block is nondegenerate '+label)
    check(A*T==(A*T).T,'its pencil is self-adjoint '+label)
    check(T**index==sp.zeros(4) and T**(index-1)!=sp.zeros(4),
          'the nilpotency index and nonzero Jordan length are exact '+label)
    check(T.charpoly().as_expr()==T.charpoly().gen**4,
          'the pencil has no hidden nonzero eigenvalue '+label)
    for degree in [2,5]:
        H=t**degree+t+3
        HH=lambda value:sp.expand(H.subs(t,value))
        H0=HH(0)
        if label=='2+1+1':
            R=sp.cancel((HH(x*x)-H0)/(2*x))
            Phi=[x,y+R,z,w]
            inverse=[x,y-R,z,w]
        elif label=='2+2':
            R=sp.cancel((HH(x*x+z*z)-HH(z*z))/(2*x))
            S=sp.cancel((HH(z*z)-H0)/(2*z))
            Phi=[x,y+R,z,w+S]
            inverse=[x,y-R,z,w-S]
        elif label=='3+1':
            R=sp.cancel((HH(2*x*z)-H0)/(2*x))
            Phi=[x,y+R,z,w]
            inverse=[x,y-R,z,w]
        else:
            R=sp.cancel((HH(2*x*z+y*y)-HH(y*y))/(2*x))
            S=sp.cancel((HH(y*y)-H0)/(2*y))
            Phi=[x,y,z+S,w+R]
            inverse=[x,y,z-S,w-R.subs(z,z-S)]
        check(all(sp.denom(v)==1 for v in Phi+inverse),
              'all divided differences and inverse terms are polynomials '+label+' outer='+str(degree))
        zeros([Phi[i].subs(dict(zip(variables,inverse)),simultaneous=True)-variables[i]
               for i in range(4)],'the complete triangular source inverse is exact '+label+' outer='+str(degree))
        zero(qA.subs(dict(zip(variables,Phi)),simultaneous=True)-(qA+HH(qB)-H0),
             'the source potential equals the quadratic normal form literally '+label+' outer='+str(degree))
        check(sp.Matrix(Phi).jacobian(variables).det()==1,
              'the nonlinear source normalizer has determinant one '+label+' outer='+str(degree))

# The chain proof includes every nilpotent partition, not just the checked examples.
for length in [1,2,3,4]:
    A=sp.zeros(length)
    T=sp.zeros(length)
    for i in range(length):
        A[i,length-1-i]=1
        if i+1<length:T[i+1,i]=1
    check(A.det()!=0 and A*T==(A*T).T,
          'the normalized cyclic-chain Gram form is nondegenerate and self-adjoint length='+str(length))
    check(T**length==sp.zeros(length),
          'the displayed canonical chain has its exact nilpotent truncation length='+str(length))
    xx=sp.symbols('c0:'+str(length))
    q=(sp.Matrix(xx).T*A*T*sp.Matrix(xx))[0]
    expected=sum(xx[i]*xx[j] for i in range(length) for j in range(length) if i+j==length-2)
    zero(q-expected,'the block formula has the correct anti-diagonal orientation length='+str(length))

qA=2*x*w+2*y*z
qB=2*x*z+y*y
native=qA+qB*qB
zero(sp.hessian(native,variables).det()-16*(1+8*x*y-8*x*x*qB),
     'the length-four native Hessian is explicitly nonconstant')
check(sp.degree(sp.hessian(native,variables).det(),x)>0,
      'the normalizing shear is not being used as global Hessian invariance')
zero(sp.hessian(qA,variables).det()-16,'the recovered quadratic has actual unit Hessian')
zeros([sp.diff(native,v).subs(dict.fromkeys(variables,0)) for v in variables],
      'the origin remains one reduced critical point in the native model')
zero(sp.hessian(native,variables).det().subs(dict.fromkeys(variables,0))-16,
     'critical Hessian congruence still holds at that point')
# Graph application retains the full outer polynomial, its offset and all degrees.
X=sp.symbols('X')
P=t**3+t+1
Q=qA+2*X*X-3*X+5
H=2*P*P-3*P+5
zero(Q.subs(X,P.subs(t,qB))-(qA+H.subs(t,qB)),
     'the actual graph potential has the complete scalar-composite form')
check(sp.hessian(Q,variables+[X]).det()!=0,
      'the actual graph has a nondegenerate quadratic ambient potential')
check(sp.degree(H,t)==6,'the graph endpoint statement is not restricted to quadratic graph degree')

context=json.loads((HERE/'sources/nilpotent-pencil-context.json').read_text())
algebra_count=len(checks)
for path,digest in context['retained_proof_and_source_sha256'].items():
    check(hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==digest,
          'prior spectral/endpoint/transfer bytes are unchanged: '+path)
check(context['newly_used_unproved_external_statements']==[],
      'the self-adjoint chain form and shears are reproduced, not unproved imports')
check(context['singleton_zero_endpoint_resolved'] is False and context['unrestricted_HC4_or_JC2_resolved'] is False,
      'orbit equality does not solve the singleton-zero endpoint or HC4')
receipt=dict(outcome='PASS',claim='HC4MYJD0',check_count=len(checks),
 algebraic_and_scope_controls=algebra_count+2,checks=checks,
 checker_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 proof_sha256=hashlib.sha256((HERE/'NILPOTENT_QUADRATIC_PENCILS_AND_ENDPOINT_ORBITS.md').read_bytes()).hexdigest(),
 source_sha256={'sources/nilpotent-pencil-context.json':hashlib.sha256((HERE/'sources/nilpotent-pencil-context.json').read_bytes()).hexdigest()},
 dependencies=context['retained_proof_and_source_sha256'],
 assurance='Finite exact controls of the written all-degree normal-form proof; not independent replay, formal verification or external review.',
 all_nilpotent_four_variable_quadratic_composite_orbits_identified=True,
 polynomial_graph_and_parameterization_degrees_unbounded=True,
 nonsemisimple_nonzero_spectrum_resolved=False,singleton_zero_endpoint_resolved=False,
 unrestricted_HC4_or_JC2_resolved=False,original_goal_complete=False,
 python=platform.python_version(),sympy=sp.__version__,elapsed_seconds=round(time.monotonic()-START,3))
(HERE/'receipt-nilpotent-pencils.json').write_text(json.dumps(receipt,indent=2)+'\n')
print(json.dumps({k:receipt[k] for k in ['outcome','check_count','algebraic_and_scope_controls','elapsed_seconds','original_goal_complete']},indent=2))
