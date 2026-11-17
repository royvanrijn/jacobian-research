#!/usr/bin/env python3
"""Narrow exact feedback, partial-factor and product-action controls; not universal proof replay."""
from pathlib import Path
from itertools import product,permutations
import ast,hashlib,json,platform,time
import sympy as s
HERE=Path(__file__).resolve().parent;ROOT=HERE.parents[1];checks=[]
def ck(label,value):
    assert value,label
    checks.append(label)
def sub(F,variables,values):return [f.subs(dict(zip(variables,values)),simultaneous=True) for f in F]
def zero(f):return s.expand(f)==0
def mul(g,h):return tuple(g[h[i]] for i in range(len(g)))
def group(gens,n):
    identity=tuple(range(n));out={identity};todo=[identity]
    while todo:
        g=todo.pop()
        for h in gens:
            a=mul(h,g)
            if a not in out:out.add(a);todo.append(a)
    return out

def primitive(G,n):
    if {g[0] for g in G}!=set(range(n)):return False
    for j in range(1,n):
        parent=list(range(n))
        def find(i):
            while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
            return i
        for g in G:
            a,b=find(g[0]),find(g[j])
            if a!=b:parent[a]=b
        if len({find(i) for i in range(n)})>1:return False
    return True

def main():
    started=time.monotonic();x,y,z,a,pi=s.symbols('x y z a pi');variables=[x,y,z]
    source=ROOT/'jacobian-followup/universal-keller-compression/verify.py';nodes=[n for n in ast.parse(source.read_text()).body if isinstance(n,ast.FunctionDef) and n.name=='source_map'];ck('reuse pure compiler only',len(nodes)==1)
    env={'x':x,'y':y,'z':z};exec(compile(ast.Module(body=nodes,type_ignores=[]),str(source),'exec'),env);compiler=env['source_map']
    K=s.Matrix(compiler(3,{}));DK=K.jacobian(variables);mu=s.Matrix([s.Rational(1,3),-1,s.Rational(2,3)]);delta=DK.adjugate()*mu
    ck('auxiliary cubic actual Jacobian one',zero(DK.det()-1))
    ck('lifted translation vector polynomial',[s.denom(s.cancel(v)) for v in delta]==[1,1,1])
    ck('lifted translation has actual target direction',all(zero(v) for v in DK*delta-mu))
    R1=K[1]+3*K[0];R2=K[2]-2*K[0]
    ck('first exact polynomial integral',zero(s.Matrix([R1]).jacobian(variables).dot(delta)))
    ck('second exact polynomial integral',zero(s.Matrix([R2]).jacobian(variables).dot(delta)))
    ck('nonzero source-coordinate response at origin',delta.subs({x:0,y:0,z:0})==s.Matrix([s.Rational(1,3),2,s.Rational(1,3)]))
    S=x/(1+x*y);relation=S+R1*S*S-R2/2+K[0]*(S**3-3*S*S-1)
    ck('rational kernel parametrization exact cleared identity',s.cancel(relation)==0)
    profile_controls=[]
    sigma=K[0];w=s.symbols('w')
    for profiles in [{4:a},{5:a},{4:a*a+a,5:a*a},{4:a*pi+a*a*pi*pi,5:a+a*a*pi}]:
        n=max(profiles);coeff={j:profiles.get(j,s.Integer(0)).subs(pi,sigma) for j in range(4,n+1)};F=s.Matrix(compiler(n,coeff));A=F.jacobian(variables).subs(y,0);V=F.diff(a).subs(y,0)
        ck('profile response determinant '+str(profiles),zero(A.det()-1))
        response=s.expand((A.adjugate()*V)[0]);B=sum(s.diff(f,a).subs(pi,z)*(x*z)**j for j,f in profiles.items());C=sum(j*(j-2)*f.subs(pi,z)*(x*z)**j for j,f in profiles.items());expected=sum((j+1)*s.diff(f,a).subs(pi,z)*(x*z)**j for j,f in profiles.items())-3*x*(x*z)*B-C*B/x
        ck('profile exact all-degree response formula '+str(profiles),s.cancel(response-expected)==0)
        laurent=s.Poly(s.expand(response.subs(z,w/x)*x**(2*n+5)),x);ck('profile nonzero loop response '+str(profiles),not laurent.is_zero)
        profile_controls.append({'profiles':{str(j):str(f) for j,f in profiles.items()},'slice_response_nonzero':True})
    # Actual cyclic dependence on both first outputs, with two polynomial processing orders.
    X=s.symbols('X:3');Y=s.symbols('Y:3');coords=list(X)+list(Y);sig1=sub(K,variables,X)[0];sig2=sub(K,variables,Y)[0]
    a1=sig1*sig1+sig2;b2=sig1;c2=sig1*sig2
    F1=sub(compiler(4,{4:a1}),variables,X)
    F2=sub(compiler(5,{4:b2,5:c2}),variables,Y)
    # The parameter expressions above use X,Y, while the compiler's x,y,z are replaced independently.
    F=F1+F2;P0=F1+list(Y);Q0=list(X)+sub(compiler(5,{4:X[0],5:X[0]*sig2}),variables,Y)
    P1=list(X)+F2;Q1=sub(compiler(4,{4:sig1*sig1+Y[0]}),variables,X)+list(Y)
    ck('coupled blocks first processing order exact',all(zero(v-q) for v,q in zip(sub(Q0,coords,P0),F)))
    ck('coupled blocks reverse processing order exact',all(zero(v-q) for v,q in zip(sub(Q1,coords,P1),F)))
    ck('coupled profiles have genuine cycle',a1.has(*Y) and c2.has(*X))
    # Whole nine-variable naive feedback: actual Jacobian fails, not merely a formal surrogate.
    blocks=[s.symbols('p'+str(i)+'_0:3') for i in range(2)];U=s.symbols('u:3');allcoords=[v for P in blocks for v in P]+list(U);L=lambda v:(3*v*v-v)/2
    Ctrace=sum(L(P[0]) for P in blocks);selector=lambda T:-T*s.prod(T-j for j in range(2,5))/6;beta=1-selector(Ctrace);KU=s.Matrix(sub(K,variables,U))
    bad=[f for P in blocks for f in sub(compiler(4,{4:U[0]}),variables,P)]+list(KU-beta*mu)
    point=dict(zip(allcoords,[1,0,1,0,0,0,0,0,0]));J=s.Matrix(bad).jacobian(allcoords).subs(point);bad_det=J.det()
    ck('actual nine-variable naive feedback Jacobian exactly43over18',bad_det==s.Rational(43,18))
    ck('actual Schur product explains Jacobian defect',bad_det==1+s.Rational(25,6)*s.Rational(1,3))
    # Every recovered control is unaffected by subtracting beta*mu.
    Vout=KU-beta*mu;auxR1=KU[1]+3*KU[0];auxR2=KU[2]-2*KU[0];h=auxR1+auxR2*auxR2;recovered=(Vout[1]+3*Vout[0])+(Vout[2]-2*Vout[0])**2
    ck('nonlinear recovered control exact reverse-factor identity',zero(recovered-h))
    ck('nonlinear recovered control annihilated by delta',zero(s.Matrix([R1+R2*R2]).jacobian(variables).dot(delta)))
    # Exact finite product-action groups distinguish arithmetic and geometric primitivity.
    group_controls=[]
    for N in [2,3]:
        points=list(product(range(3),repeat=N));index={P:i for i,P in enumerate(points)};gens=[]
        for coordinate in range(N):
            for perm in [(1,0,2),(1,2,0)]:
                gens.append(tuple(index[tuple(perm[v] if j==coordinate else v for j,v in enumerate(P))] for P in points))
        G0=group(gens,len(points));rotation=tuple(index[P[1:]+P[:1]] for P in points);G=group(gens+[rotation],len(points))
        ck(f'rank{N} geometric product group exact order',len(G0)==6**N)
        ck(f'rank{N} geometric product imprimitive',not primitive(G0,len(points)))
        ck(f'rank{N} arithmetic product-action group exact order',len(G)==N*6**N)
        ck(f'rank{N} transitive constant action makes arithmetic group primitive',primitive(G,len(points)))
        if N==2:
            def blocks_containing_origin(H):
                out=[]
                for mask in range(1<<(len(points)-1)):
                    B={0}|{j+1 for j in range(len(points)-1) if mask>>j&1}
                    if all(not (C:={g[i] for i in B})&B or C==B for g in H):out.append(sorted(B))
                return out
            geo_blocks=blocks_containing_origin(G0);arith_blocks=blocks_containing_origin(G)
            ck('all geometric intermediate blocks exact Boolean diamond',sorted(map(len,geo_blocks))==[1,3,3,9])
            ck('all arithmetic intermediate blocks primitive interval',sorted(map(len,arith_blocks))==[1,9])
        group_controls.append({'rank':N,'degree':3**N,'geometric_order':len(G0),'arithmetic_order':len(G),'geometric_primitive':False,'arithmetic_primitive':True,'top_group':'cyclic transitive order'+str(N)})
    # Actual quadratic restriction of scalars conjugates to two cubic maps over splitting constants.
    D,theta=s.symbols('D theta');Qvars=s.symbols('q:6');Qpairs=[Qvars[2*i]+theta*Qvars[2*i+1] for i in range(3)]
    field=s.QQ.frac_field(D);mod=s.Poly(theta*theta-D,theta,domain=field)
    def reduce_theta(expr):return s.rem(s.Poly(s.expand(expr),theta,domain=s.QQ.frac_field(D,*Qvars)),s.Poly(theta*theta-D,theta,domain=s.QQ.frac_field(D,*Qvars))).as_expr()
    outputs=[reduce_theta(f) for f in sub(K,variables,Qpairs)]
    coeffs=[]
    for f in outputs:coeffs.extend([s.expand(f).coeff(theta,0),s.expand(f).coeff(theta,1)])
    plus=[coeffs[2*i]+theta*coeffs[2*i+1] for i in range(3)];minus=[coeffs[2*i]-theta*coeffs[2*i+1] for i in range(3)];conjinputs=[q.subs(theta,-theta) for q in Qpairs]
    ck('quadratic Weil map actual plus split conjugacy',all(reduce_theta(f-g)==0 for f,g in zip(plus,sub(K,variables,Qpairs))))
    ck('quadratic Weil map actual minus split conjugacy',all(reduce_theta(f-g)==0 for f,g in zip(minus,sub(K,variables,conjinputs))))
    ck('quadratic Weil map has six polynomial coordinates',len(coeffs)==6 and all(s.Poly(f,*Qvars,domain=s.QQ.poly_ring(D)).as_expr()==f for f in coeffs))
    manifest=json.loads((HERE/'sources.json').read_text())
    for src in manifest['local_sources']:ck('local pin '+src['path'],hashlib.sha256((ROOT/src['path']).read_bytes()).hexdigest()==src['sha256'])
    for src in manifest['external_sources']:
        data=(HERE/src['path']).read_bytes();ck('primary pin '+src['path'],len(data)==src['bytes'] and hashlib.sha256(data).hexdigest()==src['sha256'])
    receipt={'schema_version':1,'result':'pass','check_count':len(checks),'checks':checks,'profile_response_controls':profile_controls,'actual_bad_scalar_feedback_jacobian':str(bad_det),'product_action_controls':group_controls,'quadratic_Weil_control':'Six polynomial coordinates; exact conjugacy to two cubic maps over theta^2=D. Generic arithmetic group theorem is the written Galois argument, not inferred from the permutation sample.','python':platform.python_version(),'sympy':s.__version__,'runtime_seconds':round(time.monotonic()-started,3),'sha256':{f:hashlib.sha256((HERE/f).read_bytes()).hexdigest() for f in ['PROOF.md','LITERATURE_AND_GAPS.md','verify.py','sources.json']},'assurance':'Exact cubic translation/inverse identities, unbounded-theorem response controls, literal coupled partial factors, an actual failed nine-variable scalar feedback, nonlinear output-control recovery, finite geometric/arithmetic product actions and actual quadratic restriction of scalars. No formalization of all-profile Jacobians, independent splitting-field regularity, arithmetic constant fields, Boolean reconstruction rigidity, polynomial-kernel descent or the complete feedback classification. No independent whole-proof replay, external review or priority certificate.'}
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2,ensure_ascii=False)+'\n');print(json.dumps({'result':'pass','checks':len(checks),'seconds':receipt['runtime_seconds']}))
if __name__=='__main__':main()
