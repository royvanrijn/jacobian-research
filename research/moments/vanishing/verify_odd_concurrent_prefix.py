"""Auxiliary exact centered-odd concurrent and constant-W frame seeds.
Universal all-degree reductions, full class cover and retained analytic/Riccati
arguments remain written inputs. No previous checker main is executed.
"""
from pathlib import Path
from functools import lru_cache
from itertools import combinations
import hashlib,json,time
import sympy as s
P=Path(__file__).resolve().parent
t,u,v,T,y,Z,Wc=s.symbols('t u v T y Z Wc')
sha=lambda f:hashlib.sha256(f.read_bytes()).hexdigest()
basis=[(0,0),(1,0),(0,1),(1,1)]

def phase(m,h,e,f,a,b,c,ju,jv):
    J=ju*u+jv*v;q=(a*u*u+2*c*u*v+b*v*v)/2-J*J/(2*m)
    L=e*u+f*v-h*J/m
    return s.expand(u**3+v**3-q+t*t*L+h*h*t**4/(2*m))

def universal():
    m,h,e,f,a,b,c,ju,jv=s.symbols('m h e f a b c ju jv',nonzero=True)
    J=ju*u+jv*v;q2=(a*u*u+2*c*u*v+b*v*v)/2;G=u**3+v**3+e*u+f*v
    original=t*(h*T+G)-m*T*T/2-T*J-q2
    completed=s.expand(original.subs(T,T+(t*h-J)/m))
    assert s.expand(completed-(-m*T*T/2+t*G-q2+(t*h-J)**2/(2*m)))==0
    W=phase(m,h,e,f,a,b,c,ju,jv)
    assert s.expand(t*t*(t*G-q2+(t*h-J)**2/(2*m)).subs({u:u/t,v:v/t},simultaneous=True)-W)==0
    metric=s.Matrix([[a,c,ju],[c,b,jv],[ju,jv,m]])
    schur=s.Matrix([[a-ju*ju/m,c-ju*jv/m],[c-ju*jv/m,b-jv*jv/m]])
    assert s.cancel(metric.det()-m*schur.det())==0
    old=t*s.diff(W,t)-2*W-2*t*t
    S=s.expand(old+s.Rational(2,3)*(u*s.diff(W,u)+v*s.diff(W,v)+2*t*t))
    assert s.Poly(S,t,u,v).total_degree()<=4 and s.Poly(S,u,v).total_degree()<=2
    assert s.Poly(S.subs(h,0),t,u,v).total_degree()<=3
    gauge=[]
    for i in range(5):
        for k in range(5-i):
            z=u**i*v**k;cu=s.Rational(2,3)*u*z;cv=s.Rational(2,3)*v*z
            diff=(S-old)*z+s.Rational(2,3)*t*t*(u*s.diff(z,u)+v*s.diff(z,v))
            assert s.expand(diff-t*t*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
            gauge.append([i,k])
    R=W-u**3-v**3;principals=[]
    for i in range(11):
        for k in range(11-i):
            if (i,k) in basis:continue
            cu=u**(i-2)*v**k/3 if i>=2 else s.Integer(0)
            cv=u**i*v**(k-2)/3 if i<2 else s.Integer(0)
            rem=s.expand(-cu*s.diff(R,u)-cv*s.diff(R,v)-t*t*(s.diff(cu,u)+s.diff(cv,v)))
            assert s.expand(u**i*v**k-rem-t*t*(s.diff(cu,u)+s.diff(cv,v))-s.diff(W,u)*cu-s.diff(W,v)*cv)==0
            assert s.Poly(rem,t,u,v).total_degree()<=i+k and s.Poly(rem,u,v).total_degree()<i+k
            principals.append([i,k])
    return {'Gaussian_completion':True,'Schur_determinant':True,'centered_rescaling':True,'source_joint_bound':4,'h_zero_source_joint_bound':3,'gauge_identities':gauge,'principal_identities':principals}

def binary_model(m,h):
    W=phase(s.Integer(m),s.Integer(h),5,7,11,13,1,2,3);old=t*s.diff(W,t)-2*W-2*t*t
    S=s.expand(old+s.Rational(2,3)*(u*s.diff(W,u)+v*s.diff(W,v)+2*t*t));Rlow=W-u**3-v**3
    ring=s.QQ.poly_ring(t);poly=lambda z:s.Poly(z,u,v,domain=ring);zt=s.Poly(0,t,domain=s.QQ);zero=poly(0);wp,dp=poly(W),poly(t*t)
    @lru_cache(None)
    def reduce(ex,depth):
        if ex in basis:
            out=[zt]*4;out[basis.index(ex)]=s.Poly(1,t);return tuple(out),zero,zero
        assert depth>0
        i,k=ex;cu=u**(i-2)*v**k/3 if i>=2 else s.Integer(0);cv=u**i*v**(k-2)/3 if i<2 else s.Integer(0)
        rem=poly(-cu*s.diff(Rlow,u)-cv*s.diff(Rlow,v)-t*t*(s.diff(cu,u)+s.diff(cv,v)))
        out=[zt]*4;Cu,Cv=poly(cu),poly(cv)
        for f,c in rem.terms():
            if not c:continue
            vec,x,z=reduce(f,depth-1);ct=s.Poly(c,t,domain=s.QQ)
            out=[a+ct*b for a,b in zip(out,vec)];Cu+=poly(c)*x;Cv+=poly(c)*z
        return tuple(out),Cu,Cv
    R=s.Integer(1);rows=[];cost=3 if h==0 else 4
    for j in range(1,5):
        R=s.expand(t**3*s.diff(R,t)-3*(j-1)*t*t*R+s.Rational(2,3)*t*t*(u*s.diff(R,u)+v*s.diff(R,v))+S*R)
        assert s.Poly(R,t,u,v).total_degree()<=cost*j
        out=[zt]*4;Cu,Cv=zero,zero
        for ex,c in poly(R).terms():
            if not c:continue
            vec,x,z=reduce(ex,2*j-1);ct=s.Poly(c,t,domain=s.QQ)
            out=[a+ct*b for a,b in zip(out,vec)];Cu+=poly(c)*x;Cv+=poly(c)*z
        normal=sum((poly(z.as_expr()*u**ex[0]*v**ex[1]) for z,ex in zip(out,basis)),zero)
        assert (poly(R)-normal-dp*(Cu.diff(u)+Cv.diff(v))-wp.diff(u)*Cu-wp.diff(v)*Cv).is_zero
        limits=[cost*j-sum(ex) for ex in basis];assert all(z.degree()<=b for z,b in zip(out,limits))
        rows.append({'j':j,'depth':2*j-1,'N':[str(z.as_expr()) for z in out],'H':str(t**(3*j)),'H_valuation':3*j,'row_degree_bounds':limits,'witness':{'Cu':str(Cu.as_expr()),'Cv':str(Cv.as_expr())}})
    det=s.Matrix([[11,1,2],[1,13,3],[2,3,m]]).det();assert det!=0
    return {'m':m,'h':h,'metric_determinant':int(det),'W':str(W),'source':str(S),'cost':cost,'seed_columns':rows}

def frame_phase(a,g30,g21,g12,g03,r,z):
    return s.expand(g30*a**3*t**6+g21*a*a*t**4*y+g12*a*t*t*y*y+g03*y**3+r*a*t**4+z*t*t*y-y*y/2)
def wt(expr):return max((e[0]+2*e[1] for e,c in s.Poly(s.expand(expr),t,y).terms() if c),default=-100)

def frame_universal():
    a,b,c,d,g,r,z=s.symbols('a b c d g r z',nonzero=True);W=frame_phase(a,b,c,d,g,r,z)
    old=t*s.diff(W,t)-2*W-t*t;S=s.expand(old+s.Rational(2,3)*(y*s.diff(W,y)+t*t));assert wt(S)<=6 and s.degree(S,y)<=2
    gauge=[]
    for n in range(9):
        f=y**n;C=s.Rational(2,3)*y*f
        assert s.expand((S-old)*f+s.Rational(2,3)*t*t*y*s.diff(f,y)-t*t*s.diff(C,y)-s.diff(W,y)*C)==0;gauge.append(n)
    principal=[]
    for n in range(2,13):
        C=y**(n-2)/(3*g);rem=s.expand(-C*s.diff(W-g*y**3,y)-t*t*s.diff(C,y))
        assert wt(rem)<=2*n and s.degree(rem,y)<n
        assert s.expand(y**n-rem-t*t*s.diff(C,y)-s.diff(W,y)*C)==0;principal.append(n)
    q=1-2*d*a*t*t;E=b*a**3*t**4+r*a*t*t;F=(c*a*a*t**3+z*t)**2
    numerator=s.expand(-s.diff(q,t)*q/2+s.diff(E,t)*q*q+(s.diff(F,t)*q-F*s.diff(q,t))/2)
    assert s.degree(numerator,t)<=7 and numerator.subs(t,0)==0
    assert s.cancel(s.diff(-s.log(q)/2+E+F/(2*q),t)-numerator/q**2)==0
    return {'gauge_identities':gauge,'principal_identities':principal,'weighted_source_bound':6,'Gaussian_flat_log_derivative_numerator':str(numerator),'Gaussian_flat_numerator_degree':7,'Gaussian_flat_denominator':str(q*q)}

def frame_model(a,g):
    W=frame_phase(s.Integer(a),2,3,5,s.Integer(g),7,11);S=s.expand(t*s.diff(W,t)-2*W-t*t+s.Rational(2,3)*(y*s.diff(W,y)+t*t));wp=s.Poly(W,y,domain=s.QQ.poly_ring(t));R=s.Integer(1);rows=[]
    for j in range(1,4):
        R=s.expand(t**3*s.diff(R,t)-3*(j-1)*t*t*R+s.Rational(2,3)*t*t*y*s.diff(R,y)+S*R);normal=s.expand(R);witness=s.Integer(0)
        while s.degree(normal,y)>=2:
            n=s.degree(normal,y);c=s.Poly(normal,y).coeff_monomial(y**n)/(3*g);C=c*y**(n-2);witness+=C
            normal=s.expand(normal-t*t*s.diff(C,y)-s.diff(W,y)*C)
        assert s.expand(R-normal-t*t*s.diff(witness,y)-s.diff(W,y)*witness)==0
        vals=[s.Poly(normal,y).coeff_monomial(y**k) for k in range(2)]
        assert all(s.degree(z,t)<=6*j-2*k for k,z in enumerate(vals))
        rows.append({'j':j,'N':[str(z) for z in vals],'H':str(t**(3*j)),'H_valuation':3*j,'row_degree_bounds':[6*j,6*j-2],'witness':str(witness)})
    return {'a':a,'g03':g,'W':str(W),'source':str(S),'seed_columns':rows}

def gaussian_flat_model(a,b,c,d,r,z):
    poly=a*Wc+b*Z**3+c*Z*Z*T+d*Z*T*T+r*Z+z*T
    moments=[]
    for n in range(9):
        val=0
        for (i,j,k),coef in s.Poly(s.expand(poly**n),Wc,Z,T).terms():
            if i==j and k%2==0:val+=coef*s.factorial(i)*s.factorial2(k-1)
        moments.append(val)
    q=1-2*d*a*t*t;E=b*a**3*t**4+r*a*t*t;F=(c*a*a*t**3+z*t)**2
    N=s.expand(-s.diff(q,t)*q/2+s.diff(E,t)*q*q+(s.diff(F,t)*q-F*s.diff(q,t))/2)
    I=sum(m*t**n/s.factorial(n) for n,m in enumerate(moments));defect=s.Poly(s.expand(q*q*s.diff(I,t)-N*I),t)
    assert all(defect.coeff_monomial(t**k)==0 for k in range(8)) and moments[0]==1
    assert all(moments[k]==0 for k in [1,3,5,7])
    return {'parameters':[a,b,c,d,r,z],'M0_to_M8':[str(z) for z in moments],'q':str(q),'log_derivative_numerator':str(N),'numerator_degree':int(s.degree(N,t)) if N else -1}

def arithmetic():
    out={}
    for cost,expected in [(3,26),(4,36)]:
        gates=[sum(cost*j for j in range(1,size+1))-sum(sum(basis[k]) for k in chosen) for size in range(1,5) for chosen in combinations(range(4),size) if 0 in chosen]
        assert max(gates)==expected;out[str(cost)]=max(gates)
    assert max(6,6*3-2)==16 and max(36,16,8)==36
    cover={'quadratic':7,'smooth':94,'nodal':100,'cusp':94,'transverse':98,'tangent':94,'triangle':96,'concurrent':36,'repeated':98}
    assert max(cover.values())==100
    return out,cover

def main():
    start=time.monotonic();univ=universal();fu=frame_universal();bm=[binary_model(m,h) for m,h in [(1,0),(1,1),(2,0),(2,2)]]
    fm=[frame_model(a,g) for a,g in [(0,1),(1,1),(1,2),(2,1)]]
    gm=[gaussian_flat_model(*a) for a in [(0,2,3,5,7,11),(1,2,3,5,7,11),(1,2,3,0,7,11),(2,1,2,3,4,5)]]
    gates,cover=arithmetic();cert={'binary_universal':univ,'frame_universal':fu,'binary_models':bm,'cubic_frame_models':fm,'Gaussian_flat_models':gm,'binary_gates_by_source_cost':gates,'odd_cover':cover}
    target=P/'odd_concurrent_prefix_certificate.json';target.write_text(json.dumps(cert,indent=2,default=int)+'\n')
    report={'completed':True,'checker_sha256':sha(Path(__file__)),'certificate_sha256':sha(target),'runtime_seconds':round(time.monotonic()-start,3),'sympy_version':s.__version__,'binary_source_gauge_identities':len(univ['gauge_identities']),'binary_principal_identities':len(univ['principal_identities']),'cubic_frame_gauge_identities':len(fu['gauge_identities']),'cubic_frame_principal_identities':len(fu['principal_identities']),'binary_models':len(bm),'binary_seed_identities':sum(len(m['seed_columns']) for m in bm),'cubic_frame_models':len(fm),'cubic_frame_seed_identities':sum(len(m['seed_columns']) for m in fm),'Gaussian_flat_models':len(gm),'Gaussian_flat_direct_positive_moments':32,'odd_nonisotropic_concurrent_gate':36,'odd_nonisotropic_h_zero_gate':26,'odd_isotropic_cubic_height_gate':16,'odd_isotropic_flat_gate':8,'odd_concurrent_full_gate':36,'odd_sufficient_pure_prefix':100,'global_sufficient_pure_prefix_unchanged':130,'homogeneous_sufficient_pure_prefix_unchanged':96,'least_prefix_interval':[20,130],'odd_least_prefix_interval':[20,100],'homogeneous_least_prefix_interval':[14,96],'odd_cover':cover,'optimality_proved':False,'written_analytic_or_global_proof_certified_by_checker':False,'boundary':'Exact centered completion, rescaling/Schur identity, radial Euler source and principal divisions,28 actual seed identities in eight models, flat Gaussian logarithmic derivative and32 direct moment checks, cofactor/cover arithmetic. All-degree reduction and full vertex/class/frame cover, retained nonisotropic no-constant-period and formal Airy Riccati obstruction remain written inputs. No independent replay, external review, formal verification or expanded radical cofactors.'}
    (P/'odd_concurrent_prefix_verification.json').write_text(json.dumps(report,indent=2)+'\n');print(json.dumps(report,indent=2))
if __name__=='__main__':main()
