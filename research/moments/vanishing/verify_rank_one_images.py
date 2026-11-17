"""Finite exact controls; not the general effective-field/module algorithm."""
from pathlib import Path
from itertools import product
import hashlib,json
import sympy as s
p=Path(__file__).resolve().parent
t,z,e,Z,W,T=s.symbols('t z e Z W T')
def eq(a,b):assert s.cancel(a-b)==0,(a,b)
def compose(A,B):
    out={}
    for i,ai in A.items():
        for j,bj in B.items():
            for k in range(i+1):
                d=i-k+j
                out[d]=out.get(d,0)+ai*s.binomial(i,k)*s.diff(bj,t,k)
    return {j:s.cancel(v) for j,v in out.items() if v!=0}
def clear(L):
    q=s.lcm([s.denom(s.cancel(v)) for v in L.values()])
    return {j:s.expand(s.cancel(q*v)) for j,v in L.items() if v!=0}
def apply(L,f):return s.cancel(sum(v*s.diff(f,t,j) for j,v in L.items()))
def integer_roots(poly,negative=False):
    ans=[]
    for a in s.roots(s.Poly(s.expand(poly),z),z):
        if a.is_Integer and ((a<0) if negative else (a>=0)):ans.append(int(a))
    return ans
def indicial(L,at):
    terms={}
    for j,v in L.items():
        pp=s.Poly(s.expand(v.subs(t,t+at)),t)
        if pp.is_zero:continue
        order=min(m[0] for m,cc in pp.terms() if cc!=0)
        terms[j]=(order-j,pp.nth(order))
    mu=min(v[0] for v in terms.values())
    return s.expand(sum(v[1]*s.ff(z,j) for j,v in terms.items() if v[0]==mu))
def infinity(L):
    terms={j:(s.degree(v,t)-j,s.Poly(v,t).LC()) for j,v in L.items() if v!=0}
    h=max(v[0] for v in terms.values())
    C=s.expand(sum(v[1]*s.ff(z,j) for j,v in terms.items() if v[0]==h))
    return max(integer_roots(C),default=-1)
def jet(f,n):return s.Integer(0) if n<0 else s.expand(s.series(f,t,0,n+1).removeO())
def polynomial_test(L,f):
    B=infinity(L);candidate=jet(f,B);S=apply(L,candidate)
    M=L if S==0 else clear(compose({1:S,0:-s.diff(S,t)},L))
    N=max(integer_roots(indicial(M,0)),default=-1)
    equal=jet(f-candidate,N)==0
    return {'polynomial':bool(equal),'candidate':candidate if equal else None,
            'infinity_degree_bound':B,'origin_test_through':N}
def rational_test(L,f):
    L=clear(L);r=max(L);D=s.Integer(1);poles=[]
    if r==0:
        eq(f,0);return {'rational':True,'value':s.Integer(0),'bound_denominator':D,'candidate_poles':[]}
    for a in s.roots(L[r],t):
        if a==0:continue
        B=indicial(L,a)
        k=max([-n for n in integer_roots(B,negative=True)],default=0)
        if k:D*= (t-a)**k;poles.append((str(a),k))
    D=s.simplify(s.expand(D));assert D.subs(t,0)!=0
    LD=clear(compose(L,{0:1/D}))
    decision=polynomial_test(LD,D*f)
    value=s.cancel(decision['candidate']/D) if decision['polynomial'] else None
    return {'rational':decision['polynomial'],'value':value,'bound_denominator':D,
            'candidate_poles':poles,'endpoint_test':decision}

cases=[
 ('zero',{2:s.Integer(1)},s.Integer(0),True),
 ('polynomial',{1:t,0:s.Integer(-3)},t**3,True),
 ('delayed_rational',{1:t*(1-t),0:-5+4*t},t**5/(1-t),True),
 ('logarithm',{2:1-t,1:s.Integer(-1)},s.log(1-t),False),
 ('exponential',{1:s.Integer(1),0:s.Integer(-1)},s.exp(t),False),
 ('algebraic_nonrational',{1:2*(1-t),0:s.Integer(1)},s.sqrt(1-t),False)]
rr=(1-t)**-2+(1-2*t)**-1
num,den=s.fraction(s.cancel(rr))
cases.append(('two_poles',{1:num*den,0:-(s.diff(num,t)*den-num*s.diff(den,t))},rr,True))
cases.append(('delayed_logarithm',clear(compose({2:1-t,1:s.Integer(-1)},{0:t**-5})),t**5*s.log(1-t),False))
receipts=[]
for name,L,f,expected in cases:
    eq(apply(L,f),0)
    rec=rational_test(L,f);assert rec['rational']==expected,(name,rec)
    if expected:eq(rec['value'],f)
    receipts.append({'name':name,**rec})
assert receipts[2]['bound_denominator']!=1
assert receipts[-1]['endpoint_test']['origin_test_through']>=6

# Isolate simple poles from one rational generator, then generate all higher poles.
f=3/(t-1)**3+5/(t-2)**2+t**2+1
isolation=[]
for a,k,other,l in [(1,3,2,2),(2,2,1,3)]:
    h=s.cancel((t-a)**(k-1)*(t-other)**l*f)
    nn,dd=s.fraction(h);poly,rem=s.div(nn,dd,t)
    residue=s.limit((t-a)*h,t,a);assert residue!=0
    eq((h-poly)/residue,1/(t-a))
    for j in range(5):eq(s.diff(1/(t-a),t,j),(-1)**j*s.factorial(j)/(t-a)**(j+1))
    isolation.append({'pole':a,'initial_order':k,'isolated_residue':str(residue)})
# A polynomial generator is needed: exp(t) is rank-one but fails rationality.
assert not receipts[4]['rational']

# Exact Bernstein identity for squarefree divisors with different pole counts.
bchecks=[]
for d in [1-t,(t-1)*(t-2),(t**2+1)*(t-2)]:
    A,B,g=s.gcdex(d,s.diff(d,t),t);eq(g,1)
    eq(A*d+B*s.diff(d,t),1)
    eq(A*(z+1)*d+B*(z+1)*s.diff(d,t),z+1)
    assert s.degree(s.gcd(d,s.diff(d,t)),t)==0
    bchecks.append({'divisor':str(s.expand(d)),'degree':int(s.degree(d,t)),'b':'z+1'})

# The entire deformation has one charge rule; bounded replay checks transcription.
P=W*(1+e*Z)-e*(1+s.Rational(3,2)*e*Z+s.Rational(1,2)*e**2*Z**2)*T**2
P1=P.subs(e,1)
eq(P,e*P1.subs({Z:e*Z,W:W/e},simultaneous=True))
for powers,cc in s.Poly(P,e,Z,W,T).terms():
    ee,zz,ww,tt=powers
    assert ee==zz-ww+1
cache={m:s.Poly(s.expand(P**m),Z,W,T).terms() for m in range(7)}
def moment(a,b,c,m):
    out=0
    for (zz,ww,tt),cc in cache[m]:
        if a+zz==b+ww and (c+tt)%2==0:
            out+=cc*s.factorial(a+zz)*s.factorial2(c+tt-1)
    return s.expand(out)
for m in range(1,7):
    eq(moment(0,0,0,m),0)
    eq(moment(1,0,0,m),s.factorial(m)*e**(m-1))
markers=[aa for aa in product(range(3),repeat=3) if sum(aa)<=2]
count=0
moments={}
for aa in markers:
    for m in range(1,7):
        v=moment(*aa,m);moments[(aa,m)]=v
        if v!=0:
            exponent=m+aa[1]-aa[0]
            assert exponent>=0
            assert not s.cancel(v/e**exponent).has(e)
        count+=1
ideals=[]
for N in range(6):
    forbidden=[v for (aa,m),v in moments.items() if m>N*sum(aa) and v!=0]
    gcd=s.Poly(forbidden[0],e)
    for v in forbidden[1:]:gcd=s.gcd(gcd,s.Poly(v,e))
    actual=s.expand(gcd.monic().as_expr());expected=1 if N==0 else e**N
    eq(actual,expected)
    ideals.append({'N':N,'finite_replay_gcd':str(actual),'all_order_ideal_from_written_charge_proof':str(expected)})
artin=[]
for length in range(1,7):
    witness=moment(1,0,0,length)
    assert s.rem(witness,e**length,e)!=0
    eq(witness,s.factorial(length)*e**(length-1))
    for (aa,m),v in moments.items():
        if m>length*sum(aa):eq(s.rem(v,e**length,e),0)
    artin.append({'length':length,'sharp_marker':'Z','nonzero_index':length,'exact_slope_from_written_proof':length})

receipt={'completed':True,'date':'2026-11-15','software':'SymPy '+s.__version__,
 'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
 'scalar_rationality_controls':receipts,'pole_isolation_controls':isolation,'squarefree_Bernstein_controls':bchecks,
 'rank_one_without_polynomial_seed_shortcut_detected':True,
 'Gaussian_marker_checks':count,'Gaussian_pure_checks':6,'Gaussian_marked_Z_checks':6,
 'charge_homogeneity_checked':True,'cutoff_ideal_controls':ideals,'Artin_sharpness_controls':artin,
 'general_rationality_module_solver_implemented':False,'general_rationality_module_solver_executed':False,
 'general_effective_field_scalar_solver_implemented':False,'universal_proofs_certified_by_checker':False,
 'latest_filtration_package_executed':False,'new_external_theorem_imported':False,
 'global_complex_safe_locus_uniformity_proved':False,'higher_actual_rank_classified':False,
 'new_SIC33_answer_computed':False,'new_Delta_minimum_computed':False,'numeric_cubic_frontier_changed':False,
 'operator_order_reduced':False,'independent_replay':False,'external_review':False,'formal_verification':False}
data=(json.dumps(receipt,indent=2,ensure_ascii=False,default=str)+'\n').encode()
out=p/'rank_one_image_verification.json'
if out.exists():assert out.read_bytes()==data
else:out.write_bytes(data)
print(json.dumps(receipt,indent=2,default=str))

