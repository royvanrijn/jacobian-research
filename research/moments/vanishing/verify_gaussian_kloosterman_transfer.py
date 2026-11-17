"""Exact finite transcriptions of the all-k/s Gaussian Kloosterman proof.

No CAS Galois group, large source kernel, asymptotic numerical fit or old
sealed checker is run. Classical Galois imports remain written source imports.
"""
from pathlib import Path
from fractions import Fraction
import hashlib,json,math,sys
import sympy as s
import gaussian_kloosterman_transfer as K

p=Path(__file__).resolve().parent
sha=lambda f:hashlib.sha256(Path(f).read_bytes()).hexdigest()
dump=lambda x:json.dumps(x,indent=2)+'\n'
t,theta=K.t,K.theta

def powers(data,N):
    n=data['dimension'];P={tuple(e):Fraction(c) for e,c in data['phase_terms']}
    tables=[{(0,)*(n+1):Fraction(1)}]
    for m in range(N):
        out={}
        for a,c in tables[-1].items():
            for b,d in P.items():
                e=tuple(x+y for x,y in zip(a,b));out[e]=out.get(e,Fraction(0))+c*d
        tables.append({e:c for e,c in out.items() if c})
    return tables

def direct_moment(table,marker,k,copies=1):
    # Independent contraction from the ORIGINAL polynomial, not the transfer.
    out=s.S.Zero;block=k+3; e=s.Symbol('e')
    for exp,co in table.items():
        full=[a+b for a,b in zip(exp[:-1],marker)];value=Fraction(co)
        for j in range(copies):
            v=full[j*block:(j+1)*block]
            if v[0]!=v[1] or any(x%2 for x in v[2:]):value=0;break
            value*=math.factorial(v[0])
            for x in v[2:]:value*=math.prod(range(1,x,2))
        out+=s.Rational(value.numerator,value.denominator)*e**exp[-1] if value else 0
    return s.expand(out)

def jet(expr,N):return [s.expand(expr).coeff(t,m) for m in range(N+1)]

original_controls=[];moment_count=0;marker_count=0
for k in [2,3]:
    data=K.family(k);tables=powers(data,6);a=K.coefficients(k,6)
    for aa,bb,cc,ys in [(0,0,0,[0]*k),(1,0,0,[0]*k),(2,0,0,[0]*k),
                       (0,1,0,[0]*k),(0,0,1,[0]*k),(1,1,0,[0]*k),
                       (1,2,1,[0]*k),(1,0,0,[2]+[0]*(k-1)),(0,0,0,[1]+[0]*(k-1))]:
        marker=[aa,bb,2*cc]+ys;pol,op=K.marker_operator(aa,bb,cc,ys)
        values=[direct_moment(table,marker,k)/math.factorial(m) for m,table in enumerate(tables)]
        expected=[s.expand(pol).coeff(t,m)+op.subs(theta,m)*a[m] for m in range(7)]
        assert all(s.simplify(x-y)==0 for x,y in zip(values,expected)),(k,marker,values,expected)
        original_controls.append(dict(k=k,marker=marker,original_period_coefficients=[str(x) for x in values],
            polynomial_part=str(pol),Euler_operator=str(op)))
        moment_count+=7;marker_count+=1
    assert original_controls[-8]['original_period_coefficients']==[str(x) for x in [0]+a[1:]]

recurrences=[];gauge_controls=[];connection_controls=[]
z,q=s.symbols('z q')
for k in [2,3,4]:
    a=K.coefficients(k,9);A=K.quotient_operator(k);M=K.connection(k)
    assert a[0]==1 and all(a[m+1]==(2*m+1)**k*a[m] for m in range(9))
    v=s.Matrix([1]+[sum(m**j*a[m]*t**m for m in range(10)) for j in range(k)])
    residue=v.diff(t)*t-M*v
    assert all(s.expand(f).coeff(t,m)==0 for f in residue for m in range(9))
    assert s.simplify(s.trace(M[1:,1:])+s.Rational(k,2))==0
    gauge=s.expand(-z*A.subs({t:(-1)**k/(2**k*z),theta:-q-s.Rational(1,2)}))
    assert s.simplify(gauge-(q**k-z))==0
    recurrences.append(dict(k=k,coefficients=[str(x) for x in a],operator=str(A)))
    gauge_controls.append(dict(k=k,transformed_operator=str(gauge),equals='q^k-z'))
    connection_controls.append(dict(k=k,matrix=[[str(x) for x in row] for row in M.tolist()],galois=K.galois_model(k)))

# Exact Leibniz verification of the symbolic symmetric connection, without
# inferring actual generic rank from a finite series prefix.
symmetric=[]
for k,m in [(2,2),(2,3),(3,2)]:
    basis,M=K.symmetric_connection(k,m);base=K.connection(k);v=s.symbols('v:'+str(k+1))
    monomials=s.Matrix([s.prod(x**j for x,j in zip(v,e)) for e in basis]);rhs=M*monomials
    for i,f in enumerate(monomials):
        lhs=sum(s.diff(f,v[j])*sum(base[j,l]*v[l] for l in range(k+1)) for j in range(k+1))
        assert s.expand(lhs-rhs[i])==0
    assert len(basis)==math.comb(k+m,m)
    summed=K.independent_sum_family(k,m)
    assert summed['dimension']==m*(k+3) and len(summed['phase_terms'])==5*m
    symmetric.append(dict(k=k,copies=m,rank_formula=math.comb(k+m,m),basis=[list(e) for e in basis],
        matrix=[[str(x) for x in row] for row in M.tolist()],galois=K.galois_model(k,m),
        actual_rank_uses_written_translation_independence_proof=True))

# Finite ordinary ideal controls include independent sums and preserve e.
profiles=[]
for copies in [1,2]:
    data=K.independent_sum_family(2,copies,parameterized=True);tables=powers(data,3);n=data['dimension'];e=s.Symbol('e')
    for m in range(1,4):
        for d in range(3):
            vals=[]
            for total in range(d+1):
                for alpha in K.compositions(n,total):
                    value=direct_moment(tables[m],alpha,2,copies)
                    if value:vals.append(s.Poly(value,e))
            g=s.Poly(0,e)
            for value in vals:g=s.gcd(g,value)
            if g:g=g.monic()
            expected=0 if d==0 else e**max(m-d,0)
            assert g.as_expr()==expected,(copies,m,d,g)
            profiles.append(dict(copies=copies,power=m,marker_degree=d,ordinary_ideal_generator=str(g.as_expr())))

# Same-rank origin formal-type comparison: exact chain/gauge identities.
w=s.Symbol('w',positive=True);u=s.Function('u')(w)
def Th(f):return -w*s.diff(f,w)/2
def old(f):return s.simplify(f-(4*Th(4*Th(f)+3*f)+(4*Th(f)+3*f))/(4*w**2))
def new(f):return s.simplify(f-(2*Th(2*Th(f)+f)+(2*Th(f)+f))/w**2)
assert s.simplify(old(s.sqrt(w)*u)+s.sqrt(w)*(s.diff(u,w,2)-u))==0
assert s.simplify(new(w*u)+w*(s.diff(u,w,2)+s.diff(u,w)/w-u))==0
for sign in [1,-1]:assert old(s.sqrt(w)*s.exp(sign*w))==0
formal_units=[]
for sign in [1,-1]:
    cs=[s.S.One]
    for j in range(5):cs.append(s.cancel(cs[-1]*(s.Rational(1,2)+j)**2/(2*sign*(j+1))))
    f=sum(c*w**(-j) for j,c in enumerate(cs));res=s.expand(s.diff(f,w,2)+2*sign*s.diff(f,w)+f/(4*w**2))
    assert all(res.coeff(w,-j)==0 for j in range(2,7))
    formal_units.append(dict(sign=sign,coefficients=[str(c) for c in cs]))
comparison=dict(old_operator='1-t*(4*theta+1)*(4*theta+3)',new_operator='1-4*t*(2*theta+1)^2',
    Gaussian_dimension=5,phase_degree=8,actual_rank=3,exact_Gevrey_order=2,
    common_ramification='t=1/(4*w^2)',common_quotient_formal_types=['w^(1/2)*exp(w)','w^(1/2)*exp(-w)'],
    old_identity_component_solvable=True,new_full_group='Ga^2 semidirect SL_2',
    new_formal_units=formal_units,all_global_formal_types_equal_claimed=False,Stokes_matrices_computed=False)

rejections=[]
for name,call in [('zero k',lambda:K.family(0)),('zero scale',lambda:K.family(2,0)),
                  ('zero copies',lambda:K.independent_sum_family(2,0)),
                  ('nonintegral marker',lambda:K.marker_operator(1,0,0,[s.Rational(1,2),0])),
                  ('rank-one classical-group misuse',lambda:K.galois_model(1))]:
    try:call()
    except ValueError:rejections.append(name)
    else:raise AssertionError(name)

packet=dict(original_marker_controls=original_controls,recurrences=recurrences,gauge_controls=gauge_controls,
    connection_controls=connection_controls,symmetric_controls=symmetric,ordinary_profiles=profiles,comparison=comparison,
    rejections=rejections,universal_claims_use_written_proofs=True,Galois_groups_computationally_proved=False,
    classical_Galois_import='Katz theorem,explicit primary restatement in Fresan-Sabbah-Yu1810.06454v5 Proposition2.7 proof',
    regular_transcendental_Gaussian_realization_claimed=False,all_ordinary_ideals_retain_nilpotents=True,
    independent_replay=False,formal_verification=False,external_review=False)
(p/'gaussian_kloosterman_certificate.json').write_text(dump(packet))
runtime=dict(python=sys.version,sympy=s.__version__,sealed_imports_sha256={name:sha(p/name) for name in [
    'SQUARED_SPECTATOR_INVARIANCE_AND_UNBOUNDED_GAUSSIAN_IMAGES.md','squared_spectator_geometry_source_audit.json',
    'FINITE_MULTIPLIER_MIXED_DECISION.md','finite_multiplier_source_audit.json','finite_multiplier_verification.json']},
    CAS_Galois_algorithm_run=False,old_UNKNOWN_jobs_restarted=False)
(p/'gaussian_kloosterman_runtime.json').write_text(dump(runtime))
receipt=dict(completed=True,canonical_proof_sha256=sha(p/'GAUSSIAN_KLOOSTERMAN_QUOTIENTS_AND_FAITHFUL_SYMMETRIC_IMAGES.md'),
    producer_sha256=sha(p/'gaussian_kloosterman_transfer.py'),checker_sha256=sha(Path(__file__)),
    certificate_sha256=sha(p/'gaussian_kloosterman_certificate.json'),runtime_sha256=sha(p/'gaussian_kloosterman_runtime.json'),
    original_Wick_moments=moment_count,full_marker_transfers=marker_count,recurrence_controls=len(recurrences),
    Kloosterman_gauge_controls=len(gauge_controls),symmetric_Leibniz_controls=len(symmetric),ordinary_profile_controls=len(profiles),
    origin_formal_comparison=True,rejections=len(rejections),universal_claims_use_written_proofs=True,
    Galois_groups_computationally_proved=False,actual_image_kernel_computed=False,
    independent_replay=False,formal_verification=False,external_review=False)
(p/'gaussian_kloosterman_verification.json').write_text(dump(receipt))
print(json.dumps(receipt))
