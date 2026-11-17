#!/usr/bin/env python3
"""Exact three-factor, nodal-conic and replica algebra. No analytic period theorem."""
from functools import lru_cache
from pathlib import Path
from math import factorial
import hashlib,json,time
import sympy as S


def moments(cov,nmax=8):
    @lru_cache(None)
    def w(n):
        if not any(n):return S.Integer(1)
        if sum(n)%2:return S.Integer(0)
        i=next(i for i,k in enumerate(n) if k);nn=list(n);nn[i]-=1;out=0
        for j,k in enumerate(nn):
            if k and cov[i][j]!=0:
                mm=nn.copy();mm[j]-=1;out+=k*cov[i][j]*w(tuple(mm))
        return S.expand(out)
    return [w((n,n,n)) for n in range(1,nmax+1)]


def pair_count_moment(n,cov):
    out=0
    for ij in range(n+1):
        for ik in range(n+1-ij):
            for jk in range(min(n-ij,n-ik)+1):
                rem=[n-ij-ik,n-ij-jk,n-ik-jk]
                if min(rem)<0 or any(k%2 for k in rem):continue
                loops=[k//2 for k in rem]
                den=S.prod(2**k*factorial(k) for k in loops)*factorial(ij)*factorial(ik)*factorial(jk)
                coeff=S.Rational(factorial(n)**3,den)
                out+=coeff*S.prod(cov[k][k]**loops[k] for k in range(3))*cov[0][1]**ij*cov[0][2]**ik*cov[1][2]**jk
    return S.expand(out)


def bezout(f,g,x):
    a,b,h=S.gcdex(f,g,x);assert h==1
    den=S.ilcm(*[S.denom(c) for poly in [a,b] for c in S.Poly(poly,x).all_coeffs()])
    a=S.expand(a*den);b=S.expand(b*den)
    assert S.expand(a*f+b*g)==den
    return {'variable':str(x),'f':str(f),'g':str(g),'a':str(a),'b':str(b),'constant':str(den)}


def product_checks():
    a,b,c,d,e,f=S.symbols('a b c d e f');cov=[[a,d,e],[d,b,f],[e,f,c]]
    generic=moments(cov)
    for n in [2,4,6,8]:assert S.expand(generic[n-1]-pair_count_moment(n,cov))==0
    assert all(generic[n-1]==0 for n in [1,3,5,7])
    s,q,r=S.symbols('s q r');normalized=[S.expand(g.subs({d:1,e:1,f:1})) for g in generic]
    symmetric=[]
    for n in [2,4,6,8]:
        form,rem,mapping=S.symmetrize(normalized[n-1],[a,b,c],formal=True);assert rem==0
        symmetric.append(form.subs({v:k for (v,_),k in zip(mapping,[s,q,r])}))
    assert symmetric[0]==2*s+r+8
    qval=(7*s*s+88*s+288)/20;rval=-2*s-8
    F=107*s**3+2196*s*s+11280*s+17152
    G=337*s**4+12736*s**3+109968*s*s+359424*s+406272
    red=[S.factor(v.subs(r,rval).subs(q,qval)) for v in symmetric]
    assert red[:2]==[0,0] and S.expand(red[2]+1620*F)==0 and S.expand(red[3]-381024*G)==0
    tri_cert=bezout(F,G,s)
    assert F.subs(s,-2)==2520 and S.degree(F,s)==3
    # All roots sigma of F give a nondegenerate Gram matrix after choosing
    # a,b,c to be roots of X^3-sigma X^2+qval X-rval.
    determinant=S.expand(S.det(S.Matrix([[a,1,1],[1,b,1],[1,1,c]])))
    ss,rem,mapping=S.symmetrize(determinant,[a,b,c],formal=True);assert rem==0
    ss=ss.subs({v:k for (v,_),k in zip(mapping,[s,q,r])}).subs(r,rval)
    assert S.expand(ss+3*(s+2))==0
    assert S.gcd(F,s+2)==1 and S.gcd(F,G)==1
    # Two off-diagonal edges, b normalized to1 whenever a leaf variance is nonzero.
    two=[S.expand(v.subs({d:1,e:1,f:0,b:1})).subs(a,-2*(1+c)/c) for v in generic]
    U=7*c*c-6*c+7;V=(c+1)*(2*c*c-17*c+2)
    assert S.cancel(two[1])==0 and S.cancel(two[3]+36*U)==0 and S.cancel(two[5]-32400*V)==0
    path_cert=bezout(U,V,c)
    assert U.subs(c,0)!=0 and U.subs(c,-1)!=0
    assert S.gcd(U,V)==1
    # One edge and no edges. These formulas check every remaining graph stratum.
    one=[S.expand(v.subs({d:1,e:0,f:0})) for v in generic]
    assert S.expand(one[1]-c*(a*b+2))==0
    assert S.expand(one[3]-3*c*c*(9*a*a*b*b+72*a*b+24))==0
    assert S.expand(one[3].subs({a:-2,b:1})+252*c*c)==0
    assert S.expand(generic[1].subs({d:0,e:0,f:0})-a*b*c)==0
    profiles=[{b:0,c:0,f:0},{c:0,e:0,f:0},{a:0,d:0,e:0,f:0}]
    for sub in profiles:assert all(S.expand(v.subs(sub))==0 for v in generic)
    # Exact radical's coordinate-prime decomposition in the six Gram entries.
    from itertools import combinations
    coords=[a,b,c,d,e,f]
    planes=[{a,d,e},{b,d,f},{c,e,f},{b,c,f},{a,c,e},{a,b,d}]
    minimal=[]
    for k in range(1,7):
        for selected in combinations(coords,k):
            support=set(selected)
            if all(support & plane for plane in planes) and not any(old<=support for old in minimal):minimal.append(support)
    expected=[{a,f},{b,e},{c,d},{a,b,c},{d,e,f}]
    assert {frozenset(v) for v in minimal}=={frozenset(v) for v in expected}
    for polynomial in generic:
        zz=S.Poly(polynomial,*coords,domain=S.ZZ)
        for exponents,coefficient in zz.terms():
            if coefficient==0:continue
            support={coords[i] for i,k in enumerate(exponents) if k}
            assert any(gen<=support for gen in expected)
    bad=set()
    for integer in [2,3,20,36,32400,1620,381024,252,int(tri_cert['constant']),int(path_cert['constant'])]:bad.update(S.factorint(integer))
    assert sorted(bad)==[2,3,5,7,11,23,79,107]
    # The all-order matching arguments are in the proof. Finite examples only
    # guard the profile boundary and the sharp mixed product Z^2 W.
    @lru_cache(None)
    def circular(i,j):
        if i<0 or j<0:return 0
        if i+j==0:return 1
        if i==0:return 0
        return j*circular(i-1,j-1)
    mixed=[]
    for k in range(1,7):
        for m in range(k+1,k+4):assert circular(k+m,2*m)==0
        value=circular(2*k,2*k);assert value==factorial(2*k)
        mixed.append({'degree_Q':k,'m':k,'mixed_value':value})
    return {'generic_Wick_vs_pair_counts':[2,4,6,8],
            'three_edge_symmetric_moments':list(map(str,symmetric)),
            'three_edge_elimination':{'q':str(qval),'r':str(rval),'F':str(F),'G':str(G),'determinant':str(ss)},
            'three_edge_Bezout':tri_cert,'two_edge_Bezout':path_cert,
            'sharp_nonzero_prefixes_by_edges':[2,4,6,8],
            'Gram_radical_generators':['a*f','b*e','c*d','a*b*c','d*e*f'],
            'coordinate_prime_components':6,'bad_characteristics_not_certified':sorted(bad),
            'scheme_boundary':'I8 is nonradical: a*f has Gram degree2 but every moment generator has degree at least3. No global nilpotence exponent computed.',
            'sharp_eight_witness':'sigma any complex root of F; a,b,c roots of X^3-sigma X^2+q(sigma)X-r(sigma); off-diagonals1; Gram determinant -3(sigma+2) nonzero; complex symmetric congruence realizes three linear forms.',
            'null_profile_regressions':3,'sharp_mixed_boundary':mixed}



def generic_thickness_checks():
    from itertools import product
    a,b,c,d,e,f,z=S.symbols('a b c d e f z');mm=moments([[a,d,e],[d,b,f],[e,f,c]])
    aa=-2*(d*d+e*e+4*z*d*e)/(1+2*z*z)
    row=[S.cancel(mm[k-1].subs({b:1,c:1,f:z}).subs(a,aa)) for k in [4,6,8]]
    RG=S.groebner(row[:2],d,e,domain=S.QQ.frac_field(z),order='grevlex')
    assert RG.reduce(row[2])[1]==0
    assert [S.Poly(v,d,e).total_degree() for v in row]==[4,6,8]
    star=[mm[k-1].subs({a:1,d:1,e:1}) for k in [2,4,6,8]]
    SG=S.groebner(star,b,c,f,domain=S.QQ,order='grevlex')
    def basis(G,variables):
        leading=[tuple(p.LM(order=G.order).exponents) for p in G.polys]
        bounds=[]
        for i in range(len(variables)):
            pure=[exp[i] for exp in leading if exp[i] and sum(exp)==exp[i]]
            assert pure,'No certified finite monomial bound'
            bounds.append(min(pure))
        base=[exp for exp in product(*(range(k) for k in bounds)) if not any(all(x>=y for x,y in zip(exp,lm)) for lm in leading)]
        return base
    rowbasis=basis(RG,[d,e]);starbasis=basis(SG,[b,c,f])
    assert len(rowbasis)==24 and len(starbasis)==6
    row_hilbert=[sum(sum(v)==n for v in rowbasis) for n in range(9)]
    assert row_hilbert==[1,2,3,4,4,4,3,2,1]
    assert max(map(sum,rowbasis))==8
    def mons3(n):
        return [b**i*c**j*f**(n-i-j) for i in range(n+1) for j in range(n+1-i)]
    assert all(SG.reduce(v)[1]==0 for v in mons3(4))
    assert any(SG.reduce(v)[1]!=0 for v in mons3(3))
    ranks=[]
    for k in range(5):
        vectors=[]
        for n in range(k,4):
            for mon in mons3(n):
                rem=S.Poly(SG.reduce(mon)[1],b,c,f)
                vectors.append([rem.coeff_monomial(exp) for exp in starbasis])
        ranks.append(S.Matrix.hstack(*(S.Matrix(v) for v in vectors)).rank() if vectors else 0)
    star_hilbert=[ranks[i]-ranks[i+1] for i in range(4)]
    assert ranks==[6,5,3,1,0] and star_hilbert==[1,2,2,1]
    # Semantics: inputs reduce to zero in the computed exact ideals.
    assert all(RG.reduce(v)[1]==0 for v in row) and all(SG.reduce(v)[1]==0 for v in star)
    cert={'field':'characteristic zero; row over QQ(z), star over QQ',
          'row':{'normalization':'b=c=1,f=z generic; a=-2*(d^2+e^2+4*z*d*e)/(1+2*z^2)',
                 'inputs':list(map(str,row)),'basis':list(map(str,RG)),
                 'standard_monomials':rowbasis,'length':24,'Hilbert_function':row_hilbert,'maximal_ideal_nilpotence_index':9,
                 'M8_reduces_by_M4_M6':True},
          'star':{'normalization':'a=d=e=1; variables b,c,f at0',
                  'inputs':list(map(str,star)),'basis':list(map(str,SG)),
                  'standard_monomials':starbasis,'length':6,'Hilbert_function':star_hilbert,
                  'maximal_ideal_nilpotence_index':4,'maximal_ideal_power_ranks':ranks},
          'top_dimensional_degree':3*24+3*6,
          'infinite_moment_boundary':'At a generic null row, reduced M_(2m) has degree2m and all degrees>=9 vanish. At a generic star, M_(2m) has maximal-ideal order>=m and the fourth power is zero. Thus higher moments do not change these generic local algebras. Embedded and boundary thickness not classified.'}
    return cert


def nodal_checks():
    u,v,T,h=S.symbols('u v T h');coeff=S.symbols('g00 g10 g01 g20 g11 g02 g30 g21 g12 g03')
    mons=[1,u,v,u*u,u*v,v*v,u**3,u*u*v,u*v*v,v**3]
    g=sum(a*m for a,m in zip(coeff,mons));H=u*v+h;P=T*H+g
    restriction=S.expand(g.subs(v,-h/u));tangent=S.expand(u*S.diff(g,u)-v*S.diff(g,v))
    assert S.expand(tangent.subs(v,-h/u)-u*S.diff(restriction,u))==0
    numerator=S.expand(u**3*u*S.diff(restriction,u))
    assert S.Poly(numerator,u).degree()<=6
    # Polynomial divisibility when restriction on irreducible H=0 is constant.
    alpha,beta,gamma,r=S.symbols('alpha beta gamma r');L=alpha*u+beta*v+gamma
    model=S.expand(r+H*L);assert S.expand(model.subs(v,-h/u)-r)==0
    assert S.div(model-r,H,u,v)[1]==0
    # Reducible fiber: two axes plus the center, with both center derivatives required.
    A=S.expand(g.subs(v,0));B=S.expand(g.subs(u,0))
    assert S.degree(S.diff(A,u),u)<=2 and S.degree(S.diff(B,v),v)<=2
    assert S.diff(g,u).subs({u:0,v:0})==coeff[1]
    assert S.diff(g,v).subs({u:0,v:0})==coeff[2]
    assert [S.diff(P,q).subs({u:0,v:0,h:0}) for q in [u,v,T]]==[coeff[1],coeff[2],0]
    # Exact affine-coordinate Jacobian for the moving conic chart.
    t,s=S.symbols('t s',nonzero=True);vp=(s/t-h)/u
    assert S.det(S.Matrix([[S.diff(u,u),S.diff(u,s)],[S.diff(vp,u),S.diff(vp,s)]]))==1/(t*u)
    assert S.cancel(t*H.subs(v,vp)-s)==0
    qcoeff=S.symbols('A B C D E F');Q=qcoeff[0]*u*u+qcoeff[1]*u*v+qcoeff[2]*v*v+qcoeff[3]*u+qcoeff[4]*v+qcoeff[5]
    phase=S.expand(s*s/2+t*g.subs(v,vp)-Q.subs(v,vp))
    powers=[term.as_powers_dict().get(u,0) for term in S.Add.make_args(phase)]
    assert min(powers)>=-3 and max(powers)<=3 and S.degree(phase,s)<=3
    # A 1/t-prefactor does not force a two-dimensional period to decay.
    width_phase=-s*s/(2*t*t)
    normalized_source=-1/t+S.diff(width_phase,t)
    twisted_s_of_s=1+s*S.diff(width_phase,s)
    assert S.cancel(t**3*normalized_source+t*t*twisted_s_of_s)==0
    # Null, finite, and nonisolated critical-locus examples, including the center trap.
    examples=[]
    for gg,hh,label in [(u+v,1,'smooth conic, two critical points'),(H*(u+2*v)+3,1,'smooth conic, critical curve'),(u**3+v**3,0,'both axis restrictions nonconstant but center critical line'),(u+v,0,'reducible fiber, no critical points')]:
        pp=S.expand(T*(u*v+hh)+gg.subs(h,hh));gb=S.groebner([S.diff(pp,q) for q in [T,u,v]],T,u,v)
        if label=='smooth conic, two critical points':assert gb.is_zero_dimensional
        elif label in ['smooth conic, critical curve','both axis restrictions nonconstant but center critical line']:assert not gb.is_zero_dimensional and list(gb)!=[1]
        else:assert list(gb)==[1]
        examples.append({'case':label,'gradient_Groebner_basis':list(map(str,gb))})
    return {'universal_tangent_identity':'PASS','smooth_conic_finite_bound':6,'reducible_off_center_bound':4,
            'center_line_condition':['g10=0','g01=0'],'chart_jacobian':'1/(t*u)',
            'Laurent_u_support':[-3,3],'s_degree_bound':3,'critical_locus_examples':examples,
            'prefactor_negative_control':'Phi=u+1/u-s^2/(2*t^2): (1/t) times its du/u ds period is a nonzero constant; exact normalized connection source reduces to zero.',
            'boundary':'The critical fiber u*v=0 and the u=0 chart complement are not removed from the period problem. No multidimensional logarithmic estimate is asserted.'}


def replica_checks():
    rho=S.symbols('rho');R=S.eye(4)
    for j in range(1,4):R[0,j]=R[j,0]=rho
    # Coordinates in each replica: Z,W,T, with E(Z_a W_b)=E(T_a T_b)=R_ab.
    @lru_cache(None)
    def w(n,perfect=False):
        if not any(n):return S.Integer(1)
        if sum(n)%2:return S.Integer(0)
        i=next(i for i,k in enumerate(n) if k);nn=list(n);nn[i]-=1;out=0
        ai,ki=divmod(i,3)
        for j,count in enumerate(nn):
            if not count:continue
            aj,kj=divmod(j,3)
            cc=1 if perfect else R[ai,aj]
            if (ki,kj) not in [(0,1),(1,0),(2,2)] or cc==0:continue
            mm=nn.copy();mm[j]-=1;out+=count*cc*w(tuple(mm),perfect)
        return S.expand(out)
    choices=[((0,1,0),S.Integer(1)),((1,0,1),S.Integer(1)),((3,0,0),-S.Rational(1,2))]
    from itertools import product
    star=0;perfect=0
    for selected in product(choices,repeat=4):
        n=tuple(k for m,c in selected for k in m);coef=S.prod(c for m,c in selected)
        star+=coef*w(n);perfect+=coef*w(n,True)
    assert S.expand(star+3*rho**3)==0 and S.expand(perfect)==0
    assert R.subs(rho,S.Rational(1,4)).det()==S.Rational(13,16)
    assert all(R.subs(rho,S.Rational(1,4))[:k,:k].det()>0 for k in range(1,5))
    # Formal conditional Gaussian identity: Z=t, pure exponent cancels exactly.
    t=S.symbols('t');assert -t*t**3/2+(t*t)**2/2==0
    return {'polynomial':'P=W+Z*T-Z^3/2','all_order_pure_identity':'exp(-t^4/2)*exp(t^4/2)=1, by circular translation and Gaussian integration',
            'star_covariance_four_replica_moment':str(S.expand(star)),
            'rho':'1/4','positive_covariance_principal_minors':'PASS','star_value':'-3/64',
            'perfect_correlation_moment':str(S.expand(perfect)),
            'boundary':'A complex GMC polynomial has all pure moments zero but a nonzero positive-definite correlated replica probe. This does not contradict real-coefficient orbit reconstruction and does not imply failure of GMC.'}


def main():
    start=time.monotonic();p=Path(__file__).resolve().parent
    products=product_checks();thickness=generic_thickness_checks();nodes=nodal_checks();replicas=replica_checks()
    cert={'three_edges':products['three_edge_Bezout'],'two_edges':products['two_edge_Bezout']}
    (p/'three_factor_bezout.json').write_text(json.dumps(cert,indent=2)+'\n')
    (p/'three_factor_thickness.json').write_text(json.dumps(thickness,indent=2)+'\n')
    receipt={'scope':'Exact finite/all-parameter algebra for the adjacent three-factor classification, nodal critical-locus reduction and replica-data obstruction. No new analytic period theorem.',
             'products':products,'generic_thickness':{'row_length':24,'row_nilpotence_index':9,'star_length':6,'star_nilpotence_index':4,'top_dimensional_degree':90,'field':'characteristic zero; row calculation over QQ(z), not a numeric sample'},'nodal':nodes,'replicas':replicas,'sympy_version':S.__version__,
             'elapsed_seconds':time.monotonic()-start,'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
             'certificate_sha256':hashlib.sha256((p/'three_factor_bezout.json').read_bytes()).hexdigest(),
             'thickness_certificate_sha256':hashlib.sha256((p/'three_factor_thickness.json').read_bytes()).hexdigest(),
             'assurance':'All-order product safety uses the written Wick partner count; sharpness uses algebraic covariance realization. No external review or formal verification.'}
    (p/'nodal_product_verification.json').write_text(json.dumps(receipt,indent=2)+'\n')
    print(json.dumps(receipt,indent=2))


if __name__=='__main__':main()
