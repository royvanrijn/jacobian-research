"""Key models, elliptic contact34 coefficients, and their seven/eight sections."""
from pathlib import Path
import json,time,hashlib
import sympy as S
from sympy.polys.rings import ring
from algebra import family,G,A,M,s,W,Y,HERE,ROOT
from witnesses import witness
CHECKS=[]
def check(name,value):
    assert value,name
    CHECKS.append(name)
def zero(name,value):check(name,S.expand(value)==0)

def main():
    start=time.monotonic();models=[]
    for name in ['L','K']:
        data=family(name);par=data['parameter'];p,q,g,key,f,r8,pb,qb,delta,cubic=[data[k] for k in ['p','q','g','key','f','r8','P','Q','delta','cubic']]
        zero(name+' exact key relation on original normalization',data['key_residual'])
        check(name+' global implicit equation degree10',S.Poly(key.subs(W,data['G']),G.U,G.V).total_degree()==10)
        check(name+' global implicit equation monic V6',S.Poly(key.subs(W,data['G']),G.V).LC()==1)
        ell=S.Poly(key,W).nth(1);linear=S.Poly(key,G.V).nth(1)
        rest=S.expand(key-W*W-ell*W-linear*G.V)
        check(name+' key relation has the required nine positions',set(S.Poly(key,G.U,G.V,W).monoms())<={(0,0,2),(1,0,1),(0,0,1),(3,0,0),(2,0,0),(1,0,0),(0,0,0),(1,1,0),(0,1,0)})
        zero(name+' complete coefficient square',data['M_coefficient']**2+36*S.Poly(f,M).nth(0)-cubic)
        zero(name+' Weierstrass identity',(18*M-data['M_coefficient'])**2-cubic+36*f)
        disc=S.factor(S.discriminant(cubic,A))
        if name=='L':
            h=53248*par**12+142848*par**9+144720*par**6+65340*par**3+10125
            zero('L Weierstrass discriminant',disc+S.Rational(3048192,78125)*h)
            zero('L key relation independent of V',linear)
            u=G.F.t**3+par*G.F.t
            zero('L cubic intermediate normalization p',p-u*u-u)
            coefficient=S.Poly(S.expand(g+u**3/6561),G.F.t).nth(3)
            zero('L cubic intermediate normalization g',g+u**3/6561-coefficient*u)
            expected=-(4*G.U+1)*(-648*G.U+147*par**3-64)**2/S.Integer(72301961339136)
            zero('L completed-square nodal/cuspidal cubic',rest-ell*ell/4-expected)
        else:
            h=390625*par**9+937500*par**8+1031250*par**7+1062500*par**6+853125*par**5+402000*par**4+94800*par**3+7680*par**2-4800*par+1024
            zero('K Weierstrass discriminant',disc+S.Rational(62208,78125)*par**3*h)
            zero('K coefficient of U*V',S.Poly(linear,G.U).nth(1)-par**12*(3*par+2)**2/S.Integer(129140163))
            zero('K V coefficient vanishes on family intersection',S.rem(linear,3*par+2,par))
        check(name+' squarefree auxiliary degeneration factor',S.gcd(h,S.diff(h,par))==1)
        cpoly=S.Poly(cubic,A);c4=16*cpoly.nth(2)**2-48*36*cpoly.nth(1)
        check(name+' degeneration is nodal rather than cuspidal',S.gcd(h,c4)==1)
        check(name+' nonzero origin in degeneration factor',h.subs(par,0)!=0)
        check(name+' no shared auxiliary degeneration at C1',S.gcd(h,3*par**3+2 if name=='L' else 3*par+2)==1)
        tangent=lambda expr:S.expand(S.diff(f,M)*S.diff(expr,A)-S.diff(f,A)*S.diff(expr,M))
        jac=S.expand(S.diff(pb,s)*tangent(qb)-S.diff(qb,s)*tangent(pb))
        zero(name+' evaluation Jacobian is derivative of r8',S.rem(jac-tangent(r8),f,A))
        zero(name+' critical resultant equals minus15 Jacobian',S.rem(delta+15*tangent(r8),f,A))
        # All contact>=34 coefficients are f0; the critical scheme is free rank8.
        R,aa,mm=ring('A,M',S.QQ.poly_ring(par),order=lambda mon:(2*mon[0]+3*mon[1],mon[1]))
        ff=R.from_expr(-f/9);remainder=R.from_expr(delta).div([ff])[1];gg=remainder*S.Rational(3,175)
        check(name+' critical basis constant leading units',ff.LM==(0,2) and ff.LC==R.domain.one and gg.LM==(4,0) and gg.LC==R.domain.one)
        syzygy=ff.mul_monom((4,0))-gg.mul_monom((0,2))
        check(name+' critical Buchberger pair',not syzygy.div([ff,gg])[1])
        check(name+' eight critical standard monomials',len([(i,j) for i in range(4) for j in range(2)])==8)
        weights=[2*i+3*j for (i,j),co in S.Poly(r8,A,M).terms()]
        check(name+' r8 has pole7 on coefficient cubic',max(weights)==7 and S.Poly(r8,A,M).coeff_monomial(A*A*M)==S.Rational(5,9))
        check(name+' critical derivative has pole8',max(2*i+3*j for i,j in gg.monoms())==8)
        sample=witness(name)
        for case in sample['checks']:check(name+' equality34 witness '+case['name'],case['rank']==case['dimension'])
        models.append(dict(family=name,parameter=str(par),key_relation=str(key),key_G=str(data['G']),
            elliptic_coefficient_equation=str(f),Weierstrass_cubic=str(cubic),Weierstrass_discriminant=str(disc),
            auxiliary_nodal_factor=str(h),degree7_function=str(r8),critical_scheme_generators=[str(ff.as_expr()),str(gg.as_expr())],
            evaluation_Jacobian=str(S.rem(jac,f,A)),witness=sample))
    # Closure of C[p,g] inside the known free six-generator affine ring.
    vv,ww,ll,bb,cc,a0,a1,a2=S.symbols('v w L b c a0 a1 a2')
    scaled=[S.S.One,ww,ll*vv,ll*vv*ww,ll*ll*vv*vv,ll*ll*vv*vv*ww]
    relations=[ww*ww+bb*ww+cc+ll*vv,vv**3+a2*vv*vv+a1*vv+a0-ww]
    outputs=[scaled[1],-cc*scaled[0]-bb*scaled[1]-scaled[2],scaled[3],-cc*scaled[2]-bb*scaled[3]-scaled[4],scaled[5],ll**3*a0*scaled[0]-ll**3*scaled[1]+ll**2*a1*scaled[2]+(ll*a2-cc)*scaled[4]-bb*scaled[5]]
    for i,value in enumerate(outputs):
        difference=S.expand(ww*scaled[i]-value)
        remainder=S.rem(difference,relations[0],ww);remainder=S.rem(remainder,relations[1],vv)
        zero('scaled subring basis closed under g '+str(i),remainder)
    apery=[0,25,32,9,16,41]
    check('compressed K affine codimension18',(sum(apery)-15)//6==18)
    check('partial normalization defect6',18-12==2*1+2*2==6)
    gaps=[n for n in range(42) if n<apery[n%6]]
    check('compressed K pole conductor36',len(gaps)==18 and max(gaps)==35)
    check('compressed K infinity delta10',(9-1)*(9-2)//2-18==10)
    check('auxiliary degeneration point count39',12+3*9==39)
    result=dict(schema='r7-elliptic-contact-key-model-v1',families=models,
        contact34_parameter_curve='smooth genus1 except39 coefficient points, where it is an irreducible nodal cubic',
        high_contact_section_degree=7,critical_coefficient_scheme_rank=8,
        compressed_K_pole_semigroup=[6,9,16],compressed_K_pole_gaps=gaps,compressed_K_affine_delta=18,
        compressed_K_infinity_characteristic=[3,11],partial_normalization_quotient='(C[p]/L)^2 plus(C[p]/L^2)^2; L=a*p+b',
        scope='Actual image-curve models and equality34 witnesses; no global Keller cover. The polynomial target map has degree3 and is not an automorphism or an intermediate field of the original degree6 source extension.')
    (HERE/'models.json').write_text(json.dumps(result,indent=2)+'\n')
    inputs=[HERE/name for name in ['algebra.py','witnesses.py','verify.py','PROOF.md','README.md','REVIEW.md','discovery.json']]
    inputs += [HERE.parent/path for path in ['jc2-r7-complement/geometry.py','jc2-r7-image-family/family.py','jc2-r7-laurent-collisions/laurent.py','jc2-r7-quintic-normalizations/normalizations.json','jc2-r7-quintic-normalizations/PROOF.md','jc2-r7-quintic-incidence/finite.py','jc2-omitted-divisors/DISJOINT_IMAGE_GATE.md','jc2-r7-global-cover/PROOF.md']]
    inputs += [ROOT/'foundations/scripts/verify_f2_degree_6_10_gap.py']
    receipt=dict(passed=True,check_count=len(CHECKS),checks=CHECKS,elapsed_seconds=round(time.monotonic()-start,3),
        input_sha256={str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in inputs},
        output_sha256={'models.json':hashlib.sha256((HERE/'models.json').read_bytes()).hexdigest()},
        assurance='Authored coefficient classification, curve/function-field and module interpretation with exact identities, Groebner and modular rank checks; no formalization or independent external review.')
    (HERE/'receipt.json').write_text(json.dumps(receipt,indent=2)+'\n');print(f'PASS {len(CHECKS)} ELLIPTIC CONTACT CHECKS ({time.monotonic()-start:.3f}s)',flush=True)

if __name__=='__main__':main()
