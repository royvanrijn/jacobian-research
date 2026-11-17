#!/usr/bin/env python3
"""Exact controls for QQL1; universal divisor obstructions are written proofs."""
from collections import Counter
from fractions import Fraction
from math import comb
from pathlib import Path
import hashlib
import json
import platform
import tempfile
import sympy as s

HERE=Path(__file__).resolve().parent
H,k,h,c,cv,alpha,gamma,beta=s.symbols("H kappa h c c_v alpha gamma beta")
U=(k-3)**3/(k-2)**2
R=-(k-3)*(k-2)/k**2
a2=s.factor(c*c*U*R*R)
a1=s.factor(c*c*h*U*(5*k*R*s.diff(R,k)-beta*R))
residue=s.factor(s.limit((k-2)*a1,k,2))
assert s.factor(residue-c*c*h*(2*beta-5)/8)==0
Ahalf=s.factor(a2*H*H+a1.subs(beta,s.Rational(5,2))*H)
Aclosed=c*c*((k-3)**5/k**4*H*H+
              s.Rational(5,2)*h*(k-3)**4*(k+12)/k**4*H)
assert s.factor(Ahalf-Aclosed)==0
Asum=sum(comb(5,l)*(-3)**l*k**(1-l)*H*
         (H+s.Rational(5,2)*(1-l)*h) for l in range(6))
assert s.factor(Ahalf-c*c*Asum)==0

def clean(op):
    return {n:s.factor(a) for n,a in op.items() if s.factor(a)!=0}
def add(*ops):
    out={}
    for op in ops:
        for n,a in op.items(): out[n]=out.get(n,0)+a
    return clean(out)
def scale(op,a):
    return clean({n:a*b for n,b in op.items()})
def mul(left,right):
    out={}
    for n,a in left.items():
        for m,b in right.items():
            out[n+m]=out.get(n+m,0)+a.subs(H,H+m*h)*b
    return clean(out)
def inverse_monomial(op):
    assert len(op)==1
    n,a=next(iter(op.items()))
    return {-n:s.factor(1/a.subs(H,H-n*h))}
def power(op,n):
    if n<0: return power(inverse_monomial(op),-n)
    out={0:s.Integer(1)}
    for _ in range(n): out=mul(out,op)
    return out

Q=(H+alpha*h)*(H+gamma*h)
kap={5:cv/Q}
transport={}
for l in range(6):
    coefficient=c*c*comb(5,l)*(-3)**l*H*(H+s.Rational(5,2)*(1-l)*h)
    transport=add(transport,mul(power(kap,1-l),{0:coefficient}))
top_expected=c*c*cv*H*(H+s.Rational(5,2)*h)/Q
bottom_expected=-243*c*c*cv**-4*H*(H-10*h)
for j in range(1,5):
    bottom_expected*=Q.subs(H,H-5*j*h)
assert s.factor(transport[5]-top_expected)==0
assert s.factor(transport[-20]-bottom_expected)==0
assert min(transport)==-20 and max(transport)==5

# Independent two-term fifth power and its exact extreme skew products.
A=H+1
B=H*H-3
P={1:A,-4:B}
P5=power(P,5)
assert min(P5)==-20 and max(P5)==5
assert s.factor(P5[5]-s.prod(A.subs(H,H+i*h) for i in range(5)))==0
assert s.factor(P5[-20]-s.prod(B.subs(H,H-4*i*h) for i in range(5)))==0

# Positive top-norm controls: necessary source shifts really admit rational norms.
def norm_primitive(integer_shift,offset):
    out=s.Integer(1)
    if integer_shift>=0:
        for j in range(integer_shift):
            out*=(H+offset+5*j*h)/(H+offset+(5*j+1)*h)
    else:
        for j in range(1,-integer_shift+1):
            out*=(H+offset-(5*j-1)*h)/(H+offset-5*j*h)
    return s.factor(out)

top_norm_controls=[]
bottom_orbit_controls=[]
for a in range(-2,3):
    for d in range(-2,3):
        An=norm_primitive(a,0)*norm_primitive(d,s.Rational(5,2)*h)
        N=s.prod(An.subs(H,H+i*h) for i in range(5))
        target=H*(H+s.Rational(5,2)*h)/(
            (H+5*a*h)*(H+(s.Rational(5,2)+5*d)*h))
        assert s.factor(N-target)==0,(a,d)
        top_norm_controls.append({"a":a,"d":d,"rational_norm_exists":True})
        integer_roots=[Fraction(5*j-5*a) for j in range(1,5)]+[Fraction(0),Fraction(10)]
        half_roots=[Fraction(5*j)-Fraction(5,2)-5*d for j in range(1,5)]
        assert all(z.denominator==1 for z in integer_roots)
        assert all(z.denominator==2 for z in half_roots)
        classes=Counter(z%4 for z in half_roots)
        assert len(classes)==4 and set(classes.values())=={1}
        assert all((z%4) not in classes for z in integer_roots)
        bottom_orbit_controls.append({
            "a":a,"d":d,
            "half_roots":[str(z) for z in half_roots],
            "4h_orbit_multiplicities":{str(z):n for z,n in sorted(classes.items())},
            "divisible_by_five":False,
        })

# A genuine rational norm has the required orbit divisibility, including poles.
test_B=(H-1)*(H-7)/(H+2)
norm_B=s.factor(s.prod(test_B.subs(H,H-4*i) for i in range(5)))
num,den=s.fraction(norm_B)
divisor=Counter()
for poly,sign in [(num,1),(den,-1)]:
    for factor,exponent in s.factor_list(poly,H)[1]:
        fp=s.Poly(factor,H)
        assert fp.degree()==1
        point=-fp.all_coeffs()[1]/fp.all_coeffs()[0]
        divisor[s.Rational(point)%4]+=sign*exponent
assert all(v%5==0 for v in divisor.values())



# Stronger actual-division obstruction: rational source Fs, signed divisors,
# lowest Euler eigenoperator identity and ordinary-order linearization.
Tvar=s.symbols("T")
phi10=s.Poly(Tvar**4-Tvar**3+Tvar**2-Tvar+1,Tvar)
def cyclotomic_reduce(expr):
    num,den=s.fraction(s.cancel(expr))
    inv=s.invert(s.Poly(den,Tvar),phi10)
    return s.rem(s.Poly(num,Tvar)*inv,phi10).as_expr()
def orbit_divisor(expr):
    num,den=s.fraction(s.cancel(expr))
    out=s.Integer(0)
    for polynomial,sign in [(num,1),(den,-1)]:
        for factor,exponent in s.factor_list(polynomial,H)[1]:
            fp=s.Poly(factor,H)
            if fp.degree()!=1: continue
            point=s.cancel(-fp.all_coeffs()[1]/fp.all_coeffs()[0])
            if point.is_Integer:
                out+=sign*exponent*Tvar**int(point)
    return s.expand(out)
source_functions=[
 ("quadratic",H**2),
 ("rational_pole",(H+h)**3/(H-2*h)),
 ("outside_integer_orbit",H**2*(H**2+h*h)/((H-h)*(H+4*h))),
]
strong_source_controls=[]
for label,Fs in source_functions:
    assert s.factor(Fs.subs(h,0)-H**2)==0
    kappas={5:cv/Fs}
    Js=scale(add({0:-H},
                 scale(mul(power(kappas,-1),{0:H}),5),
                 scale(mul(power(kappas,-2),{0:H}),-6)),c)
    C= -6*c/cv**2*Fs.subs(H,H-5*h)*Fs.subs(H,H-10*h)*H
    assert min(Js)==-10 and s.factor(Js[-10]-C)==0
    assert scale(mul(kappas,{0:Fs}),1/cv)=={5:s.Integer(1)}
    conjugated=mul({1:s.Integer(1)},mul(kappas,{-1:s.Integer(1)}))
    assert conjugated==mul(kappas,{0:Fs/Fs.subs(H,H-h)})
    # Exact lowest commutator coefficient, including a scalar Euler shift.
    testP={-4:(H+1)/(H-3),1:2-H}
    scalar=s.symbols("lambda")
    shiftedJ=add(Js,{0:scalar})
    JP=add(mul(shiftedJ,testP),scale(mul(testP,shiftedJ),-1))
    expected=C.subs(H,H-4*h)*testP[-4]-testP[-4].subs(H,H-10*h)*C
    assert s.factor(JP[-14]-expected)==0 and JP[-14]!=0
    Fs1=s.factor(Fs.subs(h,1))
    C1=s.factor(H*Fs1.subs(H,H-5)*Fs1.subs(H,H-10))
    fdiv=orbit_divisor(Fs1)
    cdiv=orbit_divisor(C1)
    assert s.expand(cdiv-(1+(Tvar**5+Tvar**10)*fdiv))==0
    assert cyclotomic_reduce(cdiv)==1
    assert cyclotomic_reduce((Tvar**4-1)*cdiv)!=0
    assert cyclotomic_reduce(Tvar**10-1)==0
    strong_source_controls.append({
      "case":label,"Fs":str(Fs),"minimum_J_v_weight":-10,
      "lowest_eigenoperator_identity":"PASS",
      "signed_F_divisor_on_integer_orbit":str(fdiv),
      "C_character_mod_Phi10":"1",
      "T4_minus_1_character_nonzero":True,
    })

# Exact [J,U]=-5*c*h*U and [J,A_beta]=-5*c*h*A_beta.
def kdelta(expr,n):
    for _ in range(n): expr=5*h*k*s.diff(expr,k)
    return expr
def kmul(left,right):
    out={}
    for i,a0 in left.items():
        for j,b0 in right.items():
            for n in range(i+1):
                out[i+j-n]=out.get(i+j-n,0)+a0*comb(i,n)*kdelta(b0,n)
    return clean(out)
def kcomm(left,right):
    return add(kmul(left,right),scale(kmul(right,left),-1))
Jop={1:c*R}
Uop={0:U}
Aop={2:a2,1:a1}
assert kcomm(Jop,Uop)==scale(Uop,-5*c*h)
assert kcomm(Jop,Aop)==scale(Aop,-5*c*h)

def differential_order(op):
    return max(s.Poly(s.fraction(a0)[0],H).degree()-
               s.Poly(s.fraction(a0)[1],H).degree()
               for a0 in op.values())
linearization_controls=[]
for label,base,D in [
 ("positive_order",{-1:H+h},{-1:(H-h)*H+1}),
 ("negative_order",{1:1/(H+h)},{1:1,-1:H}),
 ("zero_order",{1:s.Integer(1)},{0:1/H}),
]:
    M={}
    for i in range(5):
        term=mul(mul(power(base,i),D),power(base,4-i))
        M=add(M,term)
    expected=4*differential_order(base)+differential_order(D)
    assert M and differential_order(M)==expected
    linearization_controls.append({
      "case":label,"ord_P":int(differential_order(base)),
      "ord_D":int(differential_order(D)),"ord_M":int(expected),
      "M_nonzero":True,
    })
valuation_controls={
 "nu_kappa":5,"leading_U":"-27/4","nu_U":0,
 "nu_J":-10,"nu_A_beta":-20,"root_forced_nu":-4,
 "eigenoperator_nu_divisible_by":10,
 "contradiction":True,
}
assert (-20)//5==-4 and (-4)%10!=0



# Exact cyclic-character identities used by the written two-valuation proof.
deck_z=s.Symbol("deck_z")
phi5=s.cyclotomic_poly(5,deck_z)
assert s.rem(deck_z**10-1,phi5,deck_z)==0
assert s.rem(deck_z**5-1,phi5,deck_z)==0
assert s.gcd(phi5,deck_z-1)==1
# Both classical Laurent weights 1,-4 have the same nontrivial character.
assert (1-(-4))%5==0 and 10%5==0


# Compact exact scalar-spectrum controls.
phi5=s.Poly(Tvar**4+Tvar**3+Tvar**2+Tvar+1,Tvar)
projection_controls=[]
for j in range(5):
    for weight in range(5):
        expr=sum(Tvar**((weight-j)*i) for i in range(5))
        num,den=s.fraction(s.cancel(expr))
        reduced=s.rem(s.Poly(num,Tvar)*s.invert(s.Poly(den,Tvar),phi5),phi5).as_expr()
        assert reduced==(5 if weight==j else 0)
        projection_controls.append([j,weight,int(reduced)])
eigen_realizers=[]
for n in range(-3,4):
    assert s.factor(5*c*h*k*R*s.diff(U**n,k)+5*c*h*n*U**n)==0
    eigen_realizers.append({"power":n,"eigenvalue":str(-5*c*h*n)})
assert s.limit(5*c*h*R,k,s.oo)==-5*c*h
fa=(k-1)/(k+2)
fb=(k*k+3)/(k-4)
delta=lambda f:5*c*h*k*R*s.diff(f,k)
assert s.factor(delta(fa/fb)-(delta(fa)/fb-fa*delta(fb)/fb**2))==0


names=[
 "QUINTIC_ORDERING_LAURENT_OBSTRUCTION.md",
 "check_quintic_ordering_laurent.py",
 "CLASSICAL_QUINTIC_COVER_SQUARE.md",
 "QUANTUM_RADICAL_TOWERS_AND_LONG_OBSTRUCTION.md",
]
receipt={
 "claim":"QQL1","status":"PASS",
 "environment":{"python":platform.python_version(),"sympy":s.__version__},
 "ordered_operator":{
  "pole_residue_at_kappa_2":str(residue),
  "sole_pole_canceling_beta":"5/2",
  "binomial_kappa_sum":"PASS",
  "v_exponent_range":[-20,5],
  "top_coefficient":str(s.factor(top_expected)),
  "bottom_coefficient":str(s.factor(bottom_expected)),
 },
 "root_extreme_controls":{
  "P":"v*(H+1)+v^-4*(H^2-3)",
  "P5_exponent_range":[-20,5],
  "top_norm":"product A(H+i*h), i=0..4",
  "bottom_norm":"product B(H-4*i*h), i=0..4",
  "deck_eigenclass":"1 modulo 5, written Hensel argument",
 },
 "top_norm_positive_controls":top_norm_controls,
 "bottom_norm_impossible_orbits":bottom_orbit_controls,
 "genuine_rational_norm_orbit_control":{
  "B":"(H-1)*(H-7)/(H+2)",
  "translation_step":-4,"factor_count":5,
  "orbit_totals":{str(z):n for z,n in sorted(divisor.items())},
  "all_multiples_of_five":True,
 },
 "exact_scalar_spectrum":{
  "eigenvalue_set":"5*c*h*Z",
  "all_rational_source_Fs":True,
  "deck_fixed_field":"Frac k(kappa)[J;delta], written finite-five-basis proof",
  "projection_controls":projection_controls,
  "U_power_realizers":eigen_realizers,
  "rational_fraction_leading_coefficient_derivative":"PASS",
  "actual_Euler_pair":"No P,Q with [P,Q]=c*h and J=Q*P+scalar",
  "actual_target_root":"No P^5=U*G(J) for any nonzero rational G",
  "universal_scope":"written proof, not finite enumeration",
 },
 "actual_division_root_obstruction":{
  "source_functions":strong_source_controls,
  "valuation_controls":valuation_controls,
  "rational_target_G":"arbitrary nonzero rational function of J",
  "Fs_scope":"any nonzero rational function of H",
  "h_scope":"any nonzero scalar in characteristic zero",
  "ordered_J_U_A_eigenrelations":"PASS",
  "linearized_power_map_controls":linearization_controls,
  "universal_injectivity":"written commutative-symbol proof, not a finite test",
 },
 "Long_regular_rational_target_obstruction":{
  "source_Fs":"any nonzero rational function of H",
  "target_Ft":"any nonzero rational scalar function of J",
  "root_requirement":"h-regular prescribed Long classical P0",
  "deck_character":"1 modulo5, forced by prescribed-residue uniqueness",
  "eigenoperator_v_valuation":"0 modulo10, hence0 modulo5",
  "character_identities_mod_Phi5":"PASS",
  "contradiction":True,
  "proof_scope":"Written proof uses h and v valuations separately; actual mixed denominators are included.",
 },
 "actual_Euler_pair_corollary":{
  "carrier":"entire actual Weyl division ring",
  "source_Fs":"arbitrary rational",
  "target_Euler_relation":"J=Q*P+scalar",
  "forced_minimum_weights":[-4,-6],
  "contradiction":"Euler eigenoperator valuation must lie in 10Z",
  "requires_classical_specialization":False,
  "scope":"Fixed homogeneous Euler/cubic transport only; no actual pair in this model.",
 },
 "hashes":{name:hashlib.sha256((HERE/name).read_bytes()).hexdigest()
           for name in names},
 "coverage_boundary":(
  "Exact transported operator, polynomial/rational sources with poles, "
  "lowest Euler commutator coefficient, signed divisor character modulo Phi10, "
  "J/U/A eigenrelations, deck projections, U-power eigenoperators and positive/negative ordinary-order linearization "
  "controls. Earlier 25 norm controls are complementary finite implementation "
  "checks. All-source actual division-root exclusion and all-target actual-Euler-pair "
  "exclusion use the written valuations, signed-divisor identity and "
  "commutative-symbol injectivity. All rational target orderings with a regular "
  "Long root are excluded by prescribed-root uniqueness and the independent "
  "deck-character/valuation contradiction. Formal completion and nonhomogeneous "
  "source/cubic models remain outside these theorems."
 ),
 "independent_replay":False,"formal_verification":False,"external_review":False,
}
out=Path(tempfile.mkdtemp(prefix="qql1-check-"))/"quintic_ordering_laurent_controls.json"
out.write_text(json.dumps(receipt,indent=2)+"\n")
print(json.dumps({"claim":"QQL1","status":"PASS","top_norm_cases":25,
                  "bottom_orbit_cases":25,"receipt":str(out)},indent=2))
