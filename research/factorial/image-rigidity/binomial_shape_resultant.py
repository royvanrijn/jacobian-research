"""Exact replay for BINOMIAL_SHAPE_RESULTANT.tex.

The uniform normalization and repeated-part reductions have written proofs.
Finite symbolic controls below do not extrapolate coefficient positivity
or squarefreeness. The one beta-Hurwitz certificate is verified completely
by a degree-bounded integer resultant identity and two exact matrix methods.
"""
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import product
from math import comb, factorial, prod
from pathlib import Path
import json
import platform
import sys
import sympy as sp
from flint import fmpz_mat, fmpz_poly

sys.set_int_max_str_digits(0)
HERE=Path(__file__).resolve().parent
TEX="BINOMIAL_SHAPE_RESULTANT.tex"
CERT="binomial_shape_resultant_beta_hurwitz_certificate.json"
HKK="sources/holtz-khrushchev-kushel-2015-v1.pdf"
HKK_HASH="05aaeaeddc525d095461d161bab6417415c36e22537f9322e5760a69c7289db5"
beta,t,u=sp.symbols("beta t u")

def rising(z,k):
    value=1
    for j in range(k):
        value*=z+j
    return value

@lru_cache(None)
def P(n,a,d):
    return sum(comb(n,k)*rising(beta+a*n,d*k)*t**k
               for k in range(n+1))

def symbolic_controls():
    cases=extra=0
    for n,a,d in product(range(1,5),range(3),range(1,4)):
        r=sp.Poly(sp.resultant(P(n,a,d),P(n+1,a,d),t),beta,domain=sp.ZZ)
        deg=(2*d-1)*n*(n+1)//2
        lc=d**(n*(n+1))*prod(j**j for j in range(1,n+1))
        dd=sp.Poly(prod(rising(beta+a*(n+1),max(d*j-a,0))
                        for j in range(1,n+1)),beta,domain=sp.ZZ)
        q,rem=r.div(dd)
        assert rem.is_zero
        assert r.degree()==deg and r.LC()==lc
        assert q.degree()==deg-sum(max(d*j-a,0) for j in range(1,n+1))
        assert q.LC()>0
        if a==1:
            e=sp.Poly((beta+(d+1)*n+d)**n,beta)
            assert sp.gcd(dd,e).degree()==0
            assert q.rem(e).is_zero
            extra+=1
        cases+=1
    for n,d in product(range(1,7),range(1,5)):
        value=-(d+1)*n-d
        assert sp.expand(P(n+1,1,d).subs(beta,value)
                         -P(n,1,d).subs(beta,value))==0
        extra+=1
    return {"symbolic_resultant_normalization_cases":cases,
            "lower_one_extra_factor_controls":extra}

def gaussian_controls():
    cases=0
    for a,m,d,k in product([0,2],[1,4],range(1,5),range(9)):
        numerator=fmpz_poly([0])
        for j in range(k+1):
            term=fmpz_poly([1])
            for i in range(d*j):
                term*=fmpz_poly([1,a*m+i])
            numerator+=(-1)**(k-j)*comb(k,j)*term
        for h in range((k+1)//2):
            assert numerator[h]==0
        if k%2==0:
            gaussian=factorial(k)//(2**(k//2)*factorial(k//2))
            assert numerator[k//2]==d**k*gaussian
        cases+=1
    return {"finite_formal_gaussian_moment_controls":cases}

def reciprocal_at(n,a,d,shape):
    return fmpz_poly([comb(n,n-j)*rising(shape+a*n,d*(n-j))
                     for j in range(n+1)])

def evaluate_integer_poly(desc,z):
    value=0
    for c in desc:
        value=value*z+c
    return value

def verify_counter_polynomial(certificate,coefficients):
    assert certificate["parameters"]=={"n":7,"a":0,"d":5}
    assert len(coefficients)==113 and coefficients[0]!=0
    content=int(certificate["positive_content"])
    assert content>0
    for shape in range(1,254):
        f=reciprocal_at(7,0,5,shape)
        g=reciprocal_at(8,0,5,shape)
        actual=int(f.resultant(g))
        dd=prod(rising(shape,5*j) for j in range(1,8))
        proposed=content*dd*evaluate_integer_poly(coefficients,shape)
        assert actual==proposed,shape

def counter_certificate():
    cert=json.loads((HERE/CERT).read_text())
    co=[int(x) for x in cert["primitive_coefficients_descending"]]
    verify_counter_polynomial(cert,co)
    assert all(x>0 for x in co)
    assert cert["resultant_degree"]==252
    assert cert["divisor_degree"]==140
    assert cert["primitive_degree"]==112
    content=int(cert["positive_content"])
    assert content*co[0]==5**56*prod(j**j for j in range(1,8))
    primitive_content=0
    from math import gcd
    for c in co:
        primitive_content=gcd(primitive_content,c)
    assert primitive_content==1
    corrupted=co.copy()
    corrupted[-1]+=1
    try:
        verify_counter_polynomial(cert,corrupted)
    except AssertionError:
        pass
    else:
        raise AssertionError("Corrupted coefficient certificate was accepted")
    h=cert["negative_minor_order"]
    assert h==84
    def cc(j):
        return co[j] if 0<=j<len(co) else 0
    H=fmpz_mat([[cc(2*j-i) for j in range(1,h+1)]
                              for i in range(1,h+1)])
    determinant=int(H.det())
    assert determinant<0
    assert len(str(abs(determinant)))==cert["negative_minor_decimal_digits"]==7555
    det_hash=sha256(str(determinant).encode()).hexdigest()
    assert det_hash==cert["negative_minor_sha256"]
    # Independent rational row elimination, rather than the determinant algorithm.
    fco=[Fraction(c) for c in co]
    width=(len(fco)+1)//2
    def pad(row):
        return row+[Fraction(0)]*(width-len(row))
    first,second=pad(fco[::2]),pad(fco[1::2])
    first_column=[first[0],second[0]]
    for row in range(2,85):
        assert second[0]!=0
        nxt=[(second[0]*first[j+1]-first[0]*second[j+1])/second[0]
             for j in range(width-1)]+[Fraction(0)]
        first_column.append(nxt[0])
        first,second=second,nxt
    assert all(c>0 for c in first_column[:84])
    assert first_column[84]<0
    assert sha256((HERE/HKK).read_bytes()).hexdigest()==HKK_HASH
    return {"complete_beta_hurwitz_certificate":{
        "parameters":cert["parameters"],
        "primitive_degree":112,
        "integer_resultant_identity_nodes":253,
        "identity_degree_bound":252,
        "all_coefficients_strictly_positive":True,
        "negative_hurwitz_minor_order":84,
        "negative_hurwitz_minor_decimal_digits":7555,
        "negative_hurwitz_minor_sha256":det_hash,
        "rational_routh_first_negative_row":84,
        "corrupted_coefficient_rejection":True,
        "right_half_plane_shape_zero":True,
        "positive_real_shape_zero_asserted":False,
        "SFC_counterexample_asserted":False},
        "primary_hurwitz_source":{"file":HKK,"sha256":HKK_HASH,
             "scope":"Theorem1, p1, classical real-polynomial Routh-Hurwitz criterion"}}

def replay():
    controls={}
    controls.update(symbolic_controls())
    controls.update(gaussian_controls())
    controls.update(counter_certificate())
    return {"kind":"internal exact replay; universal proofs are the TeX source",
      "scope":{
        "natural_divisor":"all integer a>=0,d>=1,n>=1",
        "exact_shape_resultant_degree":"(2d-1)n(n+1)/2",
        "exact_positive_leading_coefficient":"d^(n(n+1))*prod_(j1..n)j^j",
        "additional_lower_one_divisor":"(beta+(d+1)n+d)^n, all n,d>=1",
        "nonnegative_positive_shape_resultant":"all beta>0; every positive shape zero is repeated",
        "residual_coefficient_positivity_universally_proved":False,
        "residual_squarefreeness_universally_proved":False,
        "beta_hurwitz_shortcut":"false at a0,d5,n7 by complete exact certificate",
        "all_ordinary_binomial_sfc":False},
      "controls":controls,
      "python":platform.python_version(),"sympy":sp.__version__,
      "source_sha256":{name:sha256((HERE/name).read_bytes()).hexdigest()
          for name in [TEX,Path(__file__).name,CERT,HKK]},
      "independent_replay":False,"formal_verification":False,
      "external_review":False,"passed":True}

if __name__=="__main__":
    report=replay()
    report["command"]="python research/factorial/image-rigidity/"+Path(__file__).name
    (HERE/"binomial_shape_resultant_verification.json").write_text(
        json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
