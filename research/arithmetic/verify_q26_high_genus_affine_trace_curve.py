#!/usr/bin/env python3
"""Exact Q26 accessory binding, S12 witnesses and modular genus certificate.

The modular discriminant is checked both by bivariate FLINT and by 3169
independent ordinary univariate Euclidean resultants in F_(29^3). The
latter is a complete polynomial identity test using a proven degree bound.
"""
from __future__ import annotations
import argparse
from functools import reduce
import hashlib
import json
import math
from pathlib import Path
import platform
import flint
import sympy as sp

HERE=Path(__file__).resolve().parent
ROOT=HERE.parents[1]
t=sp.Symbol('t')

def require(test,message):
    if not test:raise AssertionError(message)

def digest(path):return hashlib.sha256(path.read_bytes()).hexdigest()

def qqpoly(expression):
    f=sp.Poly(sp.sympify(expression),t,domain=sp.QQ)
    return flint.fmpq_poly([flint.fmpq(int(f.nth(j).p),int(f.nth(j).q))
                           for j in range(f.degree()+1)])

def primitive_curve():
    src=json.loads((HERE/'odd_affine_critical_trace_curve_inputs.json').read_text())
    rows=src['curves']['6']['accessory_s_coefficients']
    ns=[qqpoly(rows[str(j)]['numerator']) for j in range(13)]
    ds=[qqpoly(rows[str(j)]['denominator']) for j in range(13)]
    common=flint.fmpq_poly([1])
    for d in ds:
        common=(common*d)//common.gcd(d)
        common=common/common[common.degree()]
    coeff=[]
    for n,d in zip(ns,ds):
        q=common//d
        require(q*d==common,'Common denominator exact division')
        coeff.append(q*n)
    content=reduce(lambda a,b:a.gcd(b),coeff)
    require(content==flint.fmpq_poly([1]),'Raw LCD-cleared vertical component')
    require(common.degree()==144,'Raw least-common-denominator degree')
    coeff=[a//content for a in coeff]
    # This characteristic-zero check is essential; a modular gcd is insufficient.
    require(reduce(lambda a,b:a.gcd(b),coeff)==flint.fmpq_poly([1]),
            'Nonconstant Q[t] coefficient content remains')
    den=math.lcm(*(int(a.denom()) for a in coeff))
    coeff=[a*den for a in coeff]
    integers=[[int(c) for c in a] for a in coeff]
    ic=math.gcd(*(c for a in integers for c in a))
    integers=[[c//ic for c in a] for a in integers]
    if integers[-1][-1]<0:integers=[[-c for c in a] for a in integers]
    require([len(a)-1 for a in integers]==list(range(144,131,-1)),
            'Individual primitive coefficient degrees')
    require(math.gcd(*(c for a in integers for c in a))==1,'Integer content')
    require(reduce(lambda a,b:a.gcd(b),[flint.fmpq_poly(a) for a in integers])==
            flint.fmpq_poly([1]),'Primitive coefficient gcd')
    return {'degree_s':12,'degree_t':max(len(a)-1 for a in integers),
            'coefficients_ascending_s':integers}

def evaluate_int(coeff,value):
    out=0
    for a in reversed(coeff):out=out*value+a
    return out

def prime_divisors(n):
    ans=[];p=2
    while p*p<=n:
        if n%p==0:
            ans.append(p)
            while n%p==0:n//=p
        p+=1
    if n>1:ans.append(n)
    return ans

def frobenius_certificate(curve,data):
    raw=[evaluate_int(a,2) for a in curve['coefficients_ascending_s']]
    ic=math.gcd(*raw)
    primitive=[c//ic for c in raw]
    if primitive[-1]<0:primitive=[-c for c in primitive]
    require(primitive==data['primitive_coefficients_ascending'],'t2 primitive binding')
    require(len(primitive)==13 and primitive[-1]!=0,'t2 degree12')
    results={}
    expected={'irreducible12':[12],'cycle11_plus1':[1,11],
              'single_transposition_power':[1,2,9]}
    for name,row in data['witnesses'].items():
        p=row['prime']; f=flint.nmod_poly(primitive,p)
        require(f.degree()==12,'Degree drop at Frobenius prime')
        require(f.gcd(f.derivative()).degree()==0,'Repeated modular roots')
        product=flint.nmod_poly([row['scalar']],p)
        degrees=[]
        for cs in row['factors_ascending']:
            g=flint.nmod_poly(cs,p); x=flint.nmod_poly([0,1],p)
            degree=g.degree();degrees.append(degree)
            require(g[degree]==1,'Nonmonic factor')
            require(pow(x,p**degree,g)==x%g,'Rabin terminal power')
            for ell in prime_divisors(degree):
                require(g.gcd(pow(x,p**(degree//ell),g)-x).degree()==0,
                        'Rabin proper subfield gcd')
            product*=g
        require(product==f,'Exact factor product')
        require(sorted(degrees)==expected[name],'Unexpected Frobenius type')
        results[name]={'prime':p,'factor_degrees':degrees,
                       'leading_coefficient_unit':True,'separable':True,
                       'factor_products_and_Rabin_irreducibility':True}
    # a^11 f(X/a) is monic integral and preserves factor types at these units.
    a=primitive[-1]
    monic=[primitive[j]*a**(11-j) for j in range(12)]+[1]
    require(all(isinstance(c,int) for c in monic),'Monic integral transform')
    require(all(a%row['prime'] for row in data['witnesses'].values()),
            'Monic scaling leading coefficient not a unit')
    return {'specialization_t':2,'primitive_coefficients_ascending':primitive,
            'monic_integral_scaling_leading':a,'witnesses':results,
            'arithmetic_generic_group':'S12',
            'transposition_power_for_type_1_2_9':9}

def modular_certificate(curve,data):
    p=data['prime'];cs=curve['coefficients_ascending_s']
    ctx=flint.nmod_mpoly_ctx.get(('s','t'),p,'lex')
    F=ctx.from_dict({(i,j):int(c)%p for i,a in enumerate(cs)
                     for j,c in enumerate(a) if int(c)%p})
    raw=F.discriminant('s').to_dict()
    require(all(i==0 for i,j in raw),'Resultant still depends on s')
    ds=[0]*(max(j for i,j in raw)+1)
    for (i,j),c in raw.items():ds[j]=int(c)
    require(ds==data['disc_coefficients'],'Full bivariate discriminant binding')
    D=flint.nmod_poly(ds,p);G=D.gcd(D.derivative())
    A=D//G; A=A//A.gcd(G)
    LC=flint.nmod_poly(cs[-1],p)
    good=A//A.gcd(LC)
    require(good.degree()==2124,'Good simple-root count')
    require(good.gcd(good.derivative()).degree()==0,'Simple-root polynomial')
    require(good.gcd(LC).degree()==0,'Simple root at leading zero')
    require([int(c) for c in A]==data['multiplicity_one_coefficients'],
            'Multiplicity-one input mismatch')
    require([int(c) for c in good]==data['good_simple_coefficients'],
            'Good simple-root input mismatch')
    require([int(c) for c in LC]==data['LC_coefficients'],'Leading input mismatch')
    return {'prime':p,'discriminant_degree':D.degree(),
            'gcd_with_derivative_degree':G.degree(),
            'multiplicity_one_degree':A.degree(),'good_simple_degree':good.degree(),
            'leading_coefficient_degree':LC.degree(),
            'simple_root_extractor':'(D/G)/gcd(D/G,G), G=gcd(D,Dprime)',
            'simple_roots_coprime_to_leading':True}

def interpolation_certificate(curve,data,modular):
    p=data['prime'];required=data['required_good_evaluations']
    require(data['degree_bound']==22*curve['degree_t']==3168,'Degree bound')
    require(required==data['degree_bound']+1,'Insufficient exact evaluations')
    modulus=flint.nmod_poly(data['modulus_ascending'],p)
    x=flint.nmod_poly([0,1],p)
    require(modulus.degree()==3 and
            pow(x,p**3,modulus)==x%modulus and
            modulus.gcd(pow(x,p,modulus)-x).degree()==0,'Extension modulus')
    ctx=flint.fq_default_ctx(modulus=flint.fmpz_mod_poly_ctx(p)(data['modulus_ascending']),
                             var='z',fq_type='FQ_ZECH')
    P=flint.fq_default_poly_ctx(ctx);zero=ctx(0);one=ctx(1)
    coeff=[[ctx(c) for c in a] for a in curve['coefficients_ascending_s']]
    ds=[ctx(c) for c in modular['disc_coefficients']]
    require(len(ds)-1<=data['degree_bound'],'Retained discriminant degree')
    def evaluate(a,value):
        out=zero
        for c in reversed(a):out=out*value+c
        return out
    def resultant(a,b):
        result=one
        while b.degree()>0:
            m,n=a.degree(),b.degree();rem=a%b
            if rem.is_zero():return zero
            result=result*((-1)**(m*n))*b.leading_coefficient()**(m-rem.degree())
            a,b=b,rem
        return result*b[0]**a.degree()
    checked=0;skipped=0
    for index in range(p**3):
        value=ctx([index%p,(index//p)%p,index//(p*p)])
        specialized=[evaluate(a,value) for a in coeff]
        if specialized[-1]==zero:
            skipped+=1;continue
        f=P(specialized)
        # For degree12 the discriminant sign (-1)^66 is positive.
        actual=resultant(f,f.derivative())/specialized[-1]
        require(actual==evaluate(ds,value),'Independent discriminant evaluation')
        checked+=1
        if checked==required:break
    require(checked==required,'Not enough distinct good field values')
    return {'prime':p,'extension_degree':3,'modulus_ascending':data['modulus_ascending'],
            'degree_bound':data['degree_bound'],'distinct_good_evaluations':checked,
            'leading_zero_skips':skipped,'last_base_p_index':index,
            'method':'Ordinary univariate Euclidean resultant, independent of bivariate discriminant',
            'complete_polynomial_identity_by_degree_bound':True}

def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--write-receipt',type=Path,
                        default=HERE/'q26_high_genus_affine_trace_curve_verification.json')
    args=parser.parse_args()
    data=json.loads((HERE/'q26_high_genus_affine_trace_curve_inputs.json').read_text())
    curve=primitive_curve()
    require(curve==data['curve'],'Exact primitive accessory curve mismatch')
    require(curve['degree_t']==144,'Unexpected accessory t-degree')
    frobenius=frobenius_certificate(curve,data['specialization'])
    modular=modular_certificate(curve,data['modular_discriminant'])
    independent=interpolation_certificate(curve,data['independent_interpolation'],
                                          data['modular_discriminant'])
    bound=(modular['good_simple_degree']-22)//2
    require(bound==data['genus_lower_bound']==1051,'Genus arithmetic')
    paths=[Path(__file__).resolve(),HERE/'Q26_HIGH_GENUS_AFFINE_TRACE_CURVE.md',
           HERE/'q26_high_genus_affine_trace_curve_inputs.json',
           HERE/'odd_affine_critical_trace_curve_inputs.json',
           HERE/'ODD_AFFINE_CRITICAL_TRACE_CURVES.md',
           HERE/'sources/milne-ant300.pdf',HERE/'sources/faltings-1983-endlichkeitssaetze.pdf']
    result={'schema_version':1,'result':'PASS','python':platform.python_version(),
            'sympy':sp.__version__,'python_flint':flint.__version__,
            'source_sha256':{str(p.relative_to(ROOT)):digest(p) for p in paths},
            'curve':{'degree_s':12,'degree_t':144,'Q_t_coefficient_content':1,
                     'raw_LCD_degree':144,'raw_LCD_cleared_coefficient_content':1,
                     'integer_content':1,'primitive_input_exactly_rebuilt':True},
            'arithmetic_frobenius':frobenius,'modular_discriminant':modular,
            'independent_interpolation':independent,
            'geometric_splitting_group':'S12','genus_lower_bound':bound,
            'imported_finiteness':'Faltings Satz7, printed page365',
            'scope':'Finite good rational points in the normalized odd-numerator affine-trace slice; no effective point list, rational existence/nonexistence determination, exact genus, or all-pencil classification.',
            'formal_verification':False,'independent_external_replay':False,
            'external_human_review':False}
    args.write_receipt.write_text(json.dumps(result,indent=2)+'\n')
    print('PASS: exact Q26 accessory binding, S12, 2124 simple branches and genus >=1051.')

if __name__=='__main__':main()
