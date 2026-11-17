#!/usr/bin/env python3
"""Exact generic repeated-profile cubic sphere calculation.

p=y*(t-r)^2+B(t)+x*A(t)+x^2*D(t)+E*x^3, with deg(B,A,D)=(2,2,1).
The t^3 coefficient has been removed by an orthogonal null rotation.
For transcendental s=r^2 and r!=0, scale B1=r*b1,A1=r*a1,D1=r*d1.
The first twelve moments force B=A=D=E=0 over QQ(s). This does NOT
certify every algebraic specialization r: an additional specialization
locus, beyond the explicitly handled triangular pivots, remains possible.

The proof uses the positive-weight specialization lemma twice: first
Gauss valuation at p=32003, then s->1. The checker verifies integrality
at both steps, and the exact finite-field homogeneous ideal has only the
origin by explicit pure-power normal forms. See the canonical note for
this mathematical lifting argument. Requires SymPy 1.14.0, Singular 4.3.2.
"""
from functools import lru_cache
from math import factorial, comb
from pathlib import Path
import hashlib
import json
import subprocess
import sympy as S
from sympy.polys.rings import ring
from sympy.polys.domains import QQ

PRIME=32003
SPECIALIZATION=1
WEIGHTS=(1,1,1,2,2,2,3,3,4)
PHASE=(0,0,0,1,1,1,2,2,3)
TPOWER=(0,1,2,0,1,2,0,1,0)
ODD=(0,1,0,0,1,0,0,1,0)


def generate():
    K=QQ.frac_field('s')
    ss=K.gens[0]
    R,*variables=ring('B0,b1,B2,A0,a1,A2,D0,d1,E',K)
    B0,b1,B2,A0,a1,A2,D0,d1,E=variables
    @lru_cache(None)
    def df(n):
        return 1 if n<=0 else n*df(n-2)
    def exponents(n,index=0,prefix=()):
        if index==len(WEIGHTS):
            if n==0:
                yield prefix
            return
        for count in range(n//WEIGHTS[index]+1):
            yield from exponents(n-count*WEIGHTS[index],index+1,prefix+(count,))
    @lru_cache(None)
    def integral(y_count,tp,od):
        value=K.zero
        # Expand (t-r)^(2*y_count), with B1=r*b1 etc supplying r^od.
        for ell in range(2*y_count+1):
            if (tp+ell)%2:
                continue
            r_power=2*y_count-ell+od
            assert r_power%2==0
            value+=K(QQ((-1)**ell*comb(2*y_count,ell)*df(tp+ell-1),df(2*y_count+tp+ell+1)))*ss**(r_power//2)
        return value
    def moment(n):
        result={}
        for powers in exponents(n):
            y_count=sum(a*b for a,b in zip(powers,PHASE))
            tp=sum(a*b for a,b in zip(powers,TPOWER))
            od=sum(a*b for a,b in zip(powers,ODD))
            denominator=1
            for count in powers:
                denominator*=factorial(count)
            coefficient=K(QQ(factorial(n),denominator))*integral(y_count,tp,od)
            if coefficient:
                result[powers]=coefficient
        return R.from_dict(result)
    eliminate=(B0,A0,D0,E)
    pivots=(K.one,2*(5*ss+1)/15,2*(21*ss**2+18*ss+1)/105,8*(429*ss**3+715*ss**2+195*ss+5)/15015)
    solutions=[]
    moments=[]
    for n in range(1,13):
        equation=moment(n)
        for variable,value in zip(eliminate,solutions):
            equation=equation.compose(variable,value)
        assert all(sum(a*b for a,b in zip(powers,WEIGHTS))==n for powers in equation)
        if n<=4:
            variable=eliminate[n-1]
            constant=equation.compose(variable,R.zero)
            solution=-constant/pivots[n-1]
            assert equation==pivots[n-1]*(variable-solution)
            solutions.append(solution)
        else:
            moments.append(equation)
    return K,R,variables,pivots,solutions,moments,moment


def coefficient_mod(value,prime=PRIME,sval=SPECIALIZATION):
    def eval_poly(poly):
        result=0
        for powers,coefficient in poly.items():
            numerator=int(coefficient.numerator)
            denominator=int(coefficient.denominator)
            assert denominator%prime!=0
            result=(result+numerator*pow(denominator,-1,prime)*pow(sval,powers[0],prime))%prime
        return result
    # This also establishes Gauss p-integrality: every rational scalar
    # denominator is a p-unit, and denominator mod p is not identically zero.
    assert any(int(c.numerator)%prime for c in value.denom.values())
    numerator=eval_poly(value.numer)
    denominator=eval_poly(value.denom)
    assert denominator!=0
    return numerator*pow(denominator,-1,prime)%prime


def singular_polynomial(poly):
    terms=[]
    for powers,coefficient in sorted(poly.items(),reverse=True):
        value=coefficient_mod(coefficient)
        if not value:
            continue
        factors=[str(value)]
        for name,power in zip(('B0','b1','B2','A0','a1','A2','D0','d1','E'),powers):
            if power:
                factors.append(name if power==1 else name+'^'+str(power))
        terms.append('*'.join(factors))
    return '+'.join(terms) if terms else '0'


def verify_parameter_minor(moments):
    """Define an exact nonzero D(s) by a retained weight-18 matrix minor."""
    from math import lcm,gcd
    from flint import nmod_mat,__version__ as flint_version
    assert flint_version=='0.9.0'
    weights=(1,1,2,2,3)
    @lru_cache(None)
    def monomials(weight,index=0):
        if index==5:
            return ((),) if weight==0 else ()
        result=[]
        for count in range(weight//weights[index]+1):
            result.extend((count,)+tail for tail in monomials(weight-count*weights[index],index+1))
        return tuple(result)
    inputs=[]
    denominators=[]
    for polynomial in moments:
        denominator,numerator=polynomial.clear_denoms()
        denominators.append(str(denominator))
        # No finite poles beyond the three explicitly treated pivots.
        parameter=denominator.ring.gens[0]
        factors=(5*parameter+1,21*parameter**2+18*parameter+1,429*parameter**3+715*parameter**2+195*parameter+5)
        rest=denominator
        for factor in factors:
            while rest.degree()>=factor.degree():
                quotient,remainder=rest.div(factor)
                if remainder:
                    break
                rest=quotient
        assert rest.degree()==0 and rest!=0
        assert all(coefficient.denom.degree()==0 for coefficient in numerator.values())
        integer_denominator=lcm(*(int((value/next(iter(coefficient.denom.values()))).denominator) for coefficient in numerator.values() for value in coefficient.numer.values()))
        entries=[]
        values=[]
        for powers,coefficient in sorted(numerator.items()):
            assert all(powers[i]==0 for i in (0,3,6,8))
            projected=tuple(powers[i] for i in (1,2,4,5,7))
            scalar_denominator=next(iter(coefficient.denom.values()))
            coefficients={exponent[0]:int(value/scalar_denominator*integer_denominator) for exponent,value in coefficient.numer.items()}
            values.extend(abs(value) for value in coefficients.values())
            entries.append([list(projected),coefficients])
        content=gcd(*values)
        for powers,coefficients in entries:
            for exponent in coefficients:
                coefficients[exponent]//=content
        inputs.append(entries)
    canonical=json.dumps(inputs,separators=(',',':'),sort_keys=True)
    path=Path(__file__).resolve().parents[1]/'extended-geometry'/'two_pair_image_cubic_repeated_minor_selector.json'
    selector=json.loads(path.read_text())
    assert selector['format']=='two-pair-repeated-profile-macaulay-minor-v1'
    assert selector['target_weight']==18
    assert selector['variable_order']==['b1','B2','a1','A2','d1']
    assert selector['weights']==list(weights)
    assert selector['moment_orders']==list(range(5,13))
    assert selector['integer_input_sha256']==hashlib.sha256(canonical.encode()).hexdigest()
    assert selector['clearing_denominators']==denominators
    rows=monomials(18)
    row_indices={powers:i for i,powers in enumerate(rows)}
    labels=[]
    for order in range(5,13):
        labels.extend([order,list(powers)] for powers in monomials(18-order))
    assert selector['matrix_shape']==[len(rows),len(labels)]==[874,1187]
    selected=selector['selected_column_indices']
    assert len(selected)==len(set(selected))==len(rows)
    assert all(isinstance(j,int) and 0<=j<len(labels) for j in selected)
    assert selector['selected_column_labels']==[labels[j] for j in selected]
    specialization=selector['specialization']
    assert specialization['prime']==PRIME and specialization['s']==1
    matrix=[[0]*len(rows) for unused in rows]
    row_degrees=[-1]*len(rows)
    for column,index in enumerate(selected):
        order,multiplier=labels[index]
        for powers,coefficients in inputs[order-5]:
            exponent=tuple(a+b for a,b in zip(multiplier,powers))
            row=row_indices[exponent]
            matrix[row][column]=sum(coefficients.values())%PRIME
            row_degrees[row]=max(row_degrees[row],max(coefficients))
    determinant=int(nmod_mat(matrix,PRIME).det())
    assert determinant==specialization['determinant']==42
    assert min(row_degrees)>=0
    degree_bound=sum(row_degrees)
    assert selector['determinant_degree_upper_bound']==degree_bound==19448
    pure_powers={'b1':18,'B2':18,'a1':9,'A2':9,'d1':6}
    assert selector['coordinate_pure_powers']==pure_powers
    for index,power in enumerate(pure_powers.values()):
        exponent=tuple(power if j==index else 0 for j in range(5))
        assert exponent in row_indices
    # A corrupted matrix with a zero selected column loses the certificate.
    for row in matrix:
        row[0]=0
    assert nmod_mat(matrix,PRIME).det()==0
    return {
      'selector':path.name,
      'selector_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),
      'flint_version':flint_version,
      'target_weight':18,
      'matrix_shape':[874,1187],
      'selected_square_size':874,
      'integer_input_sha256':selector['integer_input_sha256'],
      'definition':'D(s) is the exact determinant of the retained 874 selected columns of the integer-polynomial weight-18 moment coefficient matrix.',
      'nonvanishing_witness':'D(1) mod 32003 = 42',
      'degree_upper_bound':19448,
      'pure_powers':pure_powers,
      'covered_parameters':'s*(5s+1)*(21s^2+18s+1)*(429s^3+715s^2+195s+5)*D(s) != 0',
      'remaining_specializations':'D(s)=0 away from s=0 and the three fixed pivot polynomials; these roots are not classified by this determinant certificate.',
      'negative_control':'Zeroing one selected matrix column makes its determinant vanish.',
    }

def main():
    assert S.__version__=='1.14.0'
    version=subprocess.run(['Singular','--version'],capture_output=True,text=True,check=True).stdout.splitlines()[0]
    assert 'version 4.3.2 ' in version
    K,R,variables,pivots,solutions,moments,moment=generate()
    # Independent direct sphere expansion at s=r=1 for the first moments.
    x,y,t=S.symbols('x y t')
    B0,b1,B2,A0,a1,A2,D0,d1,E=S.symbols('B0 b1 B2 A0 a1 A2 D0 d1 E')
    p=y*(t-1)**2+B0+b1*t+B2*t*t+x*(A0+a1*t+A2*t*t)+x*x*(D0+d1*t)+E*x**3
    for n in range(1,4):
        direct=0
        for (i,j,k),coefficient in S.Poly(S.expand(p**n),x,y,t).terms():
            if i==j and k%2==0:
                direct+=coefficient*S.factorial(i)*(S.factorial2(k-1) if k else 1)/S.factorial2(2*i+k+1)
        assert S.expand(moment(n).as_expr().subs(S.Symbol('s'),1)-direct)==0
    # Every coefficient, including the elimination expressions, is regular
    # for the two specializations used by the positive-weight lemma.
    for polynomial in [R(p) for p in pivots]+solutions+moments:
        for coefficient in polynomial.values():
            coefficient_mod(coefficient)
    polynomials=[singular_polynomial(p) for p in moments]
    script='ring r='+str(PRIME)+',(b1,B2,a1,A2,d1),wp(1,1,2,2,3);\n'
    script+='ideal I='+','.join(polynomials)+';\nideal G=std(I);\n'
    script+='print("DIM "+string(dim(G)));print("VDIM "+string(vdim(G)));\n'
    for name,exponent in [('b1',17),('B2',17),('a1',9),('A2',9),('d1',6)]:
        script+='print("POWER '+name+' "+string(reduce('+name+'^'+str(exponent)+',G)));\n'
    script+='ideal Bad=0;print("NEGATIVE "+string(reduce(b1^17,Bad)!=0));\nquit;\n'
    result=subprocess.run(['Singular','-q'],input=script,capture_output=True,text=True,check=True)
    lines=[line.strip() for line in result.stdout.splitlines() if line.strip()]
    expected=['DIM 0','VDIM 797','POWER b1 0','POWER B2 0','POWER a1 0','POWER A2 0','POWER d1 0','NEGATIVE 1']
    assert lines==expected,(lines,result.stderr)
    parameter_minor=verify_parameter_minor(moments)
    record={
      'format':'two-pair-image-cubic-repeated-generic-v1',
      'checker_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
      'command':'.venv/bin/python research/foundations/scripts/verify_two_pair_image_cubic_repeated_generic.py',
      'software':{'sympy':S.__version__,'singular':version},
      'family':'p=y*(t-r)^2+B(t)+x*A(t)+x^2*D(t)+E*x^3; deg(B,A,D)=(2,2,1)',
      'field':'QQ(s), s=r^2 transcendental; extension to algebraic closure permitted',
      'coefficient_scaling':'B1=r*b1,A1=r*a1,D1=r*d1; hence r!=0',
      'moment_cutoff':12,
      'exact_parameter_minor':parameter_minor,
      'weights':{'b1':1,'B2':1,'a1':2,'A2':2,'d1':3},
      'triangular_pivots':[str(p) for p in pivots],
      'triangular_solutions':[str(p.as_expr()) for p in solutions],
      'residual_moment_term_counts':[len(p) for p in moments],
      'prime':PRIME,'s_specialization':SPECIALIZATION,
      'finite_field_input':polynomials,
      'finite_field_result':lines,
      'pure_power_exponents':{'b1':17,'B2':17,'a1':9,'A2':9,'d1':6},
      'lifting_argument':'Weighted projective properness for the Gauss p-valuation followed by the s-1 valuation; both coefficient regularity conditions checked.',
      'claim':'Outside s=0, the three pivot polynomials and the explicitly determinant-defined finite set D(s)=0, moments 1..12 force G=0 after the uniform null rotation removing t^3.',
      'specialization_limit':'The remaining finite locus is explicitly D(s)=0 away from the fixed pivot loci; its non-pivot roots have not been classified. D is retained by its exact selected determinant definition, without expansion.',
      'full_repeated_profile_coverage':False,
      'full_sixteen_parameter_coverage':False,
      'assurance':{'internal_exact_replay':True,'independent_replay':False,'formal_verification':False,'external_review':False},
    }
    output=Path(__file__).resolve().parents[1]/'extended-geometry'/'two_pair_image_cubic_repeated_generic_verification.json'
    output.write_text(json.dumps(record,indent=2)+'\n')
    print('PASS generic repeated-profile sphere cubic over QQ(r): first twelve moments force G=0')
    print('PASS coefficient regularity, weighted homogeneity, exact finite-field basis and five pure-power normal forms')
    print('PASS exact weight-18 parameter determinant: 874 selected columns, D(1)=42 mod32003, degree at most19448')
    print('The non-pivot roots of the explicit determinant D(s) remain unclassified.')


if __name__=='__main__':
    main()
