"""Finite actual-period jet matching for a holonomic D1 cover.

This checks selected periods of a SUPPLIED presentation. Its application to
ALL Gaussian multipliers needs a proved full-amplitude cover map. The map's
validity/coverage is NOT inferred from finite moment agreement.
"""
from __future__ import annotations
from functools import lru_cache
from math import factorial
from pathlib import Path
import importlib.util
import sympy as sp

_spec=importlib.util.spec_from_file_location("gaussian_hom_side",Path(__file__).with_name("gaussian_polynomial_hom_sidecheck_20261116.py"))
H=importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(H)


def origin_certificate(operator):
    if not operator:
        raise ValueError("Nonzero annihilator required")
    shift=min(a-b for a,b in operator)
    p=sp.Poly(sp.expand(sum(c*sp.prod(H.Z-j for j in range(b))
                           for (a,b),c in operator.items() if a-b==shift)),H.Z,domain=sp.QQ)
    assert not p.is_zero
    roots=sorted(int(x) for x in sp.polys.polytools.ground_roots(p) if x.is_Integer and x>=0)
    return dict(shift=shift,origin_indicial_polynomial=str(p.as_expr()),
                nonnegative_integer_roots=roots,jet_degree=max(roots,default=-1))


def match_actual_jets(data,jets,hom_result=None):
    result=hom_result if hom_result is not None else H.solve_presentation(data)
    if result['outcome']!='COMPLETE':
        return dict(outcome='UNKNOWN',reason='Polynomial Hom computation incomplete')
    certs=data.get('coordinate_annihilators',result.get('generated_coordinate_annihilators'))
    origin=[origin_certificate(H.parse_op(c['operator'])) for c in certs]
    if len(jets)!=data['generator_count']:
        raise ValueError('One actual jet is needed per generator')
    for jet,cert in zip(jets,origin):
        if len(jet)<cert['jet_degree']+1:
            raise ValueError('Insufficient actual jets for the certified origin bound')
    keys=[(i,n) for i,cert in enumerate(origin) for n in range(cert['jet_degree']+1)]
    basis=result['polynomial_hom_basis']
    columns=len(basis)
    rows=[]
    actual=[]
    for i,n in keys:
        rows.append([sp.Rational(b[i]['coefficients'][n]) if n<len(b[i]['coefficients']) else sp.S.Zero for b in basis])
        actual.append(sp.Rational(jets[i][n]))
    matrix=sp.Matrix(rows) if rows and columns else sp.zeros(len(keys),columns)
    target=sp.Matrix(len(keys),1,actual)
    # These jets inject EVERY formal solution of the scalar annihilators,
    # hence in particular the supplied exact polynomial-Hom basis.
    if matrix.rank()!=columns:
        raise AssertionError('Origin jets fail the proved injectivity check')
    if matrix.row_join(target).rank()==columns:
        if columns:
            coeffs,_=matrix.gauss_jordan_solve(target)
        else:
            coeffs=[]
        polynomials=[]
        for i in range(data['generator_count']):
            length=max((len(b[i]['coefficients']) for b in basis),default=0)
            poly=[sum(coeffs[j]*(sp.Rational(basis[j][i]['coefficients'][n]) if n<len(basis[j][i]['coefficients']) else 0) for j in range(columns)) for n in range(length)]
            polynomials.append(H.serialize_polynomial(poly))
        return dict(outcome='MATCH',selected_polynomial_images=polynomials,
                    origin_certificates=origin,jet_coordinates=keys,
                    actual_jet_match=True,
                    conditional_on_proved_actual_map=True,
                    full_amplitude_coverage_verified=False,
                    gaussian_safety_claimed=False)
    witness=next(v for v in matrix.T.nullspace() if (v.T*target)[0]!=0)
    assert witness.T*matrix==sp.zeros(1,columns)
    return dict(outcome='NO_MATCH',origin_certificates=origin,jet_coordinates=keys,
                polynomial_map_jet_obstruction=[str(c) for c in witness],
                nonzero_original_jet_value=str((witness.T*target)[0]),
                conditional_on_proved_actual_map=True,
                full_amplitude_coverage_verified=False,
                gaussian_safety_claimed=False)


def gaussian_amplitude_jets(phase_terms,amplitude_terms,degrees,covariance):
    """Exact Q-valued Wick coefficients of E[B_i(t,X) exp(tP)].
    Amplitude exponent vectors are (time,spatial_1,...,spatial_n).
    A nondegenerate rational symmetric covariance includes split pairs.
    """
    n=len(covariance)
    C=sp.Matrix([[sp.Rational(c) for c in row] for row in covariance])
    if C.shape!=(n,n) or C!=C.T or C.det()==0:
        raise ValueError('A nondegenerate symmetric rational covariance is required')
    if len(amplitude_terms)!=len(degrees):
        raise ValueError('Amplitude count differs from degree list')
    def parse(terms,width):
        out={}
        for ex,c in terms:
            if len(ex)!=width or any(not isinstance(e,int) or e<0 for e in ex):
                raise ValueError('Bad exponent vector')
            key=tuple(ex);out[key]=out.get(key,sp.S.Zero)+sp.Rational(c)
        return {k:v for k,v in out.items() if v}
    phase=parse(phase_terms,n)
    amps=[parse(a,n+1) for a in amplitude_terms]
    @lru_cache(None)
    def wick(alpha):
        if not any(alpha):return sp.S.One
        if sum(alpha)%2:return sp.S.Zero
        i=next(i for i,a in enumerate(alpha) if a)
        beta=list(alpha);beta[i]-=1
        out=sp.S.Zero
        for j,a in enumerate(beta):
            if a and C[i,j]:
                gamma=beta[:];gamma[j]-=1
                out+=C[i,j]*a*wick(tuple(gamma))
        return out
    def multiply(left,right):
        out={}
        for a,c in left.items():
            for b,d in right.items():
                key=tuple(x+y for x,y in zip(a,b));out[key]=out.get(key,sp.S.Zero)+c*d
        return {k:v for k,v in out.items() if v}
    max_degree=max(degrees,default=-1)
    powers=[{(0,)*n:sp.S.One}]
    for m in range(max_degree):powers.append(multiply(powers[-1],phase))
    output=[]
    for amp,degree in zip(amps,degrees):
        values=[]
        for m in range(degree+1):
            value=sp.S.Zero
            for ex,c in amp.items():
                j=ex[0]
                if j>m:continue
                for alpha,d in powers[m-j].items():
                    value+=c*d*wick(tuple(a+b for a,b in zip(ex[1:],alpha)))/factorial(m-j)
            values.append(str(value))
        output.append(values)
    return output
