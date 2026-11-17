"""Polynomial integration-by-parts certificates for one original marker.

Capped enumeration returns UNKNOWN; uncapped search terminates because the
exact Gaussian source is holonomic. All returned identities are verified
by literal commutative polynomial arithmetic in the ORIGINAL phase.
"""
from itertools import product
import importlib.util
from pathlib import Path
import sympy as sp

_spec=importlib.util.spec_from_file_location('selected',Path(__file__).with_name('gaussian_selected_jets_sidecheck_20261116.py'))
S=importlib.util.module_from_spec(_spec);_spec.loader.exec_module(S)
H=S.H


def exponent_vectors(width,degree):
    for ex in product(range(degree+1),repeat=width):
        if sum(ex)<=degree:yield ex


def find(phase_terms,marker_exponent,covariance,max_degree=5,max_columns=500):
    n=len(marker_exponent)
    if n<1 or any(not isinstance(e,int) or e<0 for e in marker_exponent):
        raise ValueError('A nonnegative original spatial marker exponent is required')
    C=sp.Matrix([[sp.Rational(c) for c in row] for row in covariance])
    if C.shape!=(n,n) or C!=C.T or C.det()==0:raise ValueError('Invalid covariance')
    Cinv=C.inv()
    t=sp.Symbol('t');xs=sp.symbols('x0:'+str(n));variables=(t,*xs)
    phase=sp.S.Zero
    for ex,c in phase_terms:
        if len(ex)!=n or any(not isinstance(e,int) or e<0 for e in ex):raise ValueError('Invalid phase')
        phase+=sp.Rational(c)*sp.prod(x**e for x,e in zip(xs,ex))
    marker=sp.prod(x**e for x,e in zip(xs,marker_exponent))
    drift=[-sum(Cinv[i,j]*xs[j] for j in range(n))+t*sp.diff(phase,xs[i]) for i in range(n)]
    d=0
    while max_degree is None or d<=max_degree:
        opmons=list(exponent_vectors(2,d))
        amp_mons=list(exponent_vectors(n+1,d))
        count=len(opmons)+n*len(amp_mons)
        if max_columns is not None and count>max_columns:
            return dict(outcome='UNKNOWN',reason='IBP column budget exceeded',searched_degree=d)
        columns=[]
        for a,b in opmons:
            columns.append(sp.Poly(-t**a*phase**b*marker,*variables,domain=sp.QQ).as_dict())
        for i in range(n):
            for ex in amp_mons:
                B=sp.prod(x**e for x,e in zip(variables,ex))
                columns.append(sp.Poly(sp.diff(B,xs[i])+drift[i]*B,*variables,domain=sp.QQ).as_dict())
        equations=sorted(set().union(*(col.keys() for col in columns)))
        entries={(i,j):col[key] for j,col in enumerate(columns)
                 for i,key in enumerate(equations) if key in col}
        mat=sp.MutableSparseMatrix(len(equations),count,entries)
        selected=next((v for v in mat.nullspace() if any(v[j] for j in range(len(opmons)))),None)
        if selected is not None:
            lead=next(selected[j] for j in range(len(opmons)) if selected[j])
            selected=[v/lead for v in selected]
            A={ex:selected[j] for j,ex in enumerate(opmons) if selected[j]}
            witnesses=[]
            for i in range(n):
                start=len(opmons)+i*len(amp_mons)
                witnesses.append(sum(selected[start+j]*sp.prod(x**e for x,e in zip(variables,ex))
                                     for j,ex in enumerate(amp_mons)))
            lhs=sum(c*t**a*phase**b*marker for (a,b),c in A.items())
            rhs=sum(sp.diff(B,xs[i])+drift[i]*B for i,B in enumerate(witnesses))
            if sp.expand(lhs-rhs)!=0:raise AssertionError('False original polynomial IBP identity')
            operator=H.op_terms(A)
            scalar=dict(generator_count=1,relations=[[operator]],coordinate_annihilators=[dict(operator=operator,left_witnesses=[[[0,0,'1']]])])
            hom=H.solve(scalar)
            origin=S.origin_certificate(A)
            amp=[[[0,*marker_exponent],'1']]
            jets=S.gaussian_amplitude_jets(phase_terms,[amp],[origin['jet_degree']],covariance)
            match=S.match_actual_jets(scalar,jets,hom)
            return dict(outcome='COMPLETE',original_marker_exponent=marker_exponent,
                        scalar_annihilator=operator,
                        polynomial_ibp_witnesses=[[[list(ex),str(c)] for ex,c in sorted(sp.Poly(B,*variables,domain=sp.QQ).as_dict().items())] for B in witnesses],
                        original_ibp_identity_verified=True,
                        literal_polynomial_identity=str(sp.expand(lhs-rhs)),
                        actual_period_test=match,search_degree=d,
                        all_amplitude_coverage_claimed=False,
                        gaussian_safety_claimed=False)
        d+=1
    return dict(outcome='UNKNOWN',reason='IBP degree budget exceeded',searched_degree=max_degree)
